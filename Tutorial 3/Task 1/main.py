"""
Programmierprojekt 3 – Numerische Mathematik für Elektrotechnik (WS 2026/27)
Teil A (Aufgaben 1–3): Quadraturformeln implementieren und testen
"""

from pathlib import Path
import numpy as np
import matplotlib.pyplot as plt

ORDNER = Path(__file__).resolve().parent
TASK_1 = ORDNER.parent / "Task 1"


# =============================================================================
# 1. Schnittstellen der Quadraturformeln (vgl. Folie 3)
# =============================================================================

def _berechne_gauss_formel(N: int) -> tuple[np.ndarray, np.ndarray]:
    """
    Berechnet die N-Punkt-Gauß-Legendre-Quadraturformel auf [-1, 1]
    über das symmetrische tridiagonale Eigenwertproblem (Prop. 3.24).
    """
    if N < 1:
        raise ValueError("N muss mindestens 1 sein.")
    if N == 1:
        return np.array([0.0]), np.array([2.0])

    # Nebendiagonalelemente: beta_k = (4 - k^(-2))^(-1/2) für k = 1, ..., N-1
    k = np.arange(1, N, dtype=float)
    beta = (4.0 - k ** (-2)) ** (-0.5)

    # Symmetrische Tridiagonalmatrix J mit Nullen auf der Hauptdiagonale
    J = np.diag(beta, k=-1) + np.diag(beta, k=1)

    # Eigenwerte = Stützstellen, 2 * (1. Komponente der Eigenvektoren)^2 = Gewichte
    eigvals, eigvecs = np.linalg.eigh(J)
    x_hat = eigvals
    w_hat = 2.0 * (eigvecs[0, :] ** 2)
    return x_hat, w_hat


# Vorabberechnung für N = 1, ..., 48 (Aufgabe 1: 1..40, Aufgabe 3 benötigt bis 48)
GAUSS_CACHE: dict[int, tuple[np.ndarray, np.ndarray]] = {
    N: _berechne_gauss_formel(N) for N in range(1, 49)
}


def gauss_formel(N: int) -> tuple[np.ndarray, np.ndarray]:
    """Gibt die vorab berechnete Gauß-Legendre-Formel (x_hat, w_hat) auf [-1, 1] zurück."""
    if N not in GAUSS_CACHE:
        GAUSS_CACHE[N] = _berechne_gauss_formel(N)
    x_hat, w_hat = GAUSS_CACHE[N]
    return x_hat.copy(), w_hat.copy()


def _lagrange_baryzentrisch(x_nodes: np.ndarray, x_eval: np.ndarray) -> np.ndarray:
    """
    Wertet alle Lagrange-Basispolynome l_i(x) zu den Stützstellen x_nodes
    an den Punkten x_eval mittels der zweiten baryzentrischen Formel (PP2 / Satz 3.5) aus.
    Rückgabe: Matrix L mit L[i, j] = l_i(x_eval[j]).
    """
    N = len(x_nodes)
    M_eval = len(x_eval)

    # Baryzentrische Gewichte: lambda_i = 1 / prod_{j != i} (x_i - x_j)
    # Skalierung verhindert Über-/Unterlauf für großes N
    scale = 0.5 * (N - 1) if N > 1 else 1.0
    x_scaled = x_nodes * scale
    bary_w = np.empty(N, dtype=float)
    for i in range(N):
        diffs = x_scaled[i] - np.delete(x_scaled, i)
        bary_w[i] = 1.0 / np.prod(diffs)

    L = np.zeros((N, M_eval), dtype=float)
    for j, xq in enumerate(x_eval):
        diff = xq - x_nodes
        hit = np.where(np.abs(diff) < 1e-15)[0]
        if len(hit) > 0:
            L[hit[0], j] = 1.0
        else:
            terms = bary_w / diff
            L[:, j] = terms / np.sum(terms)

    return L


def nc_formel(N: int) -> tuple[np.ndarray, np.ndarray]:
    """
    Berechnet die abgeschlossene N-Punkt-Newton-Cotes-Formel (x_hat, w_hat) auf [-1, 1].
    Gewichte w_hat_i = int_{-1}^1 l_i(x) dx über baryzentrische Auswertung
    und exakte Gauß-Legendre-Quadratur.
    """
    if N < 1:
        raise ValueError("N muss mindestens 1 sein.")
    if N == 1:
        return np.array([0.0]), np.array([2.0])

    x_hat = np.linspace(-1.0, 1.0, N)

    # l_i hat Grad N - 1 -> Gauß-Formel mit ceil(N / 2) Punkten integriert l_i exakt
    N_gauss = max(1, (N + 1) // 2)
    x_g, w_g = gauss_formel(N_gauss)

    L = _lagrange_baryzentrisch(x_hat, x_g)
    w_hat = L @ w_g
    return x_hat, w_hat


def transformiere(
    x_hat: np.ndarray, w_hat: np.ndarray, a: float, b: float
) -> tuple[np.ndarray, np.ndarray]:
    """Überträgt (x_hat, w_hat) von [-1, 1] auf [a, b] (Lemma 3.3)."""
    mitte = 0.5 * (a + b)
    halbbreite = 0.5 * (b - a)
    x = mitte + halbbreite * x_hat
    w = halbbreite * w_hat
    return x, w


def summierte_formel(
    x_hat: np.ndarray, w_hat: np.ndarray, a: float, b: float, M: int
) -> tuple[np.ndarray, np.ndarray]:
    """
    Erzeugt die summierte Quadraturformel auf [a, b] mit M Teilintervallen.
    Gemeinsame Randpunkte (bei abgeschlossenen Formeln) werden nur einmal gespeichert
    und ihre Gewichte addiert.
    """
    if M < 1:
        raise ValueError("Anzahl der Teilintervalle M muss mindestens 1 sein.")

    N = len(x_hat)
    intervallgrenzen = np.linspace(a, b, M + 1)

    abgeschlossen = (
        N >= 2
        and np.isclose(x_hat[0], -1.0, atol=1e-14)
        and np.isclose(x_hat[-1], 1.0, atol=1e-14)
    )

    if abgeschlossen:
        gesamt_punkte = M * (N - 1) + 1
        x_ges = np.empty(gesamt_punkte, dtype=float)
        w_ges = np.zeros(gesamt_punkte, dtype=float)

        for j in range(M):
            x_sub, w_sub = transformiere(
                x_hat, w_hat, intervallgrenzen[j], intervallgrenzen[j + 1]
            )
            start = j * (N - 1)
            if j == 0:
                x_ges[0:N] = x_sub
                w_ges[0:N] = w_sub
            else:
                w_ges[start] += w_sub[0]
                x_ges[start + 1 : start + N] = x_sub[1:]
                w_ges[start + 1 : start + N] = w_sub[1:]
        return x_ges, w_ges

    # Offene Formel (z. B. Gauß-Legendre)
    x_liste, w_liste = [], []
    for j in range(M):
        x_sub, w_sub = transformiere(
            x_hat, w_hat, intervallgrenzen[j], intervallgrenzen[j + 1]
        )
        x_liste.append(x_sub)
        w_liste.append(w_sub)
    return np.concatenate(x_liste), np.concatenate(w_liste)


def quadratur(y: np.ndarray, w: np.ndarray) -> float:
    """Berechnet die Quadratursumme Q = sum_i w_i * y_i."""
    return float(np.dot(w, y))


def exaktheitsgrad(x_hat: np.ndarray, w_hat: np.ndarray, tol: float = 1e-12) -> int:
    """
    Bestimmt den gemessenen Exaktheitsgrad m durch Testen der Monome x^k (k = 0, 1, ...),
    bis |Q(x^k) - int_{-1}^1 x^k dx| > tol.
    """
    k = 0
    while k < 500:
        exakt = 2.0 / (k + 1) if (k % 2 == 0) else 0.0
        q_val = quadratur(x_hat ** k, w_hat)
        if abs(q_val - exakt) > tol:
            return k - 1
        k += 1
    return k - 1


# =============================================================================
# 2. Auswertung & Plots für Aufgaben 1, 2 und 3
# =============================================================================

def main() -> None:
    # Zielordner "Task 1" (parallel zum Skriptordner) sicherstellen
    TASK_1.mkdir(parents=True, exist_ok=True)
    tol = 1e-12

    # -------------------------------------------------------------------------
    # Aufgabe 1: Gauß-Legendre (N = 1, ..., 40)
    # -------------------------------------------------------------------------
    N_gauss_werte = np.arange(1, 41)
    exakt_gauss = np.array(
        [exaktheitsgrad(*gauss_formel(int(N)), tol=tol) for N in N_gauss_werte]
    )
    kond_gauss = np.array(
        [np.sum(np.abs(gauss_formel(int(N))[1])) for N in N_gauss_werte]
    )

    fig1, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 4.8))
    ax1.plot(N_gauss_werte, exakt_gauss, "o-", label="Gemessen (Gauß-Legendre)")
    ax1.plot(
        N_gauss_werte, 2 * N_gauss_werte - 1, "--", color="black", label="Theorie ($2N - 1$)"
    )
    ax1.set_title("Aufgabe 1: Exaktheitsgrad (Gauß-Legendre)")
    ax1.set_xlabel("Anzahl Stützstellen $N$")
    ax1.set_ylabel("Exaktheitsgrad $m$")
    ax1.grid(True, which="both", linestyle=":", alpha=0.7)
    ax1.legend()

    ax2.semilogy(N_gauss_werte, kond_gauss, "o-", label=r"Gauß-Legendre: $\sum_i |\hat{w}_i|$")
    ax2.set_title(r"Aufgabe 1: Konditionszahl $\sum_i |\hat{w}_i|$")
    ax2.set_xlabel("Anzahl Stützstellen $N$")
    ax2.set_ylabel(r"$\sum_i |\hat{w}_i|$ (log-Skala)")
    ax2.set_ylim(1, 10)
    ax2.grid(True, which="both", linestyle=":", alpha=0.7)
    ax2.legend()

    fig1.tight_layout()
    fig1.savefig(TASK_1 / "aufgabe1_gauss_legendre.png", dpi=300, bbox_inches="tight")
    plt.show()

    # -------------------------------------------------------------------------
    # Aufgabe 2: Newton-Cotes ergänzen (N = 2, ..., 40)
    # -------------------------------------------------------------------------
    N_nc_werte = np.arange(2, 41)
    exakt_nc = np.array(
        [exaktheitsgrad(*nc_formel(int(N)), tol=tol) for N in N_nc_werte]
    )
    kond_nc = np.array(
        [np.sum(np.abs(nc_formel(int(N))[1])) for N in N_nc_werte]
    )
    theorie_nc = np.where(N_nc_werte % 2 == 0, N_nc_werte - 1, N_nc_werte)

    fig2, (ax1, ax2) = plt.subplots(1, 2, figsize=(12.5, 5.0))
    ax1.plot(N_gauss_werte, exakt_gauss, "o-", markersize=4, label="Gauß-Legendre (gemessen)")
    ax1.plot(
        N_gauss_werte, 2 * N_gauss_werte - 1, "--", color="tab:blue", alpha=0.6,
        label="Gauß Theorie ($2N-1$)"
    )
    ax1.plot(N_nc_werte, exakt_nc, "s-", markersize=4, color="tab:orange", label="Newton-Cotes (gemessen)")
    ax1.plot(
        N_nc_werte, theorie_nc, "--", color="tab:orange", alpha=0.6,
        label="NC Theorie ($N-1$ bzw. $N$)"
    )
    ax1.set_title(f"Aufgabe 1 & 2: Exaktheitsgrad (tol = {tol:.0e})")
    ax1.set_xlabel("Anzahl Stützstellen $N$")
    ax1.set_ylabel("Exaktheitsgrad $m$")
    ax1.grid(True, which="both", linestyle=":", alpha=0.7)
    ax1.legend()

    ax2.semilogy(N_gauss_werte, kond_gauss, "o-", markersize=4, label="Gauß-Legendre")
    ax2.semilogy(N_nc_werte, kond_nc, "s-", markersize=4, color="tab:orange", label="Newton-Cotes")
    ax2.axhline(
        tol / np.finfo(float).eps, color="red", linestyle=":",
        label=r"Schwelle $\mathrm{tol} / \varepsilon_{\mathrm{mach}}$"
    )
    ax2.set_title(r"Aufgabe 1 & 2: Konditionszahl $\sum_i |\hat{w}_i|$")
    ax2.set_xlabel("Anzahl Stützstellen $N$")
    ax2.set_ylabel(r"$\sum_i |\hat{w}_i|$ (log-Skala)")
    ax2.grid(True, which="both", linestyle=":", alpha=0.7)
    ax2.legend()

    fig2.tight_layout()
    fig2.savefig(TASK_1 / "aufgabe2_vergleich_gauss_nc.png", dpi=300, bbox_inches="tight")
    plt.show()

    # -------------------------------------------------------------------------
    # Aufgabe 3: Summierte Quadratur auf [0, 1] mit Budget <= 48 Auswertungen
    # -------------------------------------------------------------------------
    budget = 48
    a, b = 0.0, 1.0

    funktionen = [
        (r"$f_1(x) = e^x$", lambda x: np.exp(x), np.e - 1.0),
        (
            r"$f_2(x) = \frac{1}{1 + 400(x - 0{,}3)^2}$",
            lambda x: 1.0 / (1.0 + 400.0 * (x - 0.3) ** 2),
            (np.arctan(14.0) + np.arctan(6.0)) / 20.0,
        ),
        (
            r"$f_3(x) = |\sin(5x)|$",
            lambda x: np.abs(np.sin(5.0 * x)),
            (3.0 + np.cos(5.0)) / 5.0,
        ),
    ]

    # Gauß: M * N <= 48  =>  M = floor(48 / N)
    N_g_list = np.arange(1, budget + 1)
    M_g_list = budget // N_g_list

    # Newton-Cotes (gemeinsame Randpunkte einmal): M*(N - 1) + 1 <= 48  =>  M = floor(47 / (N - 1))
    N_nc_list = np.arange(2, budget + 1)
    M_nc_list = (budget - 1) // (N_nc_list - 1)

    fig3, axes = plt.subplots(1, 3, figsize=(15.5, 4.8))
    eps_floor = 1e-16

    for ax, (titel_tex, f, exakt) in zip(axes, funktionen):
        fehler_g = []
        for N, M in zip(N_g_list, M_g_list):
            x_hat, w_hat = gauss_formel(int(N))
            x, w = summierte_formel(x_hat, w_hat, a, b, int(M))
            Q = quadratur(f(x), w)
            fehler_g.append(max(abs(Q - exakt), eps_floor))

        fehler_nc = []
        for N, M in zip(N_nc_list, M_nc_list):
            x_hat, w_hat = nc_formel(int(N))
            x, w = summierte_formel(x_hat, w_hat, a, b, int(M))
            Q = quadratur(f(x), w)
            fehler_nc.append(max(abs(Q - exakt), eps_floor))

        ax.semilogy(N_g_list, fehler_g, "o-", markersize=3.5, label="Gauß-Legendre")
        ax.semilogy(
            N_nc_list, fehler_nc, "s-", markersize=3.5, color="tab:orange", label="Newton-Cotes"
        )
        ax.set_title(titel_tex)
        ax.set_xlabel("Stützstellen $N$ pro Teilintervall")
        ax.set_ylabel("Absoluter Fehler $|Q - I|$")
        ax.set_ylim(bottom=5e-17)
        ax.grid(True, which="both", linestyle=":", alpha=0.7)
        ax.legend()

    fig3.suptitle(
        "Aufgabe 3: Summierte Quadratur auf $[0, 1]$ (Budget $\\leq 48$ Auswertungen)",
        fontsize=13,
    )
    fig3.tight_layout()
    fig3.savefig(TASK_1 / "aufgabe3_summierte_quadratur.png", dpi=300, bbox_inches="tight")
    plt.show()


if __name__ == "__main__":
    main()