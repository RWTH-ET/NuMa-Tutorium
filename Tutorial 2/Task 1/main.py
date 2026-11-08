from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np


ORDNER = Path(__file__).resolve().parent
ZIEL = 0.1  # Maximaler Umrechnungsfehler in Grad Celsius.


def vandermonde(x, y):
    """Berechne die Monomkoeffizienten und die Kondition in der Maximumsnorm."""
    # Spalten: 1, x, x**2, ..., x**n. Die Spannungen bleiben in mV.
    matrix = np.vander(x, increasing=True)
    koeffizienten = np.linalg.solve(matrix, y)
    kondition = np.linalg.cond(matrix, np.inf)
    return koeffizienten, kondition


def horner(koeffizienten, stellen):
    """Werte a_0 + a_1*x + ... + a_n*x**n mit dem Horner-Schema aus."""
    wert = 0.0
    for a in reversed(koeffizienten):
        wert = wert * stellen + a
    return wert


def baryzentrische_gewichte(x):
    """Berechne w_j = 1 / Produkt(x_j - x_k) fuer k ungleich j."""
    gewichte = np.ones(len(x))
    for j in range(len(x)):
        for k in range(len(x)):
            if j != k:
                gewichte[j] /= x[j] - x[k]

    # Ein gemeinsamer Faktor kuerzt sich in der baryzentrischen Formel heraus.
    return gewichte / np.max(np.abs(gewichte))


def baryzentrisch(x, y, gewichte, stellen):
    """Werte das Interpolationspolynom an allen gegebenen Stellen aus."""
    werte = np.zeros(len(stellen))
    for i, stelle in enumerate(stellen):
        abstaende = stelle - x
        treffer = np.flatnonzero(abstaende == 0)
        if len(treffer) > 0:
            # An einer Stuetzstelle gilt p(x_j) = y_j; nicht durch null teilen.
            werte[i] = y[treffer[0]]
        else:
            brueche = gewichte / abstaende
            werte[i] = np.sum(brueche * y) / np.sum(brueche)
    return werte


def plotte_ergebnisse(grade, fehler_v, fehler_b, konditionen, raster,
                     abweichung_v, abweichung_b):
    """Zeige und speichere die drei Diagramme aus Aufgabe 1a und 1b."""
    plt.figure(figsize=(9, 5))
    plt.semilogy(grade, fehler_v, "o-", label="Vandermonde + Horner")
    plt.semilogy(grade, fehler_b, "s-", label="Baryzentrisch")
    plt.axhline(ZIEL, color="black", linestyle="--", label="Ziel: 0,1 °C")
    plt.xlabel("Polynomgrad n")
    plt.ylabel("Maximaler absoluter Fehler [°C]")
    plt.title("Aufgabe 1: Äquidistante Stützstellen, Kalibratorwerte")
    plt.grid(True, which="both", alpha=0.3)
    plt.legend()
    plt.tight_layout()
    plt.savefig(ORDNER / "fehlervergleich.png", dpi=160)

    plt.figure(figsize=(9, 5))
    plt.semilogy(grade, konditionen, "o-")
    plt.xlabel("Polynomgrad n")
    plt.ylabel(r"Konditionszahl $\kappa_\infty(V)$")
    plt.title("Aufgabe 1a: Kondition der Vandermonde-Matrix")
    plt.grid(True, which="both", alpha=0.3)
    plt.tight_layout()
    plt.savefig(ORDNER / "kondition.png", dpi=160)

    plt.figure(figsize=(9, 5))
    plt.plot(raster, abweichung_v, label="Vandermonde + Horner")
    plt.plot(raster, abweichung_b, label="Baryzentrisch")
    # Symlog zeigt auch negative Abweichungen; nahe null ist die Achse linear.
    plt.yscale("symlog", linthresh=0.005)
    plt.axhline(0, color="black", linewidth=0.8)
    plt.axhline(ZIEL, color="gray", linestyle="--", label="±0,1 °C")
    plt.axhline(-ZIEL, color="gray", linestyle="--")
    plt.xlabel("Thermospannung E [mV]")
    plt.ylabel("Abweichung p(E) - T_exakt(E) [°C] (symlog)")
    plt.title("Aufgabe 1b: Abweichungen bei n = 20")
    plt.grid(True, which="both", alpha=0.3)
    plt.legend()
    plt.tight_layout()
    plt.savefig(ORDNER / "abweichung_n20.png", dpi=160)
    plt.show()


def main():
    # Die CSV-Dateien enthalten Leerzeichen als Trenner und eine #-Kopfzeile.
    pruefraster = np.loadtxt(ORDNER / "pruefraster.csv")
    stuetzwerte = np.loadtxt(ORDNER / "stuetzwerte_aequidistant.csv")
    raster = pruefraster[:, 0]
    exakt = pruefraster[:, 1]

    grade = list(range(2, 41))
    fehler_v = []
    fehler_b = []
    konditionen = []

    print(f"{'n':>3} {'Fehler Vandermonde [C]':>24} "
          f"{'Fehler baryzentrisch [C]':>26} {'Kondition (inf)':>17}")

    for n in grade:
        daten = stuetzwerte[stuetzwerte[:, 0] == n]
        x = daten[:, 2]  # E_mV
        y = daten[:, 4]  # T_kalibrator_C, nicht T_exakt_C!
        if len(x) != n + 1 or len(np.unique(x)) != n + 1:
            raise ValueError(f"Fuer n = {n} werden {n + 1} verschiedene Stuetzstellen benoetigt.")

        try:
            koeffizienten, kondition = vandermonde(x, y)
            p_v = horner(koeffizienten, raster)
        except np.linalg.LinAlgError:
            # Bei numerischer Singularitaet gibt es kein verlaessliches Ergebnis.
            print(f"n = {n}: Vandermonde-System numerisch singulaer.")
            p_v = np.full_like(raster, np.nan)
            kondition = np.inf

        gewichte = baryzentrische_gewichte(x)
        p_b = baryzentrisch(x, y, gewichte, raster)

        fehler_v.append(np.max(np.abs(p_v - exakt)))
        fehler_b.append(np.max(np.abs(p_b - exakt)))
        konditionen.append(kondition)
        print(f"{n:3d} {fehler_v[-1]:24.6e} {fehler_b[-1]:26.6e} {kondition:17.6e}")

        if n == 20:
            abweichung_v = p_v - exakt
            abweichung_b = p_b - exakt

    for name, fehler, abweichung in [
        ("Vandermonde", fehler_v, abweichung_v),
        ("Baryzentrisch", fehler_b, abweichung_b),
    ]:
        passende_grade = [n for n, e in zip(grade, fehler) if e < ZIEL]
        print(f"\n{name}: Grade mit Fehler < {ZIEL} C: {passende_grade}")
        if np.all(np.isfinite(abweichung)):
            index = np.argmax(np.abs(abweichung))
            print(f"Bei n = 20: groesster Fehler bei E = {raster[index]:.6f} mV, "
                  f"Abweichung = {abweichung[index]:.6f} C")

    plotte_ergebnisse(grade, fehler_v, fehler_b, konditionen, raster,
                     abweichung_v, abweichung_b)


if __name__ == "__main__":
    main()
