from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np


ORDNER = Path(__file__).resolve().parent
RUNDUNGSFEHLER = 0.005  # Halbe Anzeigestelle des Kalibrators in Grad Celsius.


def baryzentrische_gewichte(E_stuetz):
    """Berechne die Gewichte einmal fuer die gegebenen Stuetzstellen."""
    w = np.ones(len(E_stuetz))
    for j in range(len(E_stuetz)):
        produkt = 1.0
        for k in range(len(E_stuetz)):
            if j != k:
                produkt = produkt * (E_stuetz[j] - E_stuetz[k])
        w[j] = 1.0 / produkt

    # Ein gemeinsamer Faktor kuerzt sich in der baryzentrischen Formel heraus.
    w = w / np.max(np.abs(w))
    return w


def baryzentrische_basis(E_stuetz, w, E):
    """Liefere die Lagrange-Basiswerte fuer eine einzelne Spannung E."""
    differenzen = E - E_stuetz
    treffer = np.flatnonzero(differenzen == 0)
    if len(treffer) > 0:
        # An einer Stuetzstelle ist deren Basiswert 1, alle anderen sind 0.
        basis = np.zeros(len(E_stuetz))
        basis[treffer[0]] = 1.0
        return basis

    terme = w / differenzen
    basis = terme / np.sum(terme)
    return basis


def berechne_fehler(E_stuetz, T_exakt, T_kalibrator, E_raster, T_raster):
    """Vergleiche beide Polynome und berechne die Schranke auf dem Pruefraster."""
    w = baryzentrische_gewichte(E_stuetz)
    p_exakt = np.zeros(len(E_raster))
    p_kalibrator = np.zeros(len(E_raster))
    lebesgue = np.zeros(len(E_raster))

    for i, E in enumerate(E_raster):
        basis = baryzentrische_basis(E_stuetz, w, E)
        # Gleiche Stuetzstellen: Die Basis wird fuer beide Polynome verwendet.
        p_exakt[i] = np.sum(basis * T_exakt)
        p_kalibrator[i] = np.sum(basis * T_kalibrator)
        lebesgue[i] = np.sum(np.abs(basis))

    e_interp = np.max(np.abs(p_exakt - T_raster))
    e_daten = np.max(np.abs(p_kalibrator - p_exakt))
    lambda_raster = np.max(lebesgue)
    schranke = lambda_raster * RUNDUNGSFEHLER
    return e_interp, e_daten, lambda_raster, schranke


def untersuche_stuetzstellen(art):
    """Berechne die Fehler fuer n = 2 bis 40 und eine Stuetzstellenwahl."""
    E_raster, T_raster = np.loadtxt(ORDNER / "pruefraster.csv").T
    stuetzwerte = np.loadtxt(ORDNER / f"stuetzwerte_{art}.csv")
    ergebnisse = []

    print(f"\n{art}")
    print(f"{'n':>3} {'Interpolation [C]':>19} {'Datenfehler [C]':>19} "
          f"{'Lambda (Raster)':>18} {'Schranke [C]':>19}")

    for n in range(2, 41):
        daten_n = stuetzwerte[stuetzwerte[:, 0] == n]
        E_stuetz = daten_n[:, 2]
        T_exakt = daten_n[:, 3]
        T_kalibrator = daten_n[:, 4]

        if len(E_stuetz) != n + 1 or len(np.unique(E_stuetz)) != n + 1:
            raise ValueError(f"Fuer n = {n} werden {n + 1} verschiedene Stuetzstellen benoetigt.")
        if np.max(np.abs(T_kalibrator - T_exakt)) > RUNDUNGSFEHLER + 1e-10:
            raise ValueError(f"Die Kalibratorwerte fuer n = {n} verletzen die Rundungsgrenze.")

        e_interp, e_daten, lam, schranke = berechne_fehler(
            E_stuetz, T_exakt, T_kalibrator, E_raster, T_raster)
        ergebnisse.append([n, e_interp, e_daten, lam, schranke])
        print(f"{n:3d} {e_interp:19.6e} {e_daten:19.6e} "
              f"{lam:18.6e} {schranke:19.6e}")

    return np.array(ergebnisse)


def plotte_ergebnisse(art, name, ergebnisse):
    """Speichere die drei geforderten Kurven in einem gemeinsamen Diagramm."""
    n_werte = ergebnisse[:, 0]
    plt.figure(figsize=(9, 5))
    plt.semilogy(n_werte, ergebnisse[:, 1], "o-", label="Interpolationsfehler")
    plt.semilogy(n_werte, ergebnisse[:, 2], "s-", label="Datenfehler")
    plt.semilogy(n_werte, ergebnisse[:, 4], "--", label="Schranke: Λ (Prüfraster) · 0,005 °C")
    plt.xlabel("Polynomgrad n")
    plt.ylabel("Maximaler absoluter Fehler / Schranke [°C] (log)")
    plt.title(f"Aufgabe 3b: {name}")
    plt.grid(True, which="both", alpha=0.3)
    plt.legend()
    plt.tight_layout()
    plt.savefig(ORDNER / f"3b_{art}.png", dpi=160)


def main():
    for art, name in [("aequidistant", "Äquidistante Stützstellen"),
                      ("tschebyscheff", "Tschebyscheff-Stützstellen")]:
        ergebnisse = untersuche_stuetzstellen(art)
        np.savetxt(
            ORDNER / f"ergebnisse_{art}.csv", ergebnisse, delimiter=",",
            header="n,e_interp_C,e_daten_C,Lambda_pruefraster,Schranke_C",
            comments="", fmt=["%d", "%.12e", "%.12e", "%.12e", "%.12e"])
        plotte_ergebnisse(art, name, ergebnisse)

    print(f"\nDiagramme und Ergebnistabellen gespeichert in: {ORDNER}")
    plt.show()


if __name__ == "__main__":
    main()
