import importlib.util
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np


ORDNER = Path(__file__).resolve().parent
RUNDUNGSFEHLER = 0.005  # Halbe Anzeigestelle des Kalibrators in °C.


# Die Interpolationsverfahren direkt aus Task 1 wiederverwenden.
# Der Dateipfad umgeht die Leerzeichen im Ordnernamen; main() wird nicht gestartet.
spec = importlib.util.spec_from_file_location(
    "task1_interpolation", ORDNER.parent / "Task 1" / "main.py")
task1 = importlib.util.module_from_spec(spec)
spec.loader.exec_module(task1)


def berechne_fehler(x, y_exakt, y_kalibrator, E_raster, T_raster):
    """Interpolationsfehler, Datenfehler und Schranke auf dem Prüfraster."""
    gewichte = task1.baryzentrische_gewichte(x)
    p_exakt = []
    p_kalibrator = []
    lebesgue_werte = []

    for E in E_raster:
        p_exakt.append(task1.baryzentrisch(x, y_exakt, gewichte, E))
        p_kalibrator.append(task1.baryzentrisch(x, y_kalibrator, gewichte, E))

        # Nur die Lebesgue-Konstante benötigt zusätzlich die einzelnen Basiswerte.
        if E in x:
            # An einer Stützstelle gilt direkt der zugehörige Stützwert.
            basis = (x == E).astype(float)
        else:
            anteile = gewichte / (E - x)
            basis = anteile / np.sum(anteile)

        lebesgue_werte.append(np.sum(np.abs(basis)))

    p_exakt = np.array(p_exakt)
    p_kalibrator = np.array(p_kalibrator)
    e_interp = np.max(np.abs(p_exakt - T_raster))
    e_daten = np.max(np.abs(p_kalibrator - p_exakt))
    lambda_raster = np.max(lebesgue_werte)
    return e_interp, e_daten, lambda_raster, lambda_raster * RUNDUNGSFEHLER


def plotte_ergebnisse(art, titel, ergebnisse):
    n = ergebnisse[:, 0]
    plt.figure(figsize=(9, 5))
    plt.semilogy(n, ergebnisse[:, 1], "o-", label="Interpolationsfehler")
    plt.semilogy(n, ergebnisse[:, 2], "s-", label="Datenfehler")
    plt.semilogy(n, ergebnisse[:, 4], "--", label="Schranke: Λ (Prüfraster) · 0,005 °C")
    plt.xlabel("Polynomgrad n")
    plt.ylabel("Maximaler absoluter Fehler / Schranke [°C] (log)")
    plt.title(f"Aufgabe 3b: {titel}")
    plt.grid(True, which="both", alpha=0.3)
    plt.legend()
    plt.tight_layout()
    plt.savefig(ORDNER / f"3b_{art}.png", dpi=160)


def main():
    E_raster, T_raster = np.loadtxt(ORDNER / "pruefraster.csv").T

    for art, titel in [("aequidistant", "Äquidistante Stützstellen"),
                       ("tschebyscheff", "Tschebyscheff-Stützstellen")]:
        stuetzwerte = np.loadtxt(ORDNER / f"stuetzwerte_{art}.csv")
        ergebnisse = []

        for n in range(2, 41):
            daten = stuetzwerte[stuetzwerte[:, 0] == n]
            x = daten[:, 2]  # Spannungen der Stützstellen.
            y_exakt = daten[:, 3]
            y_kalibrator = daten[:, 4]

            fehler = berechne_fehler(x, y_exakt, y_kalibrator, E_raster, T_raster)
            ergebnisse.append([n, *fehler])

        ergebnisse = np.array(ergebnisse)
        np.savetxt(
            ORDNER / f"ergebnisse_{art}.csv", ergebnisse, delimiter=",",
            header="n,e_interp_C,e_daten_C,Lambda_pruefraster,Schranke_C",
            comments="", fmt=["%d", "%.12e", "%.12e", "%.12e", "%.12e"])
        plotte_ergebnisse(art, titel, ergebnisse)

    print(f"Diagramme und Ergebnistabellen gespeichert in: {ORDNER}")
    plt.show()


if __name__ == "__main__":
    main()
