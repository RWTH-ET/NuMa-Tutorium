import importlib.util
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt


ORDNER = Path(__file__).resolve().parent
TASK_1 = ORDNER.parent / "Task 1"

# Der Dateipfad umgeht die Leerzeichen im Ordnernamen.
# Der geschuetzte main()-Aufruf von Task 1 wird beim Import nicht ausgefuehrt.
spec = importlib.util.spec_from_file_location("task1_interpolation", TASK_1 / "main.py")
task1 = importlib.util.module_from_spec(spec)
spec.loader.exec_module(task1)


def main():
    # Daten laden und vorbereiten

    # Prüfraster laden (erste Zeile überspringen, Namen manuell vergeben)
    raster = pd.read_csv(TASK_1 / 'pruefraster.csv', sep=r'\s+', skiprows=1, names=['E_mV', 'T_C'])
    E_test = raster['E_mV'].values
    T_exakt = raster['T_C'].values

    # Tschebyscheff-Stützstellen laden (Spaltennamen anpassen, falls sie in der Datei leicht anders heißen)
    stuetz = pd.read_csv(TASK_1 / 'stuetzwerte_tschebyscheff.csv', sep=r'\s+', skiprows=1, names=['n', 'i', 'E_mV', 'T_exakt_C', 'T_kalibrator_C'])


    # Experimente durchführen (n = 2 bis 40)

    n_werte = range(2, 41)
    fehler_vander = []
    fehler_bary = []

    for n in n_werte:
        # Daten für das aktuelle n filtern
        daten_n = stuetz[stuetz['n'] == n]
        E_nodes = daten_n['E_mV'].values
        T_kalib = daten_n['T_kalibrator_C'].values # Kalibratorwerte als y_i

        # Verfahren 1a anwenden
        try:
            a, _ = task1.vandermonde(E_nodes, T_kalib)
            p_vander = np.array([task1.horner(a, E) for E in E_test])
        except np.linalg.LinAlgError:
            print(f"n = {n}: Vandermonde-System numerisch singulaer.")
            p_vander = np.full_like(E_test, np.nan)
        e_max_vander = np.max(np.abs(p_vander - T_exakt))
        fehler_vander.append(e_max_vander)

        # Verfahren 1b anwenden
        w = task1.baryzentrische_gewichte(E_nodes)
        p_bary = np.array([
            task1.baryzentrisch(E_nodes, T_kalib, w, E) for E in E_test
        ])
        e_max_bary = np.max(np.abs(p_bary - T_exakt))
        fehler_bary.append(e_max_bary)

    # Plotten der Ergebnisse

    plt.figure(figsize=(10, 6))
    plt.semilogy(n_werte, fehler_vander, marker='o', label='Vandermonde + Horner (1a)')
    plt.semilogy(n_werte, fehler_bary, marker='s', label='Baryzentrisch (1b)')

    plt.axhline(y=task1.ZIEL, color='black', linestyle='--', label='Ziel: 0,1 °C')

    plt.title('Aufgabe 2: Fehlervergleich mit Tschebyscheff-Stützstellen')
    plt.xlabel('Polynomgrad n')
    plt.ylabel('Maximaler absoluter Fehler [°C] (log)')
    plt.grid(True, which="both", ls="--", alpha=0.5)
    plt.legend()
    plt.tight_layout()

    plt.savefig(ORDNER / 'aufgabe2_tschebyscheff.png', dpi=300, bbox_inches='tight')

    plt.show()


if __name__ == "__main__":
    main()
