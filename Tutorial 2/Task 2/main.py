import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

# Funktionen für die Interpolationsverfahren

def vandermonde_horner(x_nodes, y_nodes, x_eval):
    """Aufgabe 1a: Vandermonde-Interpolation mit Horner-Schema"""
    # Vandermonde-Matrix V(X) erstellen (Mächte aufsteigend)
    V = np.vander(x_nodes, increasing=True)
    
    # Koeffizienten a durch Lösen des LGS bestimmen
    a = np.linalg.solve(V, y_nodes)
    
    # Horner-Schema zur Auswertung von p(x)
    p_eval = np.zeros_like(x_eval)
    for coeff in a[::-1]: # Rückwärts durch die Koeffizienten iterieren
        p_eval = p_eval * x_eval + coeff
        
    return p_eval

def baryzentrisch(x_nodes, y_nodes, x_eval):
    """Aufgabe 1b: Baryzentrische Interpolation"""
    n = len(x_nodes)
    w = np.ones(n)
    
    # Einmalige Berechnung der Gewichte w_j
    for j in range(n):
        for k in range(n):
            if j != k:
                w[j] /= (x_nodes[j] - x_nodes[k])
                
    # Auswertung von p(x)
    p_eval = np.zeros_like(x_eval, dtype=float)
    for i, x in enumerate(x_eval):
        # Gemäß Aufgabenstellung ist x != x_j im Prüfraster garantiert
        diff = x - x_nodes
        term = w / diff
        p_eval[i] = np.sum(term * y_nodes) / np.sum(term)
        
    return p_eval

# Daten laden und vorbereiten

# Prüfraster laden (erste Zeile überspringen, Namen manuell vergeben)
raster = pd.read_csv('Tutorial 2/Task 2/pruefraster.csv', sep=r'\s+', skiprows=1, names=['E_mV', 'T_C'])
E_test = raster['E_mV'].values
T_exakt = raster['T_C'].values

# Tschebyscheff-Stützstellen laden (Spaltennamen anpassen, falls sie in der Datei leicht anders heißen)
stuetz = pd.read_csv('Tutorial 2/Task 2/stuetzwerte_tschebyscheff.csv', sep=r'\s+', skiprows=1, names=['n', 'i', 'E_mV', 'T_exakt_C', 'T_kalibrator_C'])


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
    p_vander = vandermonde_horner(E_nodes, T_kalib, E_test)
    e_max_vander = np.max(np.abs(p_vander - T_exakt))
    fehler_vander.append(e_max_vander)
    
    # Verfahren 1b anwenden
    p_bary = baryzentrisch(E_nodes, T_kalib, E_test)
    e_max_bary = np.max(np.abs(p_bary - T_exakt))
    fehler_bary.append(e_max_bary)

# Plotten der Ergebnisse

plt.figure(figsize=(10, 6))
plt.semilogy(n_werte, fehler_vander, marker='o', label='Vandermonde + Horner (1a)')
plt.semilogy(n_werte, fehler_bary, marker='s', label='Baryzentrisch (1b)')

plt.axhline(y=0.1, color='black', linestyle='--', label='Ziel: 0,1 °C')

plt.title('Aufgabe 2: Fehlervergleich mit Tschebyscheff-Stützstellen')
plt.xlabel('Polynomgrad n')
plt.ylabel('Maximaler absoluter Fehler [°C] (log)')
plt.grid(True, which="both", ls="--", alpha=0.5)
plt.legend()
plt.tight_layout()

plt.savefig('Tutorial 2/Task 2/aufgabe2_tschebyscheff.png', dpi=300, bbox_inches='tight')

plt.show()