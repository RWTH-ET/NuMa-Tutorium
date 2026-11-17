from pathlib import Path
import numpy as np
import matplotlib.pyplot as plt

# Pfade definieren
ORDNER = Path(__file__).resolve().parent
TASK_2 = ORDNER.parent / "Task 2"
TASK_2.mkdir(parents=True, exist_ok=True)

# 1. WERKZEUGE AUS TEIL A

def gauss_formel(N):
    """Berechnet Stützstellen x_hat und Gewichte w_hat auf [-1, 1] via Eigenwertproblem."""
    if N == 1:
        return np.array([0.0]), np.array([2.0])
    k = np.arange(1, N)
    beta = 1.0 / np.sqrt(4.0 - k**(-2.0))
    J = np.diag(beta, k=1) + np.diag(beta, k=-1)
    eigenvalues, eigenvectors = np.linalg.eigh(J)
    x_hat = eigenvalues
    w_hat = 2.0 * (eigenvectors[0, :] ** 2)
    return x_hat, w_hat

def nc_formel(N):
    """Berechnet äquidistante Stützstellen x_hat und Gewichte w_hat auf [-1, 1]."""
    if N < 2:
        raise ValueError("Newton-Cotes benötigt N >= 2 Stützstellen.")
    x_hat = np.linspace(-1.0, 1.0, N)
    
    x_g, w_g = gauss_formel(max(N + 5, 20))
    w_hat = np.zeros(N)
    for i in range(N):
        def l_i(x):
            prod = np.ones_like(x, dtype=float)
            for j in range(N):
                if j != i:
                    prod *= (x - x_hat[j]) / (x_hat[i] - x_hat[j])
            return prod
        w_hat[i] = np.sum(w_g * l_i(x_g))
    return x_hat, w_hat

def transformiere(x_hat, w_hat, a, b):
    """Transformiert Stützstellen x_hat und Gewichte w_hat von [-1, 1] nach [a, b]."""
    x = 0.5 * (b - a) * x_hat + 0.5 * (a + b)
    w = 0.5 * (b - a) * w_hat
    return x, w

def quadratur(y, w):
    """Berechnet die Quadratur Q = sum(w_i * y_i)."""
    return np.sum(w * y)

# 2. MODELLE FÜR DIE KONDENSATORENERGIE (TEIL B)

def C_exakt(v):
    return 10.0 / (1.0 + (v / 8.0)**2)

def C_gemessen(v):
    return np.round(C_exakt(v), 2)

def f_exakt(v):
    return v * C_exakt(v)

def f_gemessen(v):
    return v * C_gemessen(v)

W_EXAKT = 320.0 * np.log(3.25)

# 3. EVALUIERUNG UND FEHLERTRENNUNG (AUFGABE 4)

def analisiere_verfahren(formel_fn, N_range, a=0.0, b=12.0):
    ergebnisse = []
    for N in N_range:
        x_hat, w_hat = formel_fn(N)
        x, w = transformiere(x_hat, w_hat, a, b)
        
        y_exakt = f_exakt(x)
        y_gemessen = f_gemessen(x)
        
        Q_f = quadratur(y_exakt, w)
        Q_tilde = quadratur(y_gemessen, w)
        
        verfahrensfehler = abs(Q_f - W_EXAKT)
        datenfehler = abs(Q_tilde - Q_f)
        gesamtfehler = abs(Q_tilde - W_EXAKT)
        
        ergebnisse.append({
            'N': N,
            'gesamtfehler': gesamtfehler,
            'verfahrensfehler': verfahrensfehler,
            'datenfehler': datenfehler
        })
    return ergebnisse

N_gauss_range = range(1, 16)
N_nc_range = range(2, 16)

res_gauss = analisiere_verfahren(gauss_formel, N_gauss_range)
res_nc = analisiere_verfahren(nc_formel, N_nc_range)

# 4. PLOTS UND SPEICHERN VIA PATHLIB

plt.figure(figsize=(12, 5))

# Gauß-Legendre Plot
plt.subplot(1, 2, 1)
plt.semilogy([r['N'] for r in res_gauss], [r['gesamtfehler'] for r in res_gauss], 'b-o', label='Gesamtfehler |Q(f_tilde) - W|')
plt.semilogy([r['N'] for r in res_gauss], [r['verfahrensfehler'] for r in res_gauss], 'g--s', label='Verfahrensfehler |Q(f) - W|')
plt.semilogy([r['N'] for r in res_gauss], [r['datenfehler'] for r in res_gauss], 'r--^', label='Datenfehler |Q(f_tilde) - Q(f)|')
plt.axhline(1.0, color='black', linestyle=':', label='Ziel-Toleranz (1 µJ)')
plt.title('Gauß-Legendre')
plt.xlabel('Anzahl Messpunkte N')
plt.ylabel('Fehler [µJ]')
plt.grid(True, which="both", ls="--")
plt.legend()

# Newton-Cotes Plot
plt.subplot(1, 2, 2)
plt.semilogy([r['N'] for r in res_nc], [r['gesamtfehler'] for r in res_nc], 'b-o', label='Gesamtfehler')
plt.semilogy([r['N'] for r in res_nc], [r['verfahrensfehler'] for r in res_nc], 'g--s', label='Verfahrensfehler')
plt.semilogy([r['N'] for r in res_nc], [r['datenfehler'] for r in res_nc], 'r--^', label='Datenfehler')
plt.axhline(1.0, color='black', linestyle=':', label='Ziel-Toleranz (1 µJ)')
plt.title('Newton-Cotes')
plt.xlabel('Anzahl Messpunkte N')
plt.ylabel('Fehler [µJ]')
plt.grid(True, which="both", ls="--")
plt.legend()

plt.tight_layout()

# Grafik mit pathlib im Ordner Task 2 speichern
plt.savefig(TASK_2 / 'aufgabe4_diagramm.png', dpi=300, bbox_inches='tight')
plt.show()

# Auswertung im Terminal
opt_gauss = next((r for r in res_gauss if r['gesamtfehler'] < 1.0), None)
opt_nc = next((r for r in res_nc if r['gesamtfehler'] < 1.0), None)

print("=== ERGEBNISSE AUFGABE 4 ===")
if opt_gauss:
    print(f"Gauß-Legendre braucht mindestens N = {opt_gauss['N']} Messpunkte für Gesamtfehler < 1 µJ (Fehler: {opt_gauss['gesamtfehler']:.4f} µJ).")
if opt_nc:
    print(f"Newton-Cotes braucht mindestens N = {opt_nc['N']} Messpunkte für Gesamtfehler < 1 µJ (Fehler: {opt_nc['gesamtfehler']:.4f} µJ).")