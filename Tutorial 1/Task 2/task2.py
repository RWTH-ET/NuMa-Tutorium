import numpy as np                          # np ist für numerische Berechnungen und Arrays zuständig
import time
import matplotlib.pyplot as plt             # plt ist für die Visualisierung der Ergebnisse zuständig

# 1. Algorithmen definieren
def dft(x):
    """Direkte DFT in O(N^2) mittels Matrixmultiplikation."""
    N = len(x)                             
    n = np.arange(N)
    k = n.reshape((N, 1))
    e = np.exp(-2j * np.pi * k * n / N)
    return np.dot(e, x)

def fft(x):
    """Rekursive Cooley-Tukey FFT in O(N log N)."""
    N = len(x)
    if N <= 1:
        return x
    even = fft(x[0::2])
    odd = fft(x[1::2])
    T = np.exp(-2j * np.pi * np.arange(N // 2) / N) * odd
    return np.concatenate([even + T, even - T])

# 2. Experiment Setup
N_values = [2**p for p in range(4, 13)] # N = 16 bis 4096
times_dft = []
times_fft = []

time_lte_dft, time_lte_fft = 0, 0
time_5g_dft, time_5g_fft = 0, 0

print(f"{'N':<6} | {'DFT Zeit (s)':<15} | {'FFT Zeit (s)':<15}")
print("-" * 42)

for N in N_values:
    # Zufälliges komplexes Signal generieren
    x = np.random.rand(N) + 1j * np.random.rand(N)
    
    # DFT Messung
    start = time.perf_counter()
    dft(x)
    t_dft = time.perf_counter() - start
    times_dft.append(t_dft)
    
    # FFT Messung
    start = time.perf_counter()
    fft(x)
    t_fft = time.perf_counter() - start
    times_fft.append(t_fft)
    
    print(f"{N:<6} | {t_dft:<15.6f} | {t_fft:<15.6f}")
    
    # Werte für OFDM-Abschätzung speichern
    if N == 2048:
        time_lte_dft, time_lte_fft = t_dft, t_fft
    elif N == 4096:
        time_5g_dft, time_5g_fft = t_dft, t_fft

# 3. Empirische Steigungen berechnen (m = delta_y / delta_x im Log-Raum)
log_N = np.log10(N_values)
slope_dft, _ = np.polyfit(log_N, np.log10(times_dft), 1)    #Legt die Steigung der DFT-Zeit im Log-Log-Raum fest
slope_fft, _ = np.polyfit(log_N, np.log10(times_fft), 1)    #_ ist eine Wegwerfvariable, die den Achsenabschnitt enthält, der hier nicht benötigt wird

# 4. OFDM-Abschätzung
print("\n=== OFDM-Symbole pro Sekunde ===")
print(f"LTE (N=2048, Ziel: 14.000 Symbole/s):")
print(f"  DFT: {int(1/time_lte_dft)} Symbole/s")            #Kehrt die Zeit in Symbole pro Sekunde um
print(f"  FFT: {int(1/time_lte_fft)} Symbole/s")

print(f"\n5G NR (N=4096, Ziel: 28.000 Symbole/s):")
print(f"  DFT: {int(1/time_5g_dft)} Symbole/s")
print(f"  FFT: {int(1/time_5g_fft)} Symbole/s")

# 5. Visualisierung (log-log Plot)
plt.figure(figsize=(9, 6))
plt.loglog(N_values, times_dft, marker='o', label=f'Direkte DFT (Steigung $\\approx$ {slope_dft:.2f})')
plt.loglog(N_values, times_fft, marker='s', label=f'FFT (Steigung $\\approx$ {slope_fft:.2f})')

plt.xlabel('Vektorlänge N')
plt.ylabel('Laufzeit in Sekunden')
plt.title('Laufzeitexperiment: DFT vs. FFT')
plt.grid(True, which="both", ls="--", alpha=0.6)
plt.legend()
plt.tight_layout()
plt.show()