import numpy as np
import time

# Konfiguration für LTE und 5G NR
systems = {
    "LTE":   {"N": 2048, "ziel": 14000},
    "5G NR": {"N": 4096, "ziel": 28000}
}

print("=== OFDM-Symbole pro Sekunde (Optimierte NumPy FFT) ===\n")

runs = 1000  # Anzahl der Durchläufe für präzisere Zeitmessung

for name, config in systems.items():
    N = config["N"]
    ziel = config["ziel"]
    
    # Zufälliges Signal generieren
    x = np.random.rand(N) + 1j * np.random.rand(N)
    
    # Einmaliger Warm-up Run (initialisiert Caching im Hintergrund)
    np.fft.fft(x)
    
    # Exakte Zeitmessung
    start = time.perf_counter()
    for _ in range(runs):
        np.fft.fft(x)
    end = time.perf_counter()
    
    # Durchschnittliche Zeit und Symbole/s berechnen
    avg_time = (end - start) / runs
    symbole_pro_sekunde = int(1 / avg_time)
    
    print(f"{name} (N={N}, Ziel: {ziel} Symbole/s):")
    print(f"  Gemessene Zeit pro FFT: {avg_time * 1e6:.2f} µs")
    print(f"  Erreicht: {symbole_pro_sekunde:,} Symbole/s\n".replace(',', '.'))