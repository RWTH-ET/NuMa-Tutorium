import cmath
import math


# Werte aus Aufgabe 1b; hier bei Bedarf anpassen.
h = [1, 0.5, 0.25, -0.25]
N = 16

def kondition(h, N):
    # DFT mit Laenge N: Fehlende Werte von h sind implizit null.
    betraege = []
    for k in range(N):
        h_dach = sum(
            wert * cmath.exp(-2j * math.pi * k * n / N) 
            for n, wert in enumerate(h)
        )
        betraege.append(abs(h_dach))

    if min(betraege) == 0:
        return math.inf  # Singulaere Matrix: Die Kondition ist unendlich.
    return max(betraege) / min(betraege)


def main():
    print(f"h = {h}, N = {N}")
    print(f"Kondition kappa: {kondition(h, N):.6g}")

    print(f"\nFamilie h = (1, a), mit N = {N}:")
    print(f"{'a':>8} {'Kondition':>12}")
    for a in [0, -0.25, -0.5, -0.75, -0.9, -0.99, -0.999, -1]:
        print(f"{a:8.3f} {kondition([1, a], N):12.6g}")



if __name__ == "__main__":
    main()
