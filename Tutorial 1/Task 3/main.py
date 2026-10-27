import cmath
import math


# Werte aus Aufgabe 1b; hier bei Bedarf anpassen.
h = [1, 0.5, 0.25, -0.25]
N = 16

def kondition(h, N):
    # Wir untersuchen die Empfindlichkeit der Entzerrung anhand von h.
    # Ein empfangenes Signal y wird hier nicht rekonstruiert.
    betraege = []

    # 1. Den Echo-Vektor in den Frequenzbereich umrechnen: ein H_k pro k.
    for k in range(N):
        # DFT mit Signallaenge N: Fehlende Eintraege von h sind implizit null.
        # H_k ist der Faktor, mit dem der Kanal diesen Frequenzanteil gewichtet.
        h_dach = sum(
            wert * cmath.exp(-2j * math.pi * k * n / N) 
            for n, wert in enumerate(h)
        )

        # 2. Betrag von H_k sammeln: Er beschreibt die Staerke des Faktors,
        # unabhaengig von seinem Vorzeichen oder seiner komplexen Phase.
        betraege.append(abs(h_dach))

    # Ein Faktor null bedeutet: Dieser Frequenzanteil wurde ausgeloescht.
    # Die Division durch H_k ist dann unmoeglich; die Umkehrung nicht eindeutig.
    if min(betraege) == 0:
        return math.inf  # Singulaere Matrix: Die Kondition ist unendlich.

    # 3. Groessten mit kleinstem Betrag vergleichen: Kappa beschreibt die
    # groesstmoegliche Verstaerkung relativer Messfehler beim Entzerren.
    return max(betraege) / min(betraege)


def main():
    # Zuerst die Kondition fuer den gegebenen Echo-Vektor aus Aufgabe 1b.
    print(f"h = {h}, N = {N}")
    print(f"Kondition kappa: {kondition(h, N):.6g}")

    # Danach vergleichen: kein Echo (a = 0) bis zur Ausloeschung bei a = -1.
    print(f"\nFamilie h = (1, a), mit N = {N}:")
    print(f"{'a':>8} {'Kondition':>12}")
    for a in [0, -0.25, -0.5, -0.75, -0.9, -0.99, -0.999, -1]:
        print(f"{a:8.3f} {kondition([1, a], N):12.6g}")



if __name__ == "__main__":
    main()
