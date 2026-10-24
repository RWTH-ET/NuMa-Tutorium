from cmath import exp, pi
from pathlib import Path


ORDNER = Path(__file__).resolve().parent
N = 16
PRAEFIX = 4
PILOTEN = [0, 4, 8, 12]


def lese_werte(dateiname):
    zeilen = (ORDNER / dateiname).read_text().splitlines()
    signal = []
    for zeile in zeilen[1:]:  # Kopfzeile "# re im" weglassen.
        re, im = zeile.split()
        signal.append(complex(float(re), float(im)))
    return signal


def dft(signal):
    # Fehlende Werte bis zur Laenge N sind implizit null.
    return [
        sum(wert * exp(-2j * pi * k * n / N)
            for n, wert in enumerate(signal))
        for k in range(N)
    ]


def b_schaetze_h(X, Y):
    # An den Piloten gilt Y[k] = H[k] * X[k].
    H_piloten = [Y[k] / X[k] for k in PILOTEN]
    h = []
    # Inverse 4-Punkt-DFT.
    for n in range(4):
        wert = sum(H_piloten[j] * exp(2j * pi * j * n / 4)
                   for j in range(4)) / 4
        h.append(wert)
    return h


def c_plotte_symbole(symbole):
    # Achsen von -4 bis 4, mit 60 Pixeln pro Einheit.
    svg = ['''<svg xmlns="http://www.w3.org/2000/svg" width="600" height="600">
        <rect width="600" height="600" fill="white"/>
        <text x="560" y="565">Re</text>
        <text x="20" y="45">Im</text>''']
    for wert in [-4, -2, 0, 2, 4]:
        x = 300 + 60 * wert
        y = 300 - 60 * wert
        svg.append(f'<path d="M {x} 60 V 540 M 60 {y} H 540" stroke="lightgray"/>')
        svg.append(f'<text x="{x}" y="565" text-anchor="middle">{wert}</text>')
        svg.append(f'<text x="40" y="{y + 5}" text-anchor="end">{wert}</text>')
    svg.append('<path d="M 60 300 H 540 M 300 60 V 540" stroke="black"/>')
    for z in symbole:
        x = 300 + 60 * z.real
        y = 300 - 60 * z.imag
        svg.append(f'<circle cx="{x}" cy="{y}" r="4"/>')
    svg.append('</svg>')
    datei = ORDNER / "symbole.svg"
    datei.write_text("\n".join(svg), encoding="utf-8")


def main():
    X_piloten = lese_werte("pilotsymbole.csv")

    # a) Praefix entfernen und Pilotwerte pruefen.
    X = dft(lese_werte("datensatz1.csv")[PRAEFIX:])
    # Toleranz fuer DFT-Rundungsfehler, auch auf unbesetzten Traegern.
    stimmen_ueberein = all(abs(X[k] - X_piloten[k]) < 1e-8 for k in range(N))
    print("Piloten korrekt:", stimmen_ueberein)

    # b) Echo aus den bekannten Pilotwerten bestimmen.
    Y = dft(lese_werte("datensatz2.csv")[PRAEFIX:])
    h = b_schaetze_h(X_piloten, Y)
    print("h:", ", ".join(f"{wert:.3f}" for wert in h))

    # c) Jeden Frequenzwert durch den Kanal-Frequenzwert teilen.
    Y = dft(lese_werte("datensatz3.csv")[PRAEFIX:])
    H = dft(h)
    symbole = [Y[k] / H[k] for k in range(N)]
    c_plotte_symbole(symbole)


if __name__ == "__main__":
    main()
