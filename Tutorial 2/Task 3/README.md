"""PP2, Aufgabe 3a–3c: python main.py (benötigt NumPy)."""
from pathlib import Path
import csv
import numpy as np
from interpolation import baryzentrische_basis, vandermonde, plot, OUT

ROOT = Path(__file__).resolve().parent
SPALTEN = ['n', 'e_interp_C', 'e_daten_C', 'Lambda_pruefraster',
           'Schranke_daten_C', 'relative_Datenaenderung',
           'relative_Koeffizientenaenderung', 'e_koeff_C', 'kappa_inf',
           'e_gesamt_baryzentrisch_C']


def berechnen():
    """Alle Maxima werden über das vorgegebene Prüfraster gebildet."""
    raster, temperatur = np.loadtxt(ROOT / 'data/pruefraster.csv').T
    ergebnisse = {}
    for art in ['aequidistant', 'tschebyscheff']:
        daten = np.loadtxt(ROOT / f'data/stuetzwerte_{art}.csv')
        zeilen = []
        for n in range(2, 41):
            block = daten[daten[:, 0] == n]
            assert len(block) == n + 1, f'Falsche Stützstellenzahl bei n={n}'
            x, exakt, kalibrator = block[:, 2:].T
            assert np.max(np.abs(kalibrator - exakt)) <= .005 + 1e-10

            # 3a: exakte Stützwerte baryzentrisch interpolieren.
            basis = baryzentrische_basis(x, raster)
            p_exakt = basis @ exakt
            e_interp = np.max(np.abs(p_exakt - temperatur))

            # 3b: Auswirkung der Kalibratorrundung und Lebesgue-Schranke.
            p_kalibrator = basis @ kalibrator
            e_daten = np.max(np.abs(p_kalibrator - p_exakt))
            lam = np.max(np.sum(np.abs(basis), axis=1))
            schranke = lam * .005
            assert e_daten <= schranke * (1 + 1e-6) + 1e-8

            # 3c gilt ausschließlich für Tschebyscheff-Stützstellen.
            relativ_a = e_koeff = kappa = float('nan')
            if art == 'tschebyscheff':
                _, a, _ = vandermonde(x, exakt, raster)
                p_vander, a_tilde, kappa = vandermonde(x, kalibrator, raster)
                relativ_a = np.max(np.abs(a_tilde - a)) / np.max(np.abs(a))
                e_koeff = np.max(np.abs(p_vander - p_kalibrator))
                assert np.allclose(
                    p_vander - temperatur,
                    (p_exakt - temperatur) + (p_kalibrator - p_exakt)
                    + (p_vander - p_kalibrator), rtol=1e-10, atol=1e-8)
            relativ_y = np.max(np.abs(kalibrator - exakt)) / np.max(np.abs(exakt))
            zeilen.append([n, e_interp, e_daten, lam, schranke, relativ_y,
                           relativ_a, e_koeff, kappa,
                           np.max(np.abs(p_kalibrator - temperatur))])
        ergebnisse[art] = np.array(zeilen)
    return ergebnisse


def main():
    OUT.mkdir(exist_ok=True)
    resultate = berechnen()
    n = np.arange(2, 41)
    a, t = resultate['aequidistant'], resultate['tschebyscheff']
    plot('3a_interpolationsfehler.svg', '3a: Interpolationsfehler', n,
         [('Äquidistant', a[:, 1]), ('Tschebyscheff', t[:, 1])])
    for art, r in resultate.items():
        plot(f'3b_{art}.svg', f'3b: {art}', n,
             [('Interpolationsfehler', r[:, 1]), ('Datenfehler', r[:, 2]),
              ('Λ auf Prüfraster · 0,005 °C', r[:, 4])])
        with (OUT / f'{art}.csv').open('w', newline='', encoding='utf-8') as f:
            writer = csv.writer(f)
            writer.writerow(SPALTEN)
            writer.writerows(r)
    plot('3c_relative_aenderung.svg', '3c: Relative Änderungen', n,
         [('Monomkoeffizienten', t[:, 6]), ('Stützwerte', t[:, 5])],
         ylabel='Relative Änderung')
    plot('3c_koeffizientenfehler.svg', '3c: Koeffizientenfehler und Kondition', n,
         [('Koeffizientenfehler', t[:, 7])], right=('κ∞(V)', t[:, 8]))
    text = ['# Messauswertung – Aufgabe 3', '',
            'Maxima auf dem vorgegebenen Prüfraster; Temperaturen in °C.', '']
    for art, r in resultate.items():
        text += [f'## {art}', '', '| n | Interpolation | Daten | Schranke |',
                 '|---|---|---|---|']
        for grad in [2, 5, 10, 15, 20, 30, 40]:
            row = r[grad - 2]
            text.append(f'| {grad} | {row[1]:.6g} | {row[2]:.6g} | {row[4]:.6g} |')
        text.append('')
    text += ['## 3c – Tschebyscheff', '',
             '| n | Relative Koeffizientenänderung | e_koeff | κ∞ |',
             '|---|---|---|---|']
    for grad in [2, 5, 10, 15, 20, 30, 40]:
        row = t[grad - 2]
        text.append(f'| {grad} | {row[6]:.6g} | {row[7]:.6g} | {row[8]:.6g} |')
    good = t[t[:, 9] < .1, 0].astype(int)
    best = np.argmin(t[:, 9])
    text += ['', f'Baryzentrisch mit Tschebyscheff: Ziel erfüllt für n={good.tolist()}.',
             f'Kleinster Gesamtfehler: {t[best,9]:.8g} °C bei n={int(t[best,0])}.']
    (OUT / 'messauswertung.md').write_text('\n'.join(text) + '\n', encoding='utf-8')
    print('3a–3c abgeschlossen. Fünf Diagramme und Messwerte liegen in:', OUT)


if __name__ == '__main__':
    main()
