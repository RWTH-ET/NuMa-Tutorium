
Wir berechen das Polynom mit den gleichen Stützwerten mit der baryzentrischen Interpolation. Für die Temperatur bneutzen wir dabei aber einmal __den exakten Wert__ und einmal __den Kalibrator, der auf 0,01 °C gerundet ist__
Dabei vergleichen wir auch __Äquidistante-__ und __Tschebyscheff-Stützstellen__ also gleichmäßig verteilte und rand-dichte.

Wir vergleichen dann die Polynome und schauen uns nur die __maximale abweichung__ an. Das ist der Datenfehler (orange)
(blaue) Kurve ist der Fehler zwischen Polynom aus exakten Stüzwerten und der tatsächlich echten Temperatur.
(grün) Worst-Case Kurve. Lebesgue-Konstante aus Stützwerten * 0,005 °C

Formel für die Lebesgue-Konstante auf dem Prüfraster:

$$
\Lambda_n^{(\mathrm{Raster})} = \max_{E\in\mathrm{Prüfraster}}
\sum_{j=0}^{n}|\ell_j(E)|.
$$

$\ell_j(E)$ ist der Basiswert der $j$-ten Stützstelle an der Spannung $E$

## Wiederverwendung der Funktionen aus Task 1

`main.py` lädt `../Task 1/main.py` mit `importlib.util` als Modul `task1`.
Die Experimente aus Task 1 werden dabei nicht gestartet, da ihr Aufruf durch
`if __name__ == "__main__"` geschützt ist.

`task1.baryzentrische_gewichte(x)` berechnet die Gewichte einmal pro Stützstellensatz.
Mit denselben Gewichten wertet `task1.baryzentrisch(...)` das Polynom sowohl für
die exakten als auch für die gerundeten Temperaturen aus. Diese Funktionen
werden direkt aus Task 1 übernommen und hier nicht erneut implementiert.
Nur die Basiswerte zur Berechnung der Lebesgue-Konstante werden zusätzlich
in Task 3b bestimmt.

## => Woher kommt der Felher? Interpolation oder Rundung.

# Aufgabe 3: Woher kommt der Fehler?

Bei der Umrechnung von Thermospannung in Temperatur entstehen drei Fehleranteile:
Das Polynom nähert die echte Kennlinie nur an (**3a: Interpolationsfehler**),
der Kalibrator rundet die bekannten Temperaturen (**3b: Datenfehler**), und
der Rechner rundet bei der Berechnung (**3c: Koeffizientenfehler**).

## Datenfehler in 3b

Der Kalibrator zeigt Temperaturen auf $0{,}01\,°C$ gerundet an. Jeder Stützwert
weicht daher höchstens um $0{,}005\,°C$ vom exakten Wert ab. Wir bilden aus den
exakten und den gerundeten Werten jeweils ein Interpolationspolynom. Ihre größte
Abweichung auf dem Prüfraster ist der **Datenfehler**:

$$
e_{\mathrm{daten}} = \max_E |p_{\mathrm{kalibrator}}(E)-p_{\mathrm{exakt}}(E)|.
$$

Die **Lebesgue-Konstante** $\Lambda_n$ beschreibt, wie stark die Lage der
Stützstellen kleine Datenfehler höchstens verstärken kann. Anders als $\pi$
hat sie keinen festen Zahlenwert: Sie hängt von der Anzahl und der Lage der
Stützstellen ab. Wir bestimmen sie näherungsweise mit der Formel oben.
Zur Orientierung ergeben sich auf unserem Prüfraster ungefähr diese Werte:

| Polynomgrad | Äquidistant | Tschebyscheff |
|---|---:|---:|
| $n=10$ | $30$ | $2{,}44$ |
| $n=20$ | $11\,000$ | $2{,}66$ |

Je größer $\Lambda_n$ ist, desto stärker *können* Rundungsfehler verstärkt
werden. Für den Datenfehler gilt auf dem Prüfraster die Schranke

$$
e_{\mathrm{daten}} \leq \Lambda_n^{(\mathrm{Raster})}\cdot 0{,}005\,°C.
$$

Wir vergleichen äquidistante und Tschebyscheff-Stützstellen. Bei vielen
äquidistanten Punkten kann der Datenfehler stark wachsen. Mit
Tschebyscheff-Stützstellen bleibt er deutlich kleiner, verschwindet durch
mehr Punkte aber nicht automatisch. Die berechnete Schranke bezieht sich
auf das Prüfraster, nicht auf jeden Punkt des gesamten Intervalls.
