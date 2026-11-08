Wir haben mehrere Referenzpunkte gegeben, für die eine spannung einer Temperatur zugeordnet ist. Wir müssen jetzt eine Kurve bauen um jeder Spannung einer Temperatur zuzuordnen. Dazu verwenden wir Vandermonde + Horner und baryzentrische Interpolation

## Vandermonde + Horner:

**1. Vandermonde:**

Bei drei Referenzpunkten suchen wir ein Polynom vom Grad höchstens 2:

$$
p(E)=a_0+a_1E+a_2E^2.
$$

Die Koeffizienten $a_0,a_1,a_2$ sind noch unbekannt. Jeder bekannte Punkt
aus Spannung $E_i$ und Temperatur $T_i$ liefert eine Gleichung:

$$
\begin{aligned}
a_0+a_1E_0+a_2E_0^2 &= T_0,\\
a_0+a_1E_1+a_2E_1^2 &= T_1,\\
a_0+a_1E_2+a_2E_2^2 &= T_2.
\end{aligned}
$$

Die Faktoren $1,E_i,E_i^2$ bilden die Zeilen der Vandermonde-Matrix V. Es gilt also **V * a = T**
Die Lösung des lin. Gleichungssystem ist a, unsere Polynomkurve.


**2. Horner (Optimierung):**

Horner gibt umwandlung des Polynoms:

$$
p(E)=a_0+a_1E+a_2E^2=a_0+E(a_1+Ea_2).
$$
Python klammert die Formel nicht automatisch nach Horner um -> unnötige Rechenleistung, da weniger multiplikationen

## Baryzentrische Interpolation:

Zuerst Gewicht für jeden Referenzpunkt:

$$
w_j=\frac{1}{\prod_{k\ne j}(E_j-E_k)}.
$$

Bei drei Referenzpunkten ist zum Beispiel:

$$
w_0=\frac{1}{(E_0-E_1)(E_0-E_2)}.
$$

Sie hängen nur von den Referenzspannungen ab und werden einmal berechnet.
Für eine gesuchte Spannung $E$ berechnen wir dann für jeden Referenzpunkt einen Faktor:

$$
b_j=\frac{w_j}{E-E_j}.
$$

Mit diesen Faktoren kombinieren wir die bekannten Temperaturen:

$$
p(E)=\frac{b_0T_0+b_1T_1+b_2T_2}{b_0+b_1+b_2}.
$$

Also: Jeden Faktor mit seiner Referenztemperatur multiplizieren, alles addieren
und durch die Summe der Faktoren teilen. Das Ergebnis ist die gesuchte Temperatur.

- Trifft $E$ genau eine Referenzspannung, geben wir direkt deren Temperatur zurück, damit keine Division durch null entsteht.
- Es entsteht mathematisch dieselbe Polynomkurve wie bei Vandermonde.

## Kurzer Vergleich

| | Vandermonde + Horner | Baryzentrisch |
|---|---|---|
| Vorteil | Liefert die Polynomkoeffizienten; danach schnelle Auswertung mit Horner | Vermeidet das empfindliche Gleichungssystem für die Monomkoeffizienten |
| Nachteil | Bei hohen Graden können Rundungsfehler stark verstärkt werden | Braucht Stützpunkte und Gewichte; mehr Divisionen pro Auswertung |
| Sinnvoll, wenn ... | der Grad klein ist oder wir die Koeffizienten weiterverwenden wollen | wir das Interpolationspolynom möglichst zuverlässig auswerten wollen |

Für die reine Interpolation ist meist **baryzentrisch** vorzuziehen.
Die Verstärkung gerundeter Daten bei vielen äquidistanten Stützstellen bleibt
aber bei beiden Verfahren bestehen. Mehr Punkte bedeuten deshalb nicht automatisch mehr Genauigkeit.
