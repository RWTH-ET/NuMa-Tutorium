# Aufgabe 1 – Temperatur aus Spannung bestimmen

Aus $n+1$ Referenzpaaren $(E_j,T_j)$ bestimmen wir ein Polynom $p$ vom Grad höchstens $n$. Es soll durch alle Stützpunkte gehen, also $p(E_j)=T_j$ erfüllen. Damit schätzen wir die Temperatur zwischen den bekannten Spannungen. Hier sind die Stützstellen gleichmäßig verteilt.

## 1a – Vandermonde und Horner

Wir schreiben $p(E)=a_0+a_1 E+\dots+a_n E^n$. Jeder bekannte Punkt liefert eine Gleichung für die unbekannten Koeffizienten. Die Zeile $j$ der Vandermonde-Matrix enthält deshalb die Potenzen von $E_j$.

$$
V_{j,k}=E_j^k,\qquad V a=T,\qquad
\kappa_\infty(V)=\|V\|_\infty\|V^{-1}\|_\infty.
$$

Nach dem Lösen des Systems setzt Horner die gesuchte Spannung durch wiederholtes Multiplizieren und Addieren ein. Bei Grad 2 wird aus $a_0+a_1 E+a_2 E^2$ zum Beispiel $a_0+E(a_1+Ea_2)$.

Die Norm $\|\cdot\|_\infty$ ist die größte Summe der Beträge in einer Matrixzeile. Eine große Konditionszahl bedeutet, dass sich kleine Änderungen der Eingangsdaten stark auf die berechneten Koeffizienten auswirken können.

## 1b – Baryzentrische Interpolation

Wir berechnen Gewichte aus den Referenzspannungen. An einer gesuchten Spannung $E$ kombinieren wir damit die bekannten Temperaturen.

$$
w_j=\frac{1}{\prod_{k\ne j}(E_j-E_k)},\qquad
b_j(E)=\frac{w_j}{E-E_j},\qquad
p(E)=\frac{\sum_j b_j(E)T_j}{\sum_j b_j(E)}.
$$

Bei $E=E_j$ geben wir direkt $T_j$ zurück. Beide Rechenwege liefern in exakter Rechnung dasselbe Polynom. Die baryzentrische Formel kommt ohne das Gleichungssystem für die Koeffizienten aus.

## Diagramme lesen

Für jeden Grad vergleichen wir $p(E)$ auf dem Prüfraster mit der Referenztemperatur $T_\mathrm{exakt}(E)$. Aufgetragen ist $e_{\max}=\max_E|p(E)-T_\mathrm{exakt}(E)|$. Beide Verfahren erreichen im Versuch bei Grad 14 rund $0{,}123\,^\circ\mathrm{C}$ als kleinsten Fehler und liegen damit über dem Ziel von $0{,}1\,^\circ\mathrm{C}$. Bei Grad 20 zeigt die Abweichungskurve große Fehler in Randnähe. [Fehlervergleich](1b_fehlervergleich.png) · [Abweichungen bei Grad 20](1b_abweichung_n20.png) · [Konditionszahl](1a_kondition.png)
