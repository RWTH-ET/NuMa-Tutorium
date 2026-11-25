# Aufgabe 3b – Wie stark wirkt sich die Rundung aus?

Wir verwenden dieselben Stützspannungen für zwei baryzentrische Kurven. $p_\mathrm{exakt}$ verwendet die exakten Temperaturen an den Stützstellen, $p_\mathrm{kal}$ die gerundeten Kalibratorwerte. So lässt sich getrennt untersuchen, was vom Polynom und was von den gerundeten Eingaben kommt.

## Zwei Fehler auf demselben Prüfraster

$$
\begin{aligned}
e_\mathrm{interp}&=\max_E |p_\mathrm{exakt}(E)-T_\mathrm{exakt}(E)|,\\
e_\mathrm{daten}&=\max_E |p_\mathrm{kal}(E)-p_\mathrm{exakt}(E)|.
\end{aligned}
$$

Der **Interpolationsfehler** entsteht beim Annähern der echten Temperaturkurve durch ein Polynom. Der **Datenfehler** entsteht durch die geänderten Temperaturen an den Stützstellen. Das Maximum nehmen wir jeweils über dasselbe Prüfraster.

## Wieso kann die Rundung so stark wirken?

Der Kalibrator rundet auf $0{,}01\,^\circ\mathrm{C}$. Die Rundungsabweichung $\Delta T_j$ beträgt also höchstens $0{,}005\,^\circ\mathrm{C}$ pro Stützstelle. Im Polynom werden diese Abweichungen mit den Lagrange-Basisfunktionen $\ell_j$ gewichtet.

$$
\begin{aligned}
p_\mathrm{kal}(E)-p_\mathrm{exakt}(E)&=\sum_j \ell_j(E)\,\Delta T_j,\\
\Lambda_\mathrm{Raster}&=\max_{E\ \text{im Prüfraster}}\sum_j|\ell_j(E)|,\\
e_\mathrm{daten}&\le\Lambda_\mathrm{Raster}\cdot0{,}005\,^\circ\mathrm{C}.
\end{aligned}
$$

Die Lebesgue-Konstante $\Lambda_\mathrm{Raster}$ beschreibt hier die größtmögliche Verstärkung auf dem Prüfraster. Aus der baryzentrischen Formel in [Aufgabe 1](../Task%201/concept.md) folgt $\ell_j(E)=b_j(E)/\sum_i b_i(E)$ für $E$ abseits der Stützstellen. Genau an einer Stützstelle zählt nur der zugehörige Temperaturwert.

## Diagramme lesen

Die gestrichelte Linie zeigt die berechnete Obergrenze. Bei Grad 20 ist $\Lambda_\mathrm{Raster}$ ungefähr $11\,000$ für gleichmäßige Abstände und $2{,}66$ für Tschebyscheff. Der tatsächlich beobachtete Datenfehler liegt bei rund $11{,}1\,^\circ\mathrm{C}$ beziehungsweise $0{,}006\,^\circ\mathrm{C}$. Daran sieht man, wie stark die Lage der Stützstellen die Wirkung derselben Rundungsauflösung verändert. [Äquidistant](3b_aequidistant.png) · [Tschebyscheff](3b_tschebyscheff.png)
