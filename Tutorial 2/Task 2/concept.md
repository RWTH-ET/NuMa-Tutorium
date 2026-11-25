# Aufgabe 2 – Tschebyscheff-Stützstellen

Wir wiederholen Aufgabe 1 mit denselben zwei Rechenwegen. Der Unterschied im Versuchsaufbau liegt in der Verteilung der $n+1$ Stützstellen. Bei Tschebyscheff liegen die Punkte näher am Rand dichter. Das dämpft starke Schwingungen der Temperaturkurve dort, wo der Fehler bei gleichmäßigen Abständen groß werden kann.

## Wo liegen die Punkte?

Die vorgegebenen Spannungen entstehen aus den Tschebyscheff-Punkten $t_j$ auf $[-1,1]$, übertragen auf den Spannungsbereich $[a,b]$.

$$
t_j=\cos\!\left(\frac{(2j+1)\pi}{2(n+1)}\right),\qquad
E_j=\frac{a+b}{2}+\frac{b-a}{2}t_j,\qquad j=0,\ldots,n.
$$

Die Winkel sind gleichmäßig verteilt. Weil sich der Kosinus nahe $0$ und $\pi$ nur langsam ändert, rücken die Punkte an den beiden Enden zusammen. Das Programm liest die Punkte und ihre gerundeten Temperaturen aus der CSV. Für Grad 2 liegen sie bei ungefähr $1{,}38$, $10{,}32$ und $19{,}26\,\mathrm{mV}$. Anschließend berechnen wir die Kurve wie in [Aufgabe 1](../Task%201/concept.md) und bestimmen den größten Betrag der Abweichung.

$$
e_{\max}(n)=\max_{E\ \text{im Prüfraster}} |p_n(E)-T_\mathrm{exakt}(E)|.
$$

## Diagramm lesen

Unter der gestrichelten Linie liegt der maximale Fehler unter $0{,}1\,^\circ\mathrm{C}$. Im Versuch erreichen beide Verfahren das Ziel erstmals bei Grad 9. Die baryzentrische Rechnung bleibt bis Grad 40 darunter. Bei Grad 40 liegt ihr Fehler bei rund $0{,}007\,^\circ\mathrm{C}$, bei Vandermonde mit Horner bei rund $0{,}31\,^\circ\mathrm{C}$. Die günstige Punktverteilung hilft beiden Verfahren. Bei hohem Grad ist die Berechnung der Vandermonde-Koeffizienten weiterhin empfindlich. [Diagramm](aufgabe2_tschebyscheff.png)
