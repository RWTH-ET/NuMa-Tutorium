# Aufgabe 4 – Kondensatorenergie aus möglichst wenigen Messungen

## Von der Kapazität zur Energie

Die Kapazität $C(v)$ hängt von der Spannung $v$ ab. Zur Energie trägt jede kleine Spannungsänderung mit $vC(v)\,dv$ bei. Für $C$ in Mikrofarad und $v$ in Volt ergibt sich.

$$
C(v)=\frac{10}{1+(v/8)^2},\qquad
W=\int_0^{12}vC(v)\,dv=320\ln(3{,}25)\,\mu\mathrm{J}
\approx377{,}17\,\mu\mathrm{J}.
$$

Für den exakten Wert setzen wir $u=1+(v/8)^2$. Damit ist $du=v\,dv/32$, und die Grenzen werden zu $1$ und $3{,}25$.

Eine Messung kostet Zeit. Daher vergleichen wir die Verfahren aus [Aufgaben 1–3](../Task%201/concept.md) bei verschiedenen Anzahlen $N$ an Messpunkten. Das Messgerät liefert $C_\mathrm{mess}(v)$, gerundet auf $0{,}01\,\mu\mathrm{F}$. Aus denselben Stützstellen und Gewichten berechnen wir einmal $Q_\mathrm{exakt}=\sum_i w_i v_i C(v_i)$ und einmal $Q_\mathrm{mess}=\sum_i w_i v_i C_\mathrm{mess}(v_i)$.

## Fehler auseinanderhalten

| Kurve | Betrag | Bedeutung |
| --- | --- | --- |
| Grün | $\lvert Q_\mathrm{exakt}-W\rvert$ | Verfahrensfehler durch die Näherung |
| Rot | $\lvert Q_\mathrm{mess}-Q_\mathrm{exakt}\rvert$ | Datenfehler durch die gerundeten Werte |
| Blau | $\lvert Q_\mathrm{mess}-W\rvert$ | Gesamtfehler gegenüber der exakten Energie |

Vor der Betragsbildung addieren sich Verfahrens- und Datenfehler. Haben sie verschiedene Vorzeichen, können sie sich teilweise aufheben.

## Ergebnis und Garantie

Die blaue Kurve unterschreitet $1\,\mu\mathrm{J}$ zuerst bei **3 Gauß-Punkten** ($0{,}425\,\mu\mathrm{J}$) und **5 Newton-Cotes-Punkten** ($0{,}366\,\mu\mathrm{J}$). [Diagramm](aufgabe4_diagramm.png)

Jeder gerundete Kapazitätswert kann um höchstens $0{,}005\,\mu\mathrm{F}$ abweichen. Deshalb ist der Datenfehler höchstens $0{,}005\sum_i |w_i v_i|$. Bei diesen beiden Punktzahlen sind die Gewichte positiv und beide Formeln integrieren $v$ exakt, also ist die Summe $72$. Die Rundung trägt somit höchstens $0{,}36\,\mu\mathrm{J}$ bei. Zusammen mit dem Verfahrensfehler liegen die Obergrenzen bei etwa $0{,}838\,\mu\mathrm{J}$ für Gauß und $0{,}801\,\mu\mathrm{J}$ für Newton-Cotes. Unter diesem Messmodell reicht das für das Ziel.
