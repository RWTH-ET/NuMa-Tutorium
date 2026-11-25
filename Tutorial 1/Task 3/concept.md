# Aufgabe 3 – Kondition der Entzerrungsmatrix

## Vom Echo zum Frequenzfaktor

Beim Echo addieren wir verzögerte, gewichtete Kopien des Signals. Bei $h=(1,a)$ kommt zu jedem Signalwert ein Anteil des vorherigen Werts dazu. Nach dem Entfernen des zyklischen Präfixes gilt $y_j=x_j+a\,x_{(j-1)\bmod N}$. In Matrixform ist das $y=Cx$.

Die DFT macht aus dieser Summe mit Verschiebung eine Multiplikation pro Frequenz. Dazu ergänzen wir den Echo-Vektor $h$ auf Länge $N$ mit Nullen.

$$
H_k=\sum_{j=0}^{N-1} h_j e^{-2\pi i kj/N},\qquad
Y_k=H_kX_k,\qquad X_k=\frac{Y_k}{H_k}\ \text{für }H_k\ne0.
$$

Die Fourierwerte $H_k$ sind die Eigenwerte der Faltungsmatrix $C$. Eine inverse DFT würde die entzerrten Zeitwerte zurückgeben. In dieser Aufgabe berechnen wir anhand von $h$, wie empfindlich die Entzerrung wäre.

## Kondition und Rauschen

$$
\kappa_2(C)=\frac{\max_k |H_k|}{\min_k |H_k|}.
$$

Diese Kondition gilt bei invertierbarem $C$ auch für die Entzerrungsmatrix $C^{-1}$. Deren Frequenzfaktoren sind die Kehrwerte $1/H_k$. Die Konditionszahl gibt an, wie viel größer der **relative** Fehler des rekonstruierten Signals im ungünstigsten Fall gegenüber dem relativen Messfehler sein kann. Bei einem einzelnen Frequenzanteil wird Rauschen genauso wie das Signal durch $H_k$ geteilt. Ist etwa $H_k=0{,}01$, wächst ein absoluter Fehler von $0{,}001$ auf $0{,}1$.

Für $h=(1,a)$ gilt $H_k=1+a\,e^{-2\pi i k/N}$. Bei geradem $N$ und $-1<a\le 0$ gehören der kleinste und größte Betrag zu $k=0$ und $k=N/2$. So ergibt sich $\kappa_2(C)=(1-a)/(1+a)$.

- Für $a=-0{,}9$ ergibt sich $\kappa_2(C)=19$.
- Bei $a=-1$ verschwindet der konstante Anteil, weil $H_0=0$. Die Entzerrung ist dann mehrdeutig und $\kappa_2(C)=\infty$.
- Für den Echo-Vektor $(1,0{,}5,0{,}25,-0{,}25)$ aus Aufgabe 1b liefert das Programm bei $N=16$ ungefähr $3{,}34$.
