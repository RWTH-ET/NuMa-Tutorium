# Aufgabe 3 – Kondition der Entzerrungsmatrix

**1 · Was macht das Echo?**

$$
\underbrace{x}_{\text{Eingang}}
\xrightarrow{\;\text{Echo }h\;}
\underbrace{y=h*x=Cx}_{\text{Ausgang}}
$$

Faltung = verzögerte, gewichtete Kopien addieren. Bei $h=(1)$ gibt es kein Echo.

$$
h=(1,a)\quad\Rightarrow\quad y_j=x_j+a\,x_{j-1}
$$

Die Position in $h$ bestimmt die Verzögerung, der Wert die Gewichtung.
Bei $a<0$ ist das Echo zusätzlich im Vorzeichen umgekehrt.

**2 · Warum eine zirkuläre Matrix?**

Nach Entfernen eines ausreichend langen zyklischen Präfixes laufen die Indizes im Kreis: $x_{-1}=x_{N-1}$.
Für $N=4$ und $h=(1,a)$:

$$
C=\begin{pmatrix}
1&0&0&a\\
a&1&0&0\\
0&a&1&0\\
0&0&a&1
\end{pmatrix}
\qquad
\begin{array}{l}
C:\ \text{Echo hinzufügen}\\
C^{-1}:\ \text{Echo entfernen}
\end{array}
$$

Die erste Zeile bedeutet $y_0=x_0+a x_3$.
Spalten werden mit h gefüllt und durchgeschoben stück für stück.

**3 · Fourier macht die Faltung einfach**

$$
\boxed{y=h*x}
\quad\xrightarrow{\text{DFT}}\quad
\boxed{\widehat y_k=\widehat h_k\widehat x_k}
\quad\xrightarrow{\text{Entzerren}}\quad
\boxed{\widehat x_k=\frac{\widehat y_k}{\widehat h_k}}
$$

Für das Zeitsignal anschließend invers transformieren.

$$
F_{k,n}=e^{-2\pi i kn/N},\qquad
C=F^{-1}\operatorname{diag}(\widehat h)F
$$

Von rechts lesen: **DFT → jede Frequenz gewichten → inverse DFT**.
$h$ wird auf die **Signallänge $N$** mit Nullen ergänzt, nicht nur auf seine eigene Länge transformiert.

**4 · Was misst die Kondition?** -> maximal schlimmst mögliche Verzerrung

$\|C\|_2$ ist der größte Streckungsfaktor von $C$.
Weil $F/\sqrt N$ Längen erhält, sind die Singularwerte von $C$ genau $|\widehat h_k|$.

$$
\boxed{\kappa_2(C)=\|C\|_2\|C^{-1}\|_2
=\frac{\max_k|\widehat h_k|}{\min_k|\widehat h_k|}}
$$

Die Faltungsmatrix $C$ und die Entzerrungsmatrix $C^{-1}$ haben dieselbe Kondition.

| Kondition | Bedeutung |
|---|---|
| $\kappa=1$ | Bestmögliche Kondition |
| $\kappa\gg1$ | Relative Fehler können stark verstärkt werden |
| $\kappa=\infty$ | Ein Fourierwert ist null: Information fehlt, keine eindeutige Umkehrung |

**5 · Beispiel: $h=(1,a)$**

$$
\widehat h_k=1+a\,e^{-2\pi i k/N},\qquad \widehat h_0=1+a
$$

Für gerades $N$ und $-1<a\leq0$ gilt $\kappa=(1-a)/(1+a)$:

| $a$ | $0$ | $-0{,}5$ | $-0{,}9$ | $-0{,}99$ | $-1$ |
|---|---|---|---|---|---|
| $\kappa$ | $1$ | $3$ | $19$ | $199$ | $\infty$ |

**$a\to-1$ → konstanter Signalanteil wird ausgelöscht → Entzerrung wird empfindlich.**

**6 · Warum stört Rauschen beim Entzerren?**

$$
\frac{\widehat y_k+\widehat\eta_k}{\widehat h_k}
=\widehat x_k+\underbrace{\frac{\widehat\eta_k}{\widehat h_k}}_{\text{mitverstärktes Rauschen}}
$$

Beispiel: $|\widehat h_k|=0{,}01$ → Rauschen dieser Frequenz wird **100-fach** verstärkt.

$$
\underbrace{\frac{\|\delta x\|_2}{\|x\|_2}}_{\text{relativer Rekonstruktionsfehler}}
\leq\kappa_2(C)\,
\underbrace{\frac{\|\eta\|_2}{\|y\|_2}}_{\text{relativer Messfehler}}
$$

Das ist eine **obere Schranke**, kein garantierter Fehler.
Für Aufgabe 1b: $h=(1,\,0{,}5,\,0{,}25,\,-0{,}25)$, $N=16$ → $\kappa\approx3{,}34227$.
