# Aufgaben 1–3 – Integrale aus wenigen Punkten berechnen

Eine Quadratur schätzt ein Integral aus Funktionswerten an Stützstellen $x_i$. Das Gewicht $w_i$ bestimmt, wie stark der jeweilige Wert eingeht.

$$
Q(f)=\sum_i w_i f(x_i),\qquad
x_i=\frac{a+b}{2}+\frac{b-a}{2}\hat x_i,\qquad
w_i=\frac{b-a}{2}\hat w_i.
$$

Die letzten beiden Formeln übertragen die Punkte und Gewichte vom Referenzintervall $[-1,1]$ auf $[a,b]$. Die Punkte werden verschoben und gestreckt, die Gewichte nur mit der halben Intervallbreite skaliert.

## 1 – Gauß-Legendre

Für $N$ Punkte bauen wir eine symmetrische Jacobi-Matrix mit Nullen auf der Hauptdiagonale und $\beta_k=(4-k^{-2})^{-1/2}$ auf den beiden Nebendiagonalen ($k=1,\ldots,N-1$). Ihre Eigenwerte sind die Stützstellen. Aus der ersten Komponente $u_{1i}$ des jeweiligen normierten Eigenvektors folgt $\hat w_i=2u_{1i}^2$.

Theoretisch ist die Formel bis Grad $m=2N-1$ exakt. Im Programm testen wir nacheinander $f(x)=x^k$. Das exakte Integral auf $[-1,1]$ ist $0$ bei ungeradem $k$ und $2/(k+1)$ bei geradem $k$. Beim ersten Unterschied über der Toleranz gilt $m_\mathrm{gemessen}=k-1$. Kleine Fehler unter der Toleranz können daher als Treffer erscheinen. [Diagramm](aufgabe1_gauss_legendre.png)

## 2 – Newton-Cotes

Wir verwenden für $N\ge2$ die gleichmäßig verteilten Punkte $\hat x_i=-1+2i/(N-1)$ für $i=0,\ldots,N-1$. Das Interpolationspolynom besteht aus den Funktionswerten mal den zugehörigen Lagrange-Basisfunktionen $\ell_i$. Beim Integrieren wird aus jeder Basisfunktion ein Gewicht.

$$
\hat w_i=\int_{-1}^{1}\ell_i(x)\,dx,\qquad
\kappa=\sum_i|\hat w_i|.
$$

Die Basisfunktion $\ell_i$ ist am eigenen Stützpunkt 1 und an den übrigen 0. Wir werten sie mit denselben baryzentrischen Faktoren wie bei der [Interpolation](../../Tutorial%202/Task%201/concept.md) aus und integrieren sie mit Gauß. Ihr Grad $N-1$ erklärt, warum $\lceil N/2\rceil$ Gauß-Punkte dafür theoretisch genügen. Newton-Cotes ist bis Grad $N-1$ bei geradem $N$ und bis Grad $N$ bei ungeradem $N$ exakt.

$\kappa$ begrenzt auf $[-1,1]$ die Fehlerverstärkung. Sind die Funktionswerte höchstens um $\varepsilon$ gestört, kann sich die Quadratur um höchstens $\kappa\varepsilon$ ändern. Bei Gauß sind alle Gewichte positiv und summieren sich zu 2. Bei Newton-Cotes kann die Summe der Beträge stark wachsen. Das erklärt die Empfindlichkeit gegenüber Rundung. [Vergleich](aufgabe2_vergleich_gauss_nc.png)

## 3 – Summierte Quadratur

Wir teilen $[0,1]$ in $M$ gleich große Stücke mit je $N$ Stützstellen. Gauß braucht $MN$ Auswertungen. Bei Newton-Cotes teilen sich Nachbarintervalle einen Randpunkt, also genügen $M(N-1)+1$ Auswertungen. Beim Budget 48 wählen wir daher $M=\lfloor48/N\rfloor$ für Gauß und $M=\lfloor47/(N-1)\rfloor$ für Newton-Cotes.

| Testfunktion | Exaktes Integral | Beste Aufteilung im Versuch $(N,M)$ |
| --- | --- | --- |
| $e^x$ | $e-1$ | Mehrere Aufteilungen nahe der Rechengenauigkeit |
| $1/(1+400(x-0{,}3)^2)$ | $(\arctan(14)+\arctan(6))/20$ | Gauß $(16,3)$, Newton-Cotes $(2,47)$ |
| $\lvert\sin(5x)\rvert$ | $(3+\cos(5))/5$ | Gauß $(23,2)$, Newton-Cotes $(2,47)$ |

Das Diagramm zeigt den absoluten Fehler $\lvert Q-I\rvert$ über $N$, wobei $I$ das exakte Integral aus der Tabelle ist. Die glatte Exponentialfunktion lässt sich gut annähern. Bei der schmalen Spitze kommt es darauf an, ob genügend Punkte nahe der Spitze liegen. Beim Knick bei $x=\pi/5$ spielt auch die Lage der Teilintervallgrenzen eine Rolle. [Fehlerdiagramm](aufgabe3_summierte_quadratur.png)
