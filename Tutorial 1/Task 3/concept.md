# Aufgabe 3 – Kondition der Entzerrungsmatrix

**Kurzdefinition: Diskrete Fourier-Transformation (DFT)**

Die DFT wandelt einen Vektor $z=(z_0,\ldots,z_{N-1})$ mit $N$ Zeitwerten in einen Vektor $Z=(Z_0,\ldots,Z_{N-1})$ mit $N$ Frequenzwerten um:

$$
\boxed{Z_k=\sum_{n=0}^{N-1}z_n\,e^{-2\pi i kn/N}},
\qquad k=0,\ldots,N-1.
$$

- $n$: Index des Zeitwerts; $k$: Index des Frequenzanteils.
- $i$: imaginäre Einheit mit $i^2=-1$. Der Exponentialterm ist ein komplexer Gewichtungsfaktor für die jeweilige Frequenz.
- Für jedes $k$ multiplizieren wir **alle Zeitwerte** mit ihren Gewichtungsfaktoren und addieren die Produkte. So entsteht jeweils **eine Zahl $Z_k$**; alle zusammen bilden den Frequenzvektor $Z$.

Für unseren Echo-Vektor heißt das $H=\operatorname{DFT}(h)$. Ist $h$ kürzer als die Signallänge $N$, ergänzen wir ihn vorher mit Nullen. Bei $h=(1,a,0,\ldots,0)$ bleiben in der Summe nur zwei Beiträge:

$$
H_k=1+a\,e^{-2\pi i k/N}.
\qquad\text{Insbesondere: }H_0=1+a.
$$

Bei $k=0$ sind alle Gewichtungsfaktoren gleich $1$: Der nullte Fourierwert ist einfach die Summe der Vektoreinträge. Deshalb ist für $a=-1$ der Wert $H_0=1-1=0$.

Kappa = Worst Case der bei Rekonstruierung auftreten kann. Dabei gilt:  
Kappa = 1 -> fehler werden nicht verstärkt
Kapp a= 100 -> Fehler können bis faktor 100 verstärkt werden  
Kappa = inf -> Signal ist verschwunden, nicht rekonstruirbar

 Rauschen ist additiv und unabhängig vom signal/echo.
Rauschen wird verstärk, abhängig vom kleinsten Betrag der Fourierwerte des Echo-Vektors 

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

Wir kennen das empfangene Signal $y$ und den Echo-Vektor $h$. Gesucht ist das ursprüngliche Signal $x$, also das Signal **vor dem Echo**. Bisher steht dafür $y=Cx$. Jetzt wollen wir verstehen, wie man diese Gleichung einfach umkehrt und wie empfindlich das Ergebnis auf Fehler in $y$ reagiert.

### 3.1 Rechnen wir beim Echo mit Multiplikation oder Addition?

**Mit beidem: Zuerst wird jede verzögerte Kopie mit ihrem Echo-Faktor multipliziert. Anschließend werden die Beiträge an jeder Zeitstelle addiert.**

Für $h=(1,0{,}5)$ heißt das:

$$
y_j=\underbrace{1\cdot x_j}_{\text{Original an Stelle }j}
    +\underbrace{0{,}5\cdot x_{j-1}}_{\text{verzögerte Kopie an derselben Stelle }j}.
$$

Beispiel mit $N=4$ und $x=(2,4,6,8)$: In unserer **zyklischen** Rechnung ist $x_{-1}=x_3=8$.

| Stelle | Rechnung | Ergebnis |
|---|---|---|
| $j=0$ | $1\cdot2+0{,}5\cdot8$ | $y_0=6$ |
| $j=1$ | $1\cdot4+0{,}5\cdot2$ | $y_1=5$ |
| $j=2$ | $1\cdot6+0{,}5\cdot4$ | $y_2=8$ |
| $j=3$ | $1\cdot8+0{,}5\cdot6$ | $y_3=11$ |

Also ist $y=(6,5,8,11)$. Genau diese Rechnungen führt auch die Matrixmultiplikation $Cx$ aus: **Zeilenweise multiplizieren und die Produkte aufsummieren.**

Allgemein, wenn $h$ insgesamt $L$ Einträge hat:

$$
\boxed{y_j=\sum_{\ell=0}^{L-1}h_\ell\,x_{(j-\ell)\bmod N}}.
$$

Diese Operation heißt **zyklische Faltung** und wird kurz $h*x$ geschrieben. Das Sternchen bezeichnet hier die ganze Summe aus gewichteten, verschobenen Werten. Es bedeutet also **keine eintragsweise Multiplikation** von $h$ und $x$. Ein Verzögerungsschritt entspricht dabei einem Abtastschritt.

### 3.2 Was hat das mit den Eigenvektoren zu tun?

Du kennst bereits die Eigenvektorgleichung $Cv_k=\lambda_k v_k$: Auf einem Eigenvektor wirkt die Matrix besonders einfach. Sie multipliziert ihn nur mit einer Zahl, dem Eigenwert $\lambda_k$.

Die Fourier-Vektoren sind Eigenvektoren der zirkulären Matrix $C$. Man kann jedes Signal als Summe dieser Vektoren schreiben:

$$
x=\sum_{k=0}^{N-1}c_kv_k
\quad\Rightarrow\quad
Cx=\sum_{k=0}^{N-1}c_k\lambda_kv_k.
$$

**Wir stellen das Signal also in einer Basis dar, in der das Echo jeden Anteil einzeln mit einer Zahl multipliziert.** Die Anteile werden in dieser Basis nicht miteinander vermischt.

Genau diesen Basiswechsel erledigt die **diskrete Fourier-Transformation (DFT)**. Ein $\widehat x_k$ beschreibt den Anteil des Signals am Fourier-Muster mit Index $k$. Das Dach bedeutet „Fourier-transformiert“. $k$ ist jetzt ein **Frequenzindex**, keine Zeitstelle. Bei der hier verwendeten Normierung gilt $c_k=\widehat x_k/N$, wenn $v_k$ die unnormierten Fourier-Vektoren sind.

Die zugehörigen Eigenwerte sind die Fourierwerte des Echo-Vektors:

$$
\boxed{\lambda_k=\widehat h_k}.
$$

Darum lautet dieselbe Echo-Rechnung nach dem Basiswechsel:

$$
\boxed{\widehat y_k=\widehat h_k\,\widehat x_k}.
$$

Hier ist es tatsächlich **eine einzelne Multiplikation pro Frequenzindex**. Die Summen und Verschiebungen aus der Zeitdarstellung stecken bereits im Basiswechsel. Beide Darstellungen beschreiben dasselbe Echo.

### 3.3 Wie berechnen wir die Fourierwerte?

Für einen Vektor $z$ der Länge $N$ lautet die DFT:

$$
\widehat z_k=\sum_{n=0}^{N-1}z_n e^{-2\pi i kn/N},
\qquad k=0,\ldots,N-1.
$$

Das wenden wir auf $x$, $y$ und $h$ an. **Auch $h$ wird auf die Signallänge $N$ mit Nullen ergänzt.** Bei $h=(1,0{,}5)$ und $N=4$ transformieren wir also $(1,0{,}5,0,0)$ mit einer 4-Punkt-DFT.

Weil die letzten beiden Einträge null sind, bleibt für dieses $h$ nur:

$$
\widehat h_k=1+0{,}5e^{-2\pi i k/4}.
$$

Für unser Beispiel ergeben sich folgende Werte. $i$ ist die imaginäre Einheit mit $i^2=-1$; die komplexen Zahlen entstehen durch die DFT.

| $k$ | $\widehat h_k$ | $\widehat x_k$ | $\widehat y_k=\widehat h_k\widehat x_k$ |
|---|---|---|---|
| $0$ | $1{,}5$ | $20$ | $30$ |
| $1$ | $1-0{,}5i$ | $-4+4i$ | $-2+6i$ |
| $2$ | $0{,}5$ | $-4$ | $-2$ |
| $3$ | $1+0{,}5i$ | $-4-4i$ | $-2-6i$ |

Zum Beispiel ist bei $k=0$ die Rechnung einfach $30=1{,}5\cdot20$. Bei $k=1$ gilt dieselbe Regel mit komplexen Zahlen: $(1-0{,}5i)(-4+4i)=-2+6i$.

### 3.4 Wie entfernen wir das Echo wieder?

Wir kennen $y$ und $h$ und berechnen daraus $\widehat y$ und $\widehat h$. Dann lösen wir jede einzelne Gleichung durch **Division** nach $\widehat x_k$ auf:

$$
\boxed{\widehat x_k=\frac{\widehat y_k}{\widehat h_k}}
\qquad\text{sofern }\widehat h_k\neq0.
$$

Im Beispiel: $\widehat x_0=30/1{,}5=20$ und $\widehat x_2=(-2)/0{,}5=-4$. Das machen wir für alle vier Frequenzindizes.

Die so erhaltenen Werte sind noch die **Fourierwerte** von $x$. Für die ursprünglichen Zeitwerte wenden wir anschließend die inverse DFT an:

$$
x_j=\frac1N\sum_{k=0}^{N-1}\widehat x_k e^{2\pi i kj/N}.
$$

Damit erhalten wir wieder $x=(2,4,6,8)$. Beispielsweise ist $x_0=(20+(-4+4i)-4+(-4-4i))/4=2$.

**Der Rechenweg zum Entzerren ist also: $y$ und $h$ transformieren → Fourierwerte paarweise dividieren → das Ergebnis invers transformieren.**

### 3.5 Was bedeutet die Matrixformel dazu?

Mit der DFT-Matrix $F$ gilt:

$$
F_{k,n}=e^{-2\pi i kn/N},\qquad
D=\operatorname{diag}(\widehat h_0,\ldots,\widehat h_{N-1}),\qquad
\boxed{C=F^{-1}DF}.
$$

$D$ ist eine Diagonalmatrix: Nur auf der Diagonale stehen die Werte $\widehat h_k$, alle anderen Einträge sind null. Deshalb multipliziert $D$ jeden Fourierwert einzeln mit seinem Echo-Faktor.

Bei $Cx=F^{-1}DFx$ lesen wir **von rechts nach links**:

1. $Fx$: das Signal in Fourier-Anteile zerlegen.
2. $D(Fx)$: jeden Anteil mit $\widehat h_k$ multiplizieren.
3. $F^{-1}(DFx)$: zurück ins Zeitsignal umrechnen.

Zum Entzerren ersetzen wir die Faktoren durch ihre Kehrwerte:

$$
C^{-1}=F^{-1}\operatorname{diag}\left(\frac1{\widehat h_0},\ldots,\frac1{\widehat h_{N-1}}\right)F.
$$

Das funktioniert nur, wenn keiner der Werte $\widehat h_k$ null ist.

**4 · Was misst die Kondition?**

Wir können das Echo jetzt theoretisch entfernen. Die nächste Frage ist: **Wie stark können kleine Messfehler beim Entzerren den berechneten Vektor $x$ verfälschen?** Das beschreibt die Kondition. Sie ist keine Angabe darüber, wie laut das Echo ist.

Jeder Fourier-Anteil wird beim Hinzufügen des Echos mit $\widehat h_k$ multipliziert. Sein Betrag wird deshalb mit $|\widehat h_k|$ skaliert. Die Betragsstriche sind nötig, weil die Fourierwerte auch negativ oder komplex sein können.

- Der größte Betrag $\max_k|\widehat h_k|$ beschreibt die stärkste Streckung durch $C$.
- Der kleinste Betrag $\min_k|\widehat h_k|$ beschreibt den Anteil, der am stärksten abgeschwächt wird.
- Beim Entzerren wird durch diese Werte dividiert. Ein sehr kleiner Betrag führt dabei zu einer sehr großen Verstärkung, auch für Messfehler.

Die Kondition vergleicht den größten mit dem kleinsten Streckungsfaktor:

$$
\boxed{\kappa_2(C)=\|C\|_2\|C^{-1}\|_2
=\frac{\max_k|\widehat h_k|}{\min_k|\widehat h_k|}}
$$

**So rechnen wir: Alle $\widehat h_k$ bestimmen → ihre Beträge bilden → größten Betrag durch kleinsten Betrag teilen.** Die volle Matrix $C$ müssen wir dafür nicht aufbauen oder invertieren.

Im Beispiel $h=(1,0{,}5)$, $N=4$ sind die Beträge $1{,}5$, $\sqrt{1{,}25}$, $0{,}5$, $\sqrt{1{,}25}$. Also ist $\kappa_2(C)=1{,}5/0{,}5=3$.

Zur Verbindung mit der Normschreibweise: $\|C\|_2$ ist der größte Streckungsfaktor, also $\max_k|\widehat h_k|$. Für die Umkehrung ist $\|C^{-1}\|_2=1/\min_k|\widehat h_k|$. Ihr Produkt ergibt die Formel oben. Mathematisch gilt das hier, weil $F/\sqrt N$ Längen erhält: Die Singularwerte von $C$ sind genau $|\widehat h_k|$.

Die Faltungsmatrix $C$ und die Entzerrungsmatrix $C^{-1}$ haben dieselbe Kondition: Beim Invertieren werden alle Streckungsfaktoren durch ihre Kehrwerte ersetzt; das Verhältnis von größtem zu kleinstem bleibt gleich.

| Kondition | Bedeutung |
|---|---|
| $\kappa=1$ | Bestmögliche Kondition: keine Verstärkung relativer Fehler in der 2-Norm |
| $\kappa\gg1$ | Relative Fehler können stark verstärkt werden |
| $\kappa=\infty$ | Ein Fourierwert ist null: Information fehlt, keine eindeutige Umkehrung |

**5 · Beispiel: $h=(1,a)$**

Jetzt verändern wir den Echo-Faktor $a$ und beobachten, wann das Entzerren schwierig wird. Wir ergänzen $h$ wieder mit Nullen auf Länge $N$. Aus der DFT-Formel folgt:

$$
\widehat h_k=1+a\,e^{-2\pi i k/N},\qquad \widehat h_0=1+a
$$

Warum? Nur die Einträge $h_0=1$ und $h_1=a$ tragen zur Summe bei: $1\cdot e^0+a\cdot e^{-2\pi i k/N}$.

Besonders anschaulich sind zwei Eigenvektoren für **gerades $N$**:

| Signal-Muster | Echo-Rechnung $y_j=x_j+a x_{j-1}$ | Faktor / Eigenwert |
|---|---|---|
| Konstant: $(1,1,1,1,\ldots)$ | $1+a\cdot1=1+a$ an jeder Stelle | $\widehat h_0=1+a$ |
| Abwechselnd: $(1,-1,1,-1,\ldots)$ | Der Vorgänger hat das entgegengesetzte Vorzeichen: $y_j=x_j-a x_j$ | $\widehat h_{N/2}=1-a$ |

Für $-1<a\leq0$ ist $1+a$ der **kleinste** Betrag und $1-a$ der **größte**. Alle anderen Beträge liegen dazwischen. Daher gilt für gerades $N$:

$$
\boxed{\kappa_2(C)=\frac{1-a}{1+a}}.
$$

Beispiel $a=-0{,}9$: Der konstante Anteil wird mit $0{,}1$ multipliziert, der abwechselnde mit $1{,}9$. Das Verhältnis ist $1{,}9/0{,}1=19$.

| $a$ | $0$ | $-0{,}5$ | $-0{,}9$ | $-0{,}99$ | $-1$ |
|---|---|---|---|---|---|
| $\kappa$ | $1$ | $3$ | $19$ | $199$ | $\infty$ |

Wenn $a$ gegen $-1$ geht, wird der konstante Anteil immer stärker abgeschwächt. Beim Entzerren müssen wir durch die immer kleinere Zahl $1+a$ teilen.

**Bei $a=-1$ wird der konstante Anteil vollständig ausgelöscht.** Zum Beispiel ergibt $x=(5,5,5,5)$ dann $y=(0,0,0,0)$, weil stets $5-5=0$ gilt. Auch $x=(2,2,2,2)$ ergibt dasselbe $y$. Aus dem Ausgang lässt sich deshalb nicht mehr eindeutig bestimmen, welches Eingangssignal vorlag. Die Matrix ist singulär; ihre Kondition wird als $\infty$ angegeben.

**6 · Warum stört Rauschen beim Entzerren?**

In der Praxis empfangen wir nicht genau $y=Cx$, sondern $y_{\text{mess}}=Cx+\eta$. Dabei ist $\eta$ der Messfehler bzw. das Rauschen. Beim Dividieren können wir Signal und Fehler nicht trennen:

$$
\frac{\widehat y_k+\widehat\eta_k}{\widehat h_k}
=\widehat x_k+\underbrace{\frac{\widehat\eta_k}{\widehat h_k}}_{\text{mitverstärktes Rauschen}}
$$

Ein konkretes Beispiel für einen einzelnen Fourier-Anteil:

1. Ursprünglicher Wert: $\widehat x_k=2$.
2. Echo-Faktor: $\widehat h_k=0{,}01$. Am Ausgang bleibt $\widehat y_k=0{,}01\cdot2=0{,}02$.
3. Ein Messfehler von $\widehat\eta_k=0{,}001$ verändert den Messwert auf $0{,}021$.
4. Entzerren: $0{,}021/0{,}01=2{,}1$. Statt $2$ erhalten wir $2{,}1$.

Der **absolute** Fehler wurde von $0{,}001$ auf $0{,}1$ vergrößert, also um den Faktor $100$. Die absolute Verstärkung dieser Frequenz ist $1/|\widehat h_k|$. Die Kondition betrachtet dagegen die mögliche Verstärkung des **relativen Gesamtfehlers** über alle Frequenzanteile.

Nennen wir das rekonstruierte Signal $x_{\text{rek}}=C^{-1}y_{\text{mess}}$. Sein Fehler ist $\delta x=x_{\text{rek}}-x=C^{-1}\eta$. Für $x\neq0$ und invertierbares $C$ gilt dann:

$$
\underbrace{\frac{\|\delta x\|_2}{\|x\|_2}}_{\text{relativer Rekonstruktionsfehler}}
\leq\kappa_2(C)\,
\underbrace{\frac{\|\eta\|_2}{\|y\|_2}}_{\text{relativer Messfehler}}
$$

Die Norm $\|\cdot\|_2$ misst hier die Größe des gesamten Vektors (Wurzel aus der Summe der quadrierten Beträge). „Relativ“ heißt: Fehlergröße geteilt durch die Größe des fehlerfreien Signals.

Bei $\kappa=19$ und einem relativen Messfehler von $1\%$ beträgt diese Schranke beispielsweise $19\cdot1\%=19\%$ relativer Rekonstruktionsfehler. Das ist eine **obere Schranke**, kein garantierter Fehler. Wie groß der tatsächliche Fehler ist, hängt auch davon ab, welche Fourier-Anteile Signal und Rauschen besitzen. Hier betrachten wir Fehler in $y$ bei festem, exakt bekanntem $h$.

### Was berechnet unser Programm für Aufgabe 3?

In `main.py` steht $h=(1,\,0{,}5,\,0{,}25,\,-0{,}25)$ und $N=16$. Die Funktion `kondition(h, N)` setzt genau die oben hergeleitete Rechnung um:

1. Für jedes $k=0,\ldots,15$ wird $\widehat h_k=\sum_{n=0}^{3}h_n e^{-2\pi i kn/16}$ berechnet. Die fehlenden Einträge bis Länge $16$ sind null und müssen deshalb nicht mit aufsummiert werden.
2. `abs(h_dach)` bildet den Betrag dieses Fourierwerts; die Liste `betraege` sammelt alle $16$ Beträge.
3. `max(betraege) / min(betraege)` liefert die Kondition, hier ungefähr $3{,}34227$. Ist der kleinste Betrag null, gibt die Funktion stattdessen unendlich zurück.

**Für diese Aufgabe rekonstruieren wir also kein konkretes Signal. Wir untersuchen allein anhand von $h$, wie empfindlich eine spätere Entzerrung auf Messfehler reagieren kann.**
