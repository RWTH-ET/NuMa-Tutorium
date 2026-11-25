# Aufgabe 2 – DFT und FFT vergleichen

Die **diskrete Fourier-Transformation** (DFT) zerlegt $N$ Signalwerte $x_j$ in $N$ Frequenzanteile $X_k$. Die FFT ist ein schneller Rechenweg für dieselbe Transformation.

## Direkte DFT

$$
X_k=\sum_{j=0}^{N-1}x_j\,e^{-2\pi i kj/N},\qquad k=0,\ldots,N-1.
$$

Für jeden Frequenzindex $k$ gewichten und addieren wir alle $N$ Eingangswerte. Das Programm schreibt die Faktoren in eine Matrix und multipliziert sie mit dem Signalvektor. $N$ Summen mit je $N$ Termen ergeben den Aufwand $O(N^2)$.

## FFT

Wir teilen das Signal nach geraden und ungeraden Positionen auf, bis die Teilvektoren nur noch einen Wert haben. Sind $A_k$ und $B_k$ die Fourierwerte der geraden und ungeraden Hälfte, setzen wir sie für $0\le k<N/2$ so zusammen.

$$
X_k=A_k+\omega_N^k B_k,\qquad
X_{k+N/2}=A_k-\omega_N^k B_k,\qquad
\omega_N=e^{-2\pi i/N}.
$$

Beide Ergebnisse nutzen dieselben Teilrechnungen. Bei Zweierpotenzen gibt es $\log_2 N$ Ebenen mit je ungefähr $N$ Rechenschritten. So ergibt sich $O(N\log N)$.

## Laufzeitdiagramm

Auf der x-Achse steht $N$, auf der y-Achse die Rechenzeit. Auf beiden Achsen ist die Skala logarithmisch. Eine steilere Kurve bedeutet, dass die Zeit mit $N$ schneller wächst. Mit $1/\text{Zeit pro Transformation}$ schätzen wir die Symbole pro Sekunde. Die Ziele im Programm sind 14.000 für LTE bei $N=2048$ und 28.000 für 5G NR bei $N=4096$. `fast.py` wiederholt die Messung mit der optimierten NumPy-FFT.
