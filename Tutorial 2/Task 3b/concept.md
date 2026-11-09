# Aufgabe 3: Woher kommt der Fehler?

Wir rechnen eine Thermospannung in eine Temperatur um. Aufgabe 3 untersucht,
warum das Ergebnis von der exakten Temperatur abweicht:

| Teil | Ursache | Was wird verglichen? |
|---|---|---|
| **3a: Interpolationsfehler** | Ein Polynom bildet die Kennlinie nur näherungsweise ab. | Polynom aus exakten Stützwerten mit der exakten Kennlinie |
| **3b: Datenfehler** | Der Kalibrator liefert gerundete Temperaturen. | Polynom aus gerundeten mit dem Polynom aus exakten Stützwerten |
| **3c: Koeffizientenfehler** | Bei der Koeffizientenberechnung und Auswertung entstehen Rundungsfehler im Rechner. | Vandermonde + Horner mit der baryzentrischen Referenz |

Dieser Ordner implementiert **3b**. Den Interpolationsfehler aus 3a berechnen
wir für den geforderten Vergleich mit.

## Zwei Polynome vergleichen

Der Kalibrator rundet auf $0{,}01\,°C$. Ein Stützwert weicht deshalb höchstens
um $0{,}005\,°C$ vom exakten Wert ab. Wir verwenden dieselben Spannungen,
aber zwei unterschiedliche Temperaturspalten:

- `T_exakt_C` ergibt das Polynom $p_{\mathrm{exakt}}$.
- `T_kalibrator_C` ergibt das Polynom $p_{\mathrm{kalibrator}}$.

Beide Polynome werden baryzentrisch ausgewertet, wie in Task 1 und 2.
Für jede Spannung im Prüfraster berechnen wir die Differenz. Der größte
Betrag ist der Datenfehler:

$$
e_{\mathrm{daten}}=\max_{E\in\mathrm{Prüfraster}}
|p_{\mathrm{kalibrator}}(E)-p_{\mathrm{exakt}}(E)|.
$$

Der Interpolationsfehler vergleicht dagegen $p_{\mathrm{exakt}}(E)$ mit
der exakten Temperatur aus `pruefraster.csv`.

## Wie stark werden die Rundungsfehler verstärkt?

Aus den baryzentrischen Gewichten berechnen wir die Lagrange-Basiswerte:

$$
\ell_j(E)=\frac{w_j/(E-E_j)}{\sum_k w_k/(E-E_k)}.
$$

Damit gilt $p(E)=\sum_j\ell_j(E)T_j$. Die Summe der Beträge dieser
Basiswerte beschreibt, wie stark sich Datenfehler im ungünstigsten Fall
auswirken. Ihr Maximum ist die Lebesgue-Konstante $\Lambda_n$:

$$
e_{\mathrm{daten}}\leq\Lambda_n\cdot0{,}005\,°C.
$$

**Beispiel:** Bei $\Lambda_n=3$ beträgt die Schranke $0{,}015\,°C$.
Der tatsächliche Datenfehler kann deutlich kleiner sein.

Im Code bestimmen wir das Maximum über die 2000 Punkte des Prüfrasters.
Die berechnete Schranke gilt damit für dieses Raster; sie ist kein Nachweis
für das ganze kontinuierliche Intervall. Trifft eine Spannung genau eine
Stützstelle, setzen wir deren Basiswert auf 1 und die übrigen auf 0.
So vermeiden wir eine Division durch null.

## Was zeigen die Ergebnisse?

- **Kleine Grade:** Der Interpolationsfehler dominiert. Das Polynom bildet
  die Kennlinie noch nicht ausreichend genau ab.
- **Große Grade, äquidistant:** Der Datenfehler wächst stark, weil die
  Stützstellenwahl die kleinen Rundungsfehler stark verstärkt. Das erklärt
  die großen baryzentrischen Fehler aus Aufgabe 1. Ein gelegentlich kleinerer
  Vandermonde-Fehler kann durch teilweise Aufhebung von Daten- und Rechenfehlern
  entstehen und ist kein verlässlicher Genauigkeitsvorteil.
- **Große Grade, Tschebyscheff:** Der Interpolationsfehler sinkt stark, der
  Datenfehler bleibt im Bereich weniger Tausendstel Grad und dominiert schließlich.
  Mehr Stützstellen machen die Kalibratorwerte nicht genauer.

Mit den vorliegenden gerundeten Daten lässt sich der Fehler also nicht
allein durch Erhöhen des Grades beliebig klein machen. Die Schranke selbst
ist dabei eine Obergrenze, keine feste Untergrenze für den tatsächlichen Fehler.

## Ausführen und Ausgaben

Im Ordner `Task 3b`:

```sh
python -m pip install -r requirements.txt
python main.py
```

Für beide Stützstellenwahlen werden die Grade $n=2,\ldots,40$ untersucht.
Alle Dateipfade beziehen sich auf den Skriptordner. Die drei Eingabedateien
sind unveränderte Kopien der Messdaten aus Task 1.

- `3b_aequidistant.png` und `3b_tschebyscheff.png` zeigen jeweils
  Interpolationsfehler, Datenfehler und Schranke mit logarithmischer Fehlerachse.
- `ergebnisse_aequidistant.csv` und `ergebnisse_tschebyscheff.csv` enthalten
  die zugehörigen Zahlen einschließlich der Lebesgue-Konstante auf dem Raster.

Bei hohen äquidistanten Graden kann auch die baryzentrische Rechnung durch
Rundungsfehler im Rechner beeinflusst werden. Die Ergebnisse sind daher
numerische Messwerte, keine exakten Fehlerwerte.
