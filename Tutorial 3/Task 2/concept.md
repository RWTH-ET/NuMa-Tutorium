Quadratur-Werkzeuge (Teil A)
gauss_formel(N): Berechnet Stützstellen und Gewichte auf $[-1, 1]$ über den Golub-Welsch-Algorithmus (Eigenwertproblem der tridiagonalen Jacobi-Matrix).
nc_formel(N): Berechnet Newton-Cotes-Gewichte für äquidistante Stützstellen durch Integration der Lagrange-Basispolynome (exakt ausgewertet via höhergradiger Gauß-Quadratur).
transformiere() & quadratur(): Skaliert Stützstellen/Gewichte auf das Intervall $[0, 12]$ und führt die gewichtete Summation $\sum w_i y_i$ aus.

Modell & Messdaten (Teil B)
f_exakt() vs. f_gemessen(): Definiert den exakten Verlauf $f(v) = v \cdot C(v)$ und eine auf 2 Nachkommastellen gerundete Variante zur Simulation von Messdatenfehlern.
W_EXAKT: Der analytisch berechnete Referenzwert ($320 \cdot \ln(3{,}25)$).

Fehleranalyse (analisiere_verfahren)
Evaluiert beide Verfahren für Stützstellenanzahlen $N \in [1, 15]$ und trennt den Gesamtfehler auf in:
Verfahrensfehler: $\vert{}Q(f_{\text{exakt}}) - W_{\text{exakt}}\vert{}$ (Diskretisierungsfehler des Quadraturverfahrens).
Datenfehler: $\vert{}Q(f_{\text{gemessen}}) - Q(f_{\text{exakt}})\vert{}$ (Auswirkung des Mess-/Rundungsfehlers).
Gesamtfehler: $\vert{}Q(f_{\text{gemessen}}) - W_{\text{exakt}}\vert{}$.

Visualisierung & Auswertung
Trägt die Fehlerverläufe logarithmisch über $N$ auf und vergleicht sie mit einer Ziel-Toleranz von $1\,\mu\text{J}$.
Speichert die Grafiken mithilfe von pathlib im Ordner Task 2/aufgabe4_diagramm.png.
Ermittelt und druckt die minimale Anzahl an Messpunkten $N$, die benötigt wird, um die Ziel-Toleranz einzuhalten.