# Code-Beschreibung: Implementierung der Interpolationsverfahren (Aufgabe 2)

Dieses Python-Skript setzt die theoretischen Interpolationsverfahren in die Praxis um. Es berechnet die Temperaturkurven für die Tschebyscheff-Stützstellen für die Polynomgrade $n = 2$ bis $40$ und vergleicht die Fehler der beiden Methoden visuell.

Das Skript gliedert sich in vier Hauptbereiche:

## 1. Definition der Interpolationsfunktionen

Das Skript definiert zwei Funktionen, die jeweils eine Liste von $x$-Werten (Spannung) und den zugehörigen $y$-Werten (Temperatur) aufnehmen und die berechneten Temperaturen für das Prüfraster (`x_eval`) zurückgeben.

*   **`vandermonde_horner` (Aufgabe 1a):**
    *   Erstellt mit `np.vander` automatisch die Vandermonde-Matrix.
    *   Löst das lineare Gleichungssystem ($V \cdot a = y$) mithilfe von `np.linalg.solve`, um den Koeffizientenvektor `a` zu erhalten.
    *   Wertet das Polynom effizient über eine `for`-Schleife aus, die das Horner-Schema anwendet (Rückwärts-Iterieren durch die Koeffizienten).
*   **`baryzentrisch` (Aufgabe 1b):**
    *   Verzichtet auf Matrix-Operationen. Stattdessen berechnet eine doppelte `for`-Schleife zunächst die festen Gewichte $w_j$.
    *   In einer zweiten Schleife wird für jeden gesuchten Spannungswert im Prüfraster die baryzentrische Formel angewendet. Brüche werden summiert und geteilt, um direkt die Temperatur auszurechnen.

## 2. Daten laden und vorbereiten

Mithilfe der Bibliothek `pandas` (`pd.read_csv`) werden die externen Daten eingelesen. 
*   Die erste Zeile (Kommentar mit `#`) wird mit `skiprows=1` übersprungen, stattdessen werden die Spaltennamen über `names=[...]` fest im Code zugewiesen.
*   **`raster` (Prüfraster):** Enthält 2000 feingliedrige Referenzpunkte (`T_exakt`), gegen die unsere berechneten Polynome später getestet werden.
*   **`stuetz` (Tschebyscheff-Stützstellen):** Enthält die vom Kalibrator gerundeten Eingabedaten für verschiedene Polynomgrade $n$.

## 3. Experimente durchführen (Die Schleife)

Das Herzstück des Skripts ist eine `for`-Schleife, die alle Polynomgrade von $n = 2$ bis $n = 40$ durchläuft:
1.  **Filtern:** Für den aktuellen Durchlauf werden nur die Stützstellen aus der Tabelle geladen, die zum aktuellen Grad $n$ gehören.
2.  **Berechnen:** Beide Funktionen (`vandermonde_horner` und `baryzentrisch`) werden mit diesen Stützstellen gefüttert. Sie berechnen die Temperaturen für alle 2000 Punkte des Prüfrasters.
3.  **Fehlerermittlung:** Der maximale absolute Fehler ($e_{max}$) wird berechnet, indem die Differenz zur exakten Temperatur gebildet wird (`np.abs(p - T_exakt)`). Der größte Wert (`np.max`) wird in den Listen `fehler_vander` bzw. `fehler_bary` für das spätere Diagramm gespeichert.

## 4. Visualisierung der Ergebnisse

Mit `matplotlib.pyplot` werden die gesammelten Fehlerwerte geplottet:
*   `plt.semilogy` sorgt für eine logarithmische y-Achse, da die Fehler bei Vandermonde exponentiell wachsen können.
*   Eine gestrichelte schwarze Linie (`plt.axhline`) markiert das Projektziel von maximal $0,1^{\circ}C$ Abweichung.
*   Zuletzt wird das fertige Diagramm mit `plt.savefig` als hochauflösende PNG-Datei (`aufgabe2_tschebyscheff.png`) im Ordner gespeichert, bevor es auf dem Bildschirm angezeigt wird.