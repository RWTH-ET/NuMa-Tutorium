# Code-Beschreibung: Implementierung der Interpolationsverfahren (Aufgabe 2)

Dieses Python-Skript setzt die theoretischen Interpolationsverfahren in die Praxis um. Es berechnet die Temperaturkurven für die Tschebyscheff-Stützstellen für die Polynomgrade $n = 2$ bis $40$ und vergleicht die Fehler der beiden Methoden visuell.

Das Skript gliedert sich in vier Hauptbereiche:

## 1. Interpolationsfunktionen aus Task 1 wiederverwenden

Das Skript lädt `../Task 1/main.py` mit `importlib.util` als Modul `task1`.
So können die vorhandenen Funktionen trotz der Leerzeichen im Ordnernamen
direkt verwendet werden. Der geschützte `main()`-Aufruf von Task 1 verhindert,
dass dessen Experimente beim Import gestartet werden.

*   **Aufgabe 1a:** `task1.vandermonde(E_nodes, T_kalib)` liefert die
    Koeffizienten und die Konditionszahl. Hier werden nur die Koeffizienten
    benötigt. `task1.horner(a, E)` wertet das Polynom für jeden Rasterpunkt aus.
*   **Aufgabe 1b:** `task1.baryzentrische_gewichte(E_nodes)` berechnet die Gewichte
    einmal pro Grad. `task1.baryzentrisch(E_nodes, T_kalib, w, E)` wertet das
    Polynom für jeden Rasterpunkt aus, einschließlich Treffern auf Stützstellen.

Task 2 enthält keine eigene Implementierung dieser Verfahren.

## 2. Daten laden und vorbereiten

Mithilfe der Bibliothek `pandas` (`pd.read_csv`) werden die bereits im Ordner
von Task 1 vorhandenen CSV-Dateien eingelesen. Die Pfade werden relativ zur
Skriptdatei bestimmt, unabhängig vom aktuellen Arbeitsordner.
*   Die erste Zeile (Kommentar mit `#`) wird mit `skiprows=1` übersprungen, stattdessen werden die Spaltennamen über `names=[...]` fest im Code zugewiesen.
*   **`raster` (Prüfraster):** Enthält 2000 feingliedrige Referenzpunkte (`T_exakt`), gegen die unsere berechneten Polynome später getestet werden.
*   **`stuetz` (Tschebyscheff-Stützstellen):** Enthält die vom Kalibrator gerundeten Eingabedaten für verschiedene Polynomgrade $n$.

## 3. Experimente durchführen (Die Schleife)

Das Herzstück des Skripts ist eine `for`-Schleife, die alle Polynomgrade von $n = 2$ bis $n = 40$ durchläuft:
1.  **Filtern:** Für den aktuellen Durchlauf werden nur die Stützstellen aus der Tabelle geladen, die zum aktuellen Grad $n$ gehören.
2.  **Berechnen:** Die Funktionen aus Task 1 berechnen mit diesen Stützstellen die Temperaturen für alle 2000 Punkte des Prüfrasters. Ist das Vandermonde-System numerisch singulär, wird dessen Ergebnis als `NaN` gespeichert; die baryzentrische Auswertung läuft weiter.
3.  **Fehlerermittlung:** Der maximale absolute Fehler ($e_{max}$) wird berechnet, indem die Differenz zur exakten Temperatur gebildet wird (`np.abs(p - T_exakt)`). Der größte Wert (`np.max`) wird in den Listen `fehler_vander` bzw. `fehler_bary` für das spätere Diagramm gespeichert.

## 4. Visualisierung der Ergebnisse

Mit `matplotlib.pyplot` werden die gesammelten Fehlerwerte geplottet:
*   `plt.semilogy` sorgt für eine logarithmische y-Achse, da die Fehler bei Vandermonde exponentiell wachsen können.
*   Eine gestrichelte schwarze Linie (`plt.axhline`) markiert das Projektziel von maximal $0,1^{\circ}C$ Abweichung.
*   Zuletzt wird das fertige Diagramm mit `plt.savefig` als hochauflösende PNG-Datei (`aufgabe2_tschebyscheff.png`) im Ordner gespeichert, bevor es auf dem Bildschirm angezeigt wird.
