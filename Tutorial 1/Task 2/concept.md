Algorithmen:

*dft(x): (Diskrete Fourier-Transformation)*
-Setzt die mathematische Formel direkt um
-Nutzt eine Matrix aus komplexen Exponenten (numpy) anstatt zwei verschachtelte Schleifen
-Multipliziert die Matrix via Skalarprodukt (ny.dot) mit dem Eingangsvektor
*Laufzeit:*   Quadratisch, da für N Elemente N berechnungen stattfinen O(N^2)


*fft(x): (Schnelle Fourier-Transformation)*
-Funktioniert nach dem Cooley-Tukey-Algorithmus (Teile und Herrsche) (Rekursiv)
-Array wird in Elemente mit geraden (x[0::2]) und ungeraden Indizes (x[1::2]) aufgeteilt
-Wenn der Vektor nur noch ein Element lang ist wird er zurückgegeben
-Hälften werden rekursiv an die Funktion selbst übergeben
-Das Aufteilen verursacht eine Phasenverschiebung im Frequenzbereich, wird durch den 
 komplexen Rotationsfaktor korrigiert
-Zum Schluss werden sie übr den komplexen Drehfaktor T wieder zusammengesetzt (Symmetrie)
*Laufzeit:*   NlogN, dank rekursive Implementierung und Symmetrie (O(NlogN))