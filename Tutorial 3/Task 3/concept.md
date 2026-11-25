# Aufgabe 5A – Was passiert bei einem festen Messfehler?

## Versuchsaufbau

Ein Gerät mit Offset zeigt bei jeder Messung denselben Betrag zu viel oder zu wenig an. Wir nutzen die Kapazitätskurve aus [Aufgabe 4](../Task%202/concept.md) und setzen $C_\mathrm{mess}(v)=C(v)+\delta$. Wir testen $\delta=-0{,}02$, $0$, $0{,}005$ und $0{,}02\,\mu\mathrm{F}$ mit Gauß-Legendre und summierter Trapezregel. $N$ ist in beiden Fällen die Anzahl verschiedener Messpunkte. Bei der Trapezregel ergeben sich daraus $N-1$ Teilintervalle.

## Warum bleibt der Fehler bestehen?

Durch den Offset ändert sich der Integrand von $vC(v)$ zu $vC(v)+\delta v$. Beide Verfahren integrieren den linearen Zusatz $\delta v$ exakt, bis auf die Rundung des Rechners.

$$
E_D=Q_\mathrm{mess}-Q_\mathrm{exakt}
=\delta\!\int_0^{12}v\,dv
=\delta \left[\frac{v^2}{2}\right]_0^{12}
=72\delta.
$$

Bei $+0{,}02\,\mu\mathrm{F}$ erhöht sich die berechnete Energie dadurch um $1{,}44\,\mu\mathrm{J}$. Bei $+0{,}005\,\mu\mathrm{F}$ sind es $0{,}36\,\mu\mathrm{J}$. Für $|E_D|<1\,\mu\mathrm{J}$ muss der Betrag des Offsets unter $1/72\,\mu\mathrm{F}$ liegen, also ungefähr unter $0{,}0139\,\mu\mathrm{F}$.

## Wie liest man das Diagramm?

Oben steht der Betrag des Gesamtfehlers. Unten sieht man für $+0{,}02\,\mu\mathrm{F}$, wie sich der Verfahrensfehler $E_Q=Q_\mathrm{exakt}-W$ und der Datenfehler $E_D$ mit Vorzeichen zu $E_\mathrm{ges}=E_Q+E_D$ addieren.

Bei zwei Gauß-Punkten ergibt sich etwa $-2{,}26+1{,}44=-0{,}82\,\mu\mathrm{J}$. Damit liegt der Betrag kurz unter dem Ziel von $1\,\mu\mathrm{J}$. Bei 40 Punkten ist der Verfahrensfehler fast weg und der Gesamtfehler liegt bei $1{,}44\,\mu\mathrm{J}$. Mehr Punkte machen die Quadratur genauer. Der feste Offset bleibt und kann das Ergebnis über die Toleranzgrenze schieben. [Diagramm](aufgabe5_offset_diagramm.png) · [Messwerte](aufgabe5_offset_ergebnisse.csv)
