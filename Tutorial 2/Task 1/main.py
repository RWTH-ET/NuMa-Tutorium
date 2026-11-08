from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np


ORDNER = Path(__file__).resolve().parent
ZIEL = 0.1  # Maximaler Temperaturfehler in °C.


def vandermonde(E_stuetz, T):
    V = []
    for E_i in E_stuetz: #Für jede Ref.Spannung E_i
        zeile = []
        for j in range(len(E_stuetz)): #Erstellt Polynom mit Grad len(E_stuetz) - 1
            potenz = float(E_i) ** j
            zeile.append(potenz)
        V.append(zeile)

    a = np.linalg.solve(V, T)
    kappa = np.linalg.cond(V, np.inf) #np.inf für unendlichnorm (größte zeilensumme der Beträge)
    return a, kappa


def horner(a, E):
    #Optimierung
    p = 0.0
    for a_j in reversed(a): #in die Form a_0+E(a_1+Ea_2)
        p = p * E
        p = p + a_j
    return p


def baryzentrische_gewichte(E_stuetz):
    n = len(E_stuetz) - 1
    w = np.zeros(n + 1)

    # w_j = 1 / Produkt(E_j - E_k), fuer alle k ungleich j.
    for j in range(n + 1):
        E_j = E_stuetz[j]
        produkt = 1.0
        for k in range(n + 1):
            if j != k:
                E_k = E_stuetz[k]
                produkt = produkt * (E_j - E_k)

        w[j] = 1.0 / produkt
    return w


def baryzentrisch(E_stuetz, T, w, E):
    """Berechne die Temperatur p(E) fuer eine einzelne Spannung E."""
    zaehler = 0.0
    nenner = 0.0

    for j in range(len(E_stuetz)):
        E_j = E_stuetz[j]
        T_j = T[j]

        # Fuer eine Referenzspannung ist die Temperatur bereits bekannt.
        if E == E_j:
            return T_j

        w_j = w[j]
        b_j = w_j / (E - E_j)

        # Zaehler = b_0*T_0 + b_1*T_1 + ... + b_n*T_n
        beitrag = b_j * T_j
        zaehler = zaehler + beitrag

        # Nenner = b_0 + b_1 + ... + b_n
        nenner = nenner + b_j

    p = zaehler / nenner
    return p


#---------------------------------
# Ergbenisse Speichern / Plotten
#---------------------------------

def speichere_diagramm(dateiname, titel, ylabel, xlabel="Polynomgrad n"):
    """Beschrifte und speichere das aktuelle Diagramm."""
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    plt.title(titel)
    plt.grid(True, which="both", alpha=0.3)
    plt.legend()
    plt.tight_layout()
    plt.savefig(ORDNER / dateiname, dpi=160)


def plotte_ergebnisse(n_werte, e_max_v, e_max_b, kappa_werte, E_raster,
                     abweichung_v, abweichung_b):
    """Zeige je zwei Diagramme fuer Aufgabe 1a und 1b."""
    fehlerachse = "Maximaler absoluter Fehler [°C]"

    plt.figure(figsize=(9, 5))
    plt.semilogy(n_werte, e_max_v, "o-", label="Vandermonde + Horner")
    plt.axhline(ZIEL, color="black", linestyle="--", label="Ziel: 0,1 °C")
    speichere_diagramm("1a_fehler_vandermonde.png",
                       "Aufgabe 1a: Fehler von Vandermonde + Horner", fehlerachse)

    plt.figure(figsize=(9, 5))
    plt.semilogy(n_werte, kappa_werte, "o-", label="Vandermonde-Matrix")
    speichere_diagramm("1a_kondition.png",
                       "Aufgabe 1a: Kondition der Vandermonde-Matrix",
                       r"Konditionszahl $\kappa_\infty(V)$")

    plt.figure(figsize=(9, 5))
    plt.semilogy(n_werte, e_max_v, "o-", label="Vandermonde + Horner (1a)")
    plt.semilogy(n_werte, e_max_b, "s-", label="Baryzentrisch (1b)")
    plt.axhline(ZIEL, color="black", linestyle="--", label="Ziel: 0,1 °C")
    speichere_diagramm("1b_fehlervergleich.png",
                       "Aufgabe 1b: Fehlervergleich beider Verfahren", fehlerachse)

    plt.figure(figsize=(9, 5))
    plt.plot(E_raster, abweichung_v, label="Vandermonde + Horner (1a)")
    plt.plot(E_raster, abweichung_b, label="Baryzentrisch (1b)")
    # Symlog erlaubt negative Abweichungen und null.
    plt.yscale("symlog", linthresh=0.005)
    plt.axhline(0, color="black", linewidth=0.8)
    plt.axhline(ZIEL, color="gray", linestyle="--", label="±0,1 °C")
    plt.axhline(-ZIEL, color="gray", linestyle="--")
    speichere_diagramm("1b_abweichung_n20.png", "Aufgabe 1b: Abweichungen bei n = 20",
                       "Abweichung p(E) - T_exakt(E) [°C] (symlog)", "Thermospannung E [mV]")
    plt.show()


def zeige_auswertung(name, n_werte, e_max, abweichung, E_raster):
    passende_grade = []
    for i in range(len(n_werte)):
        n = n_werte[i]
        if e_max[i] < ZIEL:
            passende_grade.append(n)
    print(f"\n{name}: Grade mit Fehler < {ZIEL} C: {passende_grade}")

    # Ein fehlgeschlagenes Vandermonde-System liefert NaN statt einer Zahl.
    endliche_werte = np.isfinite(abweichung)
    if np.all(endliche_werte):
        betraege = np.abs(abweichung)
        i = np.argmax(betraege)
        print(f"Bei n = 20: groesster Fehler bei E = {E_raster[i]:.6f} mV, "
              f"Abweichung = {abweichung[i]:.6f} C")


def main():
    # Die CSV-Dateien enthalten Leerzeichen als Trenner und eine #-Kopfzeile.
    pruefraster = np.loadtxt(ORDNER / "pruefraster.csv")
    E_raster = []
    T_exakt = []
    for zeile in pruefraster:
        E_raster.append(zeile[0])
        T_exakt.append(zeile[1])

    # NumPy-Arrays erlauben spaeter das Berechnen der Fehler fuer alle Punkte.
    E_raster = np.array(E_raster)
    T_exakt = np.array(T_exakt)
    stuetzwerte = np.loadtxt(ORDNER / "stuetzwerte_aequidistant.csv")

    n_werte = range(2, 41)
    e_max_v = []
    e_max_b = []
    kappa_werte = []

    print(f"{'n':>3} {'Fehler Vandermonde [C]':>24} "
          f"{'Fehler baryzentrisch [C]':>26} {'Kondition (inf)':>17}")

    for n in n_werte:
        # Nur die Referenzpunkte fuer den aktuellen Grad n auswaehlen.
        E_stuetz = []
        T = []
        for zeile in stuetzwerte:
            if zeile[0] == n:
                E_stuetz.append(zeile[2])  # Referenzspannung E_i.
                T.append(zeile[4])  # Kalibratorwert T_i.

        E_stuetz = np.array(E_stuetz)
        T = np.array(T)
        if len(E_stuetz) != n + 1 or len(np.unique(E_stuetz)) != n + 1:
            raise ValueError(
                f"Fuer n = {n} werden {n + 1} verschiedene Stuetzstellen benoetigt.")

        try:
            a, kappa = vandermonde(E_stuetz, T)
            p_v = np.zeros(len(E_raster))
            for i in range(len(E_raster)):
                E = E_raster[i]
                p_v[i] = horner(a, E)
        except np.linalg.LinAlgError:
            print(f"n = {n}: Vandermonde-System numerisch singulaer.")
            p_v = np.full_like(E_raster, np.nan)
            kappa = np.inf

        w = baryzentrische_gewichte(E_stuetz)
        p_b = np.zeros(len(E_raster))
        for i in range(len(E_raster)):
            E = E_raster[i]
            p_b[i] = baryzentrisch(E_stuetz, T, w, E)

        # NumPy fuehrt jede dieser Rechnungen fuer alle Rasterpunkte aus.
        differenz_v = p_v - T_exakt
        differenz_b = p_b - T_exakt
        betraege_v = np.abs(differenz_v)
        betraege_b = np.abs(differenz_b)
        e_max_v_n = np.max(betraege_v)
        e_max_b_n = np.max(betraege_b)

        e_max_v.append(e_max_v_n)
        e_max_b.append(e_max_b_n)
        kappa_werte.append(kappa)
        print(f"{n:3d} {e_max_v_n:24.6e} {e_max_b_n:26.6e} {kappa:17.6e}")

        if n == 20:
            abweichung_v = differenz_v
            abweichung_b = differenz_b

    zeige_auswertung("Vandermonde", n_werte, e_max_v, abweichung_v, E_raster)
    zeige_auswertung("Baryzentrisch", n_werte, e_max_b, abweichung_b, E_raster)

    plotte_ergebnisse(n_werte, e_max_v, e_max_b, kappa_werte, E_raster,
                     abweichung_v, abweichung_b)


if __name__ == "__main__":
    main()
