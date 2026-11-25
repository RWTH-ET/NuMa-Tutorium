"""Teil C, Aufgabe 5A: Konstanter Offset des Kapazitätsmessgeräts.

Start: python main.py (NumPy und Matplotlib erforderlich).
Die Ergebnisse werden neben diesem Skript gespeichert.
"""

from pathlib import Path
import csv
import importlib.util

import numpy as np
import matplotlib.pyplot as plt


ORDNER = Path(__file__).resolve().parent
W_EXAKT = 320.0 * np.log(3.25)  # Mikro-Joule bei v in Volt und C in Mikro-Farad
OFFSETS = (-0.02, 0.0, 0.005, 0.02)  # Mikro-Farad; hier keine zusätzliche Rundung
TOLERANZ = 1.0  # Mikro-Joule


def lade_werkzeuge():
    """Verwendet die Quadratur-Werkzeuge aus Task 1 ohne dessen main() zu starten."""
    datei = ORDNER.parent / "Task 1" / "main.py"
    spec = importlib.util.spec_from_file_location("pp3_teil_a", datei)
    if spec is None or spec.loader is None:
        raise ImportError(f"Quadratur-Werkzeuge konnten nicht geladen werden: {datei}")
    modul = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(modul)
    return modul


def C_exakt(v):
    """Kapazität in Mikro-Farad, Spannung v in Volt."""
    return 10.0 / (1.0 + (v / 8.0) ** 2)


def C_mit_offset(v, delta):
    """Systematischer Messfehler: derselbe Offset delta an jedem Messpunkt."""
    return C_exakt(v) + delta


def analysiere_verfahren(werkzeuge, verfahren, N_werte, offsets=OFFSETS):
    """Erfasst vorzeichenbehaftete Fehler und eine Schranke des Datenfehlers."""
    if verfahren not in ("Gauss-Legendre", "Summierte Trapezregel"):
        raise ValueError(f"Unbekanntes Verfahren: {verfahren}")
    ergebnisse = []
    for N in N_werte:
        N = int(N)
        if verfahren == "Gauss-Legendre":
            x, w = werkzeuge.transformiere(*werkzeuge.gauss_formel(N), 0.0, 12.0)
        else:
            # N verschiedene Messpunkte entsprechen N-1 Teilintervallen.
            if N < 2:
                raise ValueError("Die summierte Trapezregel benötigt N >= 2.")
            x, w = werkzeuge.summierte_formel(
                *werkzeuge.nc_formel(2), 0.0, 12.0, N - 1
            )
        Q_exakt = werkzeuge.quadratur(x * C_exakt(x), w)
        for delta in offsets:
            Q_offset = werkzeuge.quadratur(x * C_mit_offset(x, delta), w)
            ergebnisse.append({
                "verfahren": verfahren,
                "N": N,
                "offset_uF": float(delta),
                "Q_uJ": Q_offset,
                "quadraturfehler_uJ": Q_exakt - W_EXAKT,
                "datenfehler_uJ": Q_offset - Q_exakt,
                "gesamtfehler_uJ": Q_offset - W_EXAKT,
                # Für dieses Modell ist delta * Integral(v) = 72 * delta.
                "theoretischer_datenfehler_uJ": 72.0 * delta,
                "datenfehlerschranke_uJ": abs(delta) * float(np.dot(abs(w), abs(x))),
            })
    return ergebnisse


def speichere_tabelle(ergebnisse):
    ziel = ORDNER / "aufgabe5_offset_ergebnisse.csv"
    with ziel.open("w", newline="", encoding="utf-8-sig") as datei:
        writer = csv.DictWriter(datei, fieldnames=list(ergebnisse[0]))
        writer.writeheader()
        writer.writerows(ergebnisse)
    return ziel


def zeichne_experiment(ergebnisse):
    fig, axes = plt.subplots(2, 2, figsize=(13, 9), sharex="col")
    farben = ("tab:purple", "tab:blue", "tab:green", "tab:orange")
    for spalte, verfahren in enumerate(("Gauss-Legendre", "Summierte Trapezregel")):
        oben, unten = axes[0, spalte], axes[1, spalte]
        for delta, farbe in zip(OFFSETS, farben):
            daten = [r for r in ergebnisse
                     if r["verfahren"] == verfahren and r["offset_uF"] == delta]
            oben.semilogy(
                [r["N"] for r in daten],
                [max(abs(r["gesamtfehler_uJ"]), 1e-13) for r in daten],
                color=farbe, label=fr"$\delta={delta:g}\,\mu\mathrm{{F}}$",
            )
        oben.axhline(TOLERANZ, color="black", ls=":", label="Ziel: Fehler < 1 µJ")
        oben.set_title(verfahren)
        oben.set_ylabel("Absoluter Gesamtfehler [µJ]")
        oben.legend(fontsize=9)

        daten = [r for r in ergebnisse
                 if r["verfahren"] == verfahren and r["offset_uF"] == 0.02]
        for key, label, stil in (
            ("quadraturfehler_uJ", "Quadraturfehler", "--"),
            ("datenfehler_uJ", "Datenfehler", ":"),
            ("gesamtfehler_uJ", "Gesamtfehler", "-"),
        ):
            unten.plot([r["N"] for r in daten], [r[key] for r in daten], stil, label=label)
        unten.set_yscale("symlog", linthresh=0.01)
        unten.axhline(0.0, color="gray", lw=0.7)
        unten.set_title("Fehler mit Vorzeichen für δ = 0,02 µF")
        unten.set_xlabel("Anzahl verschiedener Messpunkte N")
        unten.set_ylabel("Fehler [µJ], symmetrisch logarithmisch")
        unten.legend(fontsize=9)
        oben.grid(True, which="both", ls=":", alpha=0.5)
        unten.grid(True, which="both", ls=":", alpha=0.5)

    fig.suptitle("Aufgabe 5A: Mehr Messpunkte beseitigen keinen konstanten Offset")
    fig.tight_layout()
    ziel = ORDNER / "aufgabe5_offset_diagramm.png"
    fig.savefig(ziel, dpi=200, bbox_inches="tight")
    return ziel


def drucke_auswertung(ergebnisse):
    print("=== AUFGABE 5A: KONSTANTER MESSGERÄTE-OFFSET ===")
    print(f"Referenzenergie W = {W_EXAKT:.9f} µJ")
    print("Messmodell: C_mess(v) = C(v) + delta, ohne zusätzliche Rundung.")
    print("Exakter systematischer Beitrag: E_D = 72 * delta [µJ].")
    print(f"Für |E_D| < 1 µJ muss |delta| < {1 / 72:.8f} µF sein.")
    print("Dies allein garantiert noch nicht, dass der Gesamtfehler < 1 µJ ist.\n")
    for verfahren in ("Gauss-Legendre", "Summierte Trapezregel"):
        print(verfahren)
        for delta in OFFSETS:
            daten = [r for r in ergebnisse
                     if r["verfahren"] == verfahren and r["offset_uF"] == delta]
            erster = next((r for r in daten if abs(r["gesamtfehler_uJ"]) < TOLERANZ), None)
            treffer = f"erstes N mit Fehler < 1 µJ: {erster['N']}" if erster else "kein Treffer im Testbereich"
            print(f"  delta={delta:+.3f} µF: {treffer}; "
                  f"Fehler bei N={daten[-1]['N']}: {abs(daten[-1]['gesamtfehler_uJ']):.6f} µJ")
        print()
    print("Einzelne gute Ergebnisse können durch Aufhebung der Fehler entstehen.")
    print("Die Fehler addieren sich mit Vorzeichen, ihre Beträge im Allgemeinen nicht.")
    print("Ein größerer Messumfang korrigiert den Offset nicht; dafür ist eine Kalibrierung nötig.")


def main():
    werkzeuge = lade_werkzeuge()
    ergebnisse = analysiere_verfahren(werkzeuge, "Gauss-Legendre", range(1, 41))
    ergebnisse += analysiere_verfahren(werkzeuge, "Summierte Trapezregel", range(2, 41))
    tabelle = speichere_tabelle(ergebnisse)
    grafik = zeichne_experiment(ergebnisse)
    drucke_auswertung(ergebnisse)
    print(f"Tabelle: {tabelle}")
    print(f"Grafik: {grafik}")
    plt.show()


if __name__ == "__main__":
    main()
