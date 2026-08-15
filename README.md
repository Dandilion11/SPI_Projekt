# SPI-Projekt: Struktur- und Parameteridentifikation

Modellierung eines Hefe-Fermentationsprozesses (Krämer & King 2017) mit drei
Modellen unterschiedlicher Komplexität, Parameteridentifikation,
Unsicherheitsanalyse über die Fisher-Informationsmatrix und optimale
Versuchsplanung.

---

## Was dieser Branch ändert

Kurzfassung der wichtigsten Änderungen gegenüber dem vorherigen Stand.
Details stehen jeweils im Kopf der betroffenen Datei.

### Datenaufbereitung

| Änderung | Warum |
|---|---|
| `sigma^2 = a*y + b` statt `(a*y + b)^2` | Das Varianzmodell war fälschlich quadriert. |
| Einheitenumrechnung korrigiert: Base `/1e6`, O2 `*1e4` | Varianz skaliert **quadratisch** mit einer Einheitenumrechnung, nicht linear. |
| Messreihen werden nach Zeit sortiert, NaN entfernt | `x0` nimmt `y(1)` — das muss der erste Wert *nach der Zeit* sein, nicht der erste im Speicher. |
| Kriterien (c) und (d) systematisch geprüft | Ergebnis: alle Offline-Messzeiten liegen exakt auf Probenahmezeitpunkten; beide auffälligen Glucose-Anstiege in Def03 sind bilanzkonform. Es bleibt **eine** Bereinigung im ganzen Projekt (O2-Sensorausfall in Def03). |
| Varianz aus `feld.VarParam.a/b` statt aus dem Paper |  |
### Simulation

| Änderung | Warum |
|---|---|
| Segmentgrenzen nur an echten Unstetigkeiten (Feed, Probenahme), nicht an jeder Auswertezeit | Vorher startete `ode15s` an jedem Messpunkt neu und das Ergebnis hing vom Auswerteraster ab. Gitterabhängigkeit: 3.6e-01 → 1e-08. |
| `x = max(x,0)` aus allen Modelldateien entfernt | Nicht differenzierbar und inkonsistent mit den Jacobi-Matrizen. Positivität sichert stattdessen `NonNegative` im Integrator. |
| Deterministisches Schrittlimit (`ode_steplimit`) | Pathologische Parametersätze ließen den Optimierer minutenlang hängen. Ein Schrittzähler ist reproduzierbar, ein Zeitlimit wäre es nicht. |

### Gütefunktional

| Änderung | Warum |
|---|---|
| Mittelung pro Messkanal | DOT und Base haben hunderte Punkte, Biomasse ~20. Ohne Mittelung dominieren die dicht abgetasteten Kanäle allein durch ihre Punktzahl. Dadurch ist `N_eff` = Anzahl aktiver Kanäle und `chi2/N_eff ≈ 1` die Erwartung bei korrektem Modell. |
| Base als **Inkrement** bewertet | `m_B` ist ein kumulatives Integral: benachbarte Punkte teilen ihre Integrationsgeschichte, ihre Fehler sind autokorreliert. `sigma` beschreibt nur das Rauschen *einer* Zugabe. Vorher trug Base bis zu 90 % von J und verbog `mumax`, `KS` und `YXS`. |
| Handgewichte auf DOT (`/20`, `/10`) entfernt | Willkürlich, in beiden Modellen unterschiedlich, und dadurch nicht vergleichbar. Robustheitstest: die identifizierten Parameter ändern sich um < 1 %. |
| Trapez-Zeitgewichte getestet und **verworfen** | Parameter verschieben sich um 6–14 %, die Aussagen bleiben gleich, aber der Validierungsfehler ist höher (Val/Train 1.84 statt 1.49). Zeitgewichte betonen isolierte Punkte nach langen Abtastlücken — genau dort, wo die fehlende Ethanol-Kinetik am stärksten wirkt. |

### Identifizierbarkeit — der rote Faden

Alle drei Modelle enthalten Parameterkombinationen, die die Daten festlegen,
während die einzelnen Faktoren offen bleiben. Konsequenz: solche Parameter
werden **fixiert**, gefittet wird nur, was die Daten bestimmen.

| Parameter | Behandlung | Begründung |
|---|---|---|
| `KLa` | fixiert auf 1, `YXO` wird zu `YXO_eff = KLa·YXO` | DOT ist quasistationär (tau ~ Sekunden, Abtastung ~Minuten). Drei Optima mit KLa = 200 / 7.1 / 714 ergaben identisches J bei konstantem Produkt. **Reparametrisierung, keine Annahme.** |
| `KS` (Modell 3) | fixiert auf 4.29 aus dem Batchfit | Nur identifizierbar, wo `cGlc` durch `KS` läuft — also in der Batch-Phase. Im Fed-Batch lief `KS` an die untere Grenze. Nach dem Fixieren stimmt `mumax` (0.343) mit Modell 1 (0.317) überein. |
| `KAm`, `KPh` | fixiert auf 0.01 | Nicht identifizierbar (Sensitivität < 1 %), bleiben aber im Modell: ohne sie laufen `mAm`/`mPh` negativ. |

Beobachtete Korrelationen (aus der FIM):

* Modell 1: `mumax`–`KS` = **0.97**
* Modell 2: `Y_Bam`–`Y_AmX` = **−0.95**
* Modell 3: `YXS`–`YAmX` = **−0.82**

Drei Modelle, drei Produkt-Degenerationen, drei verschiedene Ursachen —
zu schnelle DOT-Dynamik, Ammonium am Rauschgrund, gemeinsamer Faktor
`rX·mX`. Das ist die Motivation für die Versuchsplanung.

### Optimierung

| Änderung | Warum |
|---|---|
| Log-Parametrisierung (`q = log10(p)`) | Die Parameter spannen vier Größenordnungen; linear ist die Hesse-Matrix schlecht konditioniert. Das Optimum bleibt identisch, nur die Konvergenz wird besser. |
| Referenzpunkt `pRef` als Startpunkt **und** Akzeptanzschwelle | Ein früherer Lauf startete nur von `p0` und landete schlechter als ein bereits bekannter Vektor, ohne dass es auffiel. |
| LHS-Multistart | Bestätigt, dass das gefundene Optimum nicht nur ein lokaler Zufallstreffer ist. |

### Fisher-Informationsmatrix

* Jacobi-Matrizen neu erzeugt — die alten waren aus einer **doppelten**
  Monod-Kinetik abgeleitet, der Phosphat-Term fehlte (`dfdp(:,10) = 0`).
* Funktionsnamen an Dateinamen angeglichen: MATLAB wählt beim Aufruf die
  *Datei*. Vorher landete `Modell3_dfdp` still bei der alten 8-Parameter-Version.
* FIM summiert über **alle** Trainingsexperimente, mit Varianz je Messpunkt.
* Probenahmen werden berücksichtigt: dort springt nicht nur `x`, sondern
  über die Kettenregel auch `XP`.
* FIM wird auf die **freien** Parameter eingeschränkt — über alle wäre sie
  singulär.

---

## Ergebnisse im Überblick

| | Modell 1 | Modell 2 | Modell 3 |
|---|---|---|---|
| Typ | Batch | Batch + Am + Base | Fed-Batch |
| freie Parameter | 4 | 6 | 6 |
| chi2/N_eff Training | 0.11 | 1.52 | 44.0 |
| Val/Train | 1.72 | 2.87 | 1.38 |

`mumax` liegt in allen drei Modellen bei 0.32–0.38 1/h, `YXS` bei
0.08–0.14 g/g. Der Wert von ~0.1 g/g gegenüber ~0.5 g/g bei rein
respiratorischem Wachstum ist ein unabhängiger Hinweis auf Overflow-
Metabolismus (Ethanol).

**Zentraler Befund:** Modell 1 und 2 passen die Batch-Phase bis auf
Messgenauigkeit (chi2/N_eff ≈ 0.1), Modell 3 verfehlt den Fed-Batch um
Faktor 44 — bei identischem Varianzmodell. Die Fitgüte korreliert mit dem
Glucose-Regime:

| Experiment | chi2/N | cGlc > 0.5 g/L |
|---|---|---|
| Def03 | 21.4 | 76 % |
| Def04 | 29.2 | 73 % |
| Def06 | 53.7 | 39 % |
| Def07 | 60.7 | 40 % |
| Def10 | 71.9 | 21 % |

Dazu passend: bei Glucose-Überschuss **über**schätzt das Modell die
Biomasse, bei Glucose-Mangel **unter**schätzt es sie. Ein Vorzeichenwechsel
lässt sich nicht durch einen falsch getroffenen Ertragskoeffizienten
erklären, nur durch einen fehlenden Stoffwechselweg. Ergänzend: in Def07
wurden 23.16 g Ammonium zugeführt und 0.00 g blieben übrig — das Modell
kann rund 20 g nicht bilanzieren, weil bei `cGlc = 0` auch `rX = 0` gilt
und damit jede Ammoniumaufnahme aufhört.

---

## Dateien

### Datenaufbereitung

| Datei | Inhalt |
|---|---|
| `main_01_data_preprocessing.m` | Preprocessing Modell 1/2. Schneidet auf die feedfreie Batch-Phase (Def03 ≤ 6 h, Def04 ≤ 8 h), berechnet die Messvarianz, speichert `Processed_Batch_Data.mat`. |
| `main_01_data_preprocessing_m3_extended.m` | Preprocessing Modell 3 über alle Experimente. Ruft `clean_mess`, baut `x0` (Konzentration × V0 → Masse), speichert `TrainSet`/`ValSet`. Enthält Kontrollausgaben und einen Varianz-Plausibilitätscheck. |
| `clean_mess.m` | Die Bereinigungslogik. Eigene Datei, damit Fit und Bereinigungs-Plot garantiert dasselbe verwenden. Entfernt aktuell nur den O2-Sensorausfall in Def03 (t = 18.5–23.6 h). |
| `find_outlier_points.m` | Sucht systematisch Ausreißer-**Kandidaten** nach den Kriterien (a)–(d). Entfernt nichts — liefert die Belege für die Entscheidungen in `clean_mess`. |
| `plot_cleaning.m` | Abbildung für einen Datensatz: blau = behalten, rot = entfernt. |
| `create_cleaned_plots.m` | Erzeugt und speichert diese Abbildungen für alle Datensätze. |

### Modelle

| Datei | Inhalt |
|---|---|
| `Modell1.m` | Batch: `[cX; cGlc; DOT]`. |
| `Modell2.m` | Batch + Ammonium + Base: `[cX; cGlc; cAm; cBase; DOT]`. |
| `Modell3_woEtOH_10p.m` | Fed-Batch: `[V; mX; mGlc; mAm; mPh; mB; DOT]`, dreifache Monod-Kinetik. |
| `Modell3_woEtOH_XP.m` | Erweitertes System: Zustände **und** Sensitivitäten `XP = dx/dtheta`. |
| `Modell3_mgl.m` / `Modell3_dmgldx.m` | Messgleichung (Massen → Konzentrationen) und ihre Jacobi-Matrix. |
| `Modell3_dfdx_10p.m` / `Modell3_dfdp_10p.m` | **Automatisch erzeugt** von `gen_jacobians.m`. Nicht von Hand ändern. |
| `sim_m3_sample_10p.m` | Simulationswrapper: segmentweise Integration mit Probenahme-Sprüngen. |
| `probe_m3.m` | Diskreter Zustandssprung bei einer Probenahme. |

### Identifikation und Analyse

| Datei | Inhalt |
|---|---|
| `main_04_parameter_fitting_LHS.m` | Fit Modell 1 und 2 mit LHS-Multistart. Ein `spec`-Table pro Modell steuert Kanäle, Zustandsindizes und `x0`. |
| `main_04_parameter_fitting_m3_woEtOH_LHS_10p_multi.m` | Fit Modell 3 über vier Experimente gleichzeitig. Enthält zusätzlich die lokale Sensitivitätstabelle, die Auswertung über alle fünf Experimente und den Nachweis der `KLa·YXO`-Degeneration. |
| `main_05_modellvergleich_8p_vs_10p.m` | Vergleicht das alte 8-Parameter-Modell gegen die 10-Parameter-Version bei **gleicher** Zahl freier Parameter. Zeigt: das alte Modell fittet Ammonium nur besser, weil es mehr Ammonium verbraucht als vorhanden ist (`mAm` bis −23 g). |
| `main_05_uncertainty_FIM_m1.m` / `_m2.m` / `_m3_multi.m` | Parameterunsicherheit über die FIM: Standardabweichungen, Korrelationsmatrix, Konditionszahl, Unsicherheits-Ellipsoid. |
| `main_06_ovp_modell3.m` | Optimale Versuchsplanung: optimiert den Glucose-Feedverlauf eines zusätzlichen Experiments nach dem A-Kriterium. |
| `gen_jacobians.m` | Erzeugt `Modell3_dfdx_10p.m` und `Modell3_dfdp_10p.m` symbolisch. |
| `check_jacobians.m` | Prüft sie gegen finite Differenzen. **Nach jedem `gen_jacobians`-Lauf ausführen** — eine falsche Jacobi-Matrix fällt sonst nirgends auf, sie verfälscht nur still die FIM. |

---

## Reihenfolge der Ausführung

```
1. main_01_data_preprocessing.m              -> Processed_Batch_Data.mat
   main_01_data_preprocessing_m3_extended.m  -> Processed_FedBatch_..._MultiExp.mat

2. main_04_parameter_fitting_LHS.m           -> p_opt.mat, p_opt_Modell2.mat
   main_04_parameter_fitting_m3_...multi.m   -> p_opt_Modell3_...multi.mat

3. gen_jacobians.m  ->  check_jacobians.m    (nur bei Modelländerung)

4. main_05_uncertainty_FIM_m1 / _m2 / _m3_multi.m   -> FIM_Modell*.mat

5. main_06_ovp_modell3.m                     -> OVP_Modell3.mat
```

Die Schritte 4 und 5 lesen die in Schritt 2 gespeicherten Parameter. Nach
einem neuen Fit müssen sie **beide** neu laufen — `main_06` prüft das und
bricht sonst mit einer Meldung ab.

---

## Ordnerstruktur

```
Daten/
  MessDaten_SPI1_Projekt/     Rohdaten (Mess_RamScDef*.mat)
  Daten_Processed/            aufbereitete Daten
  p_opt/                      identifizierte Parameter
  FIM/                        FIM- und OVP-Ergebnisse
Modelle/                      Modell- und Jacobi-Dateien
utils/                        plot_gaussian_ellipsoid, Parameteranalyse, ...
Bilder/                       erzeugte Abbildungen (PNG 300 dpi + PDF)
```

Die Skripte erwarten `Daten/` unterhalb des Ausführungsverzeichnisses und
`Modelle/`, `utils/` eine Ebene darüber. Jedes Skript definiert dafür ganz
oben eine Variable `DATEN` bzw. `projectRoot` — bei abweichender Struktur
reicht es, die anzupassen.

---

## Offene Punkte

* **Bootstrap** als Gegenprobe zur FIM wurde für Modell 1/2 nicht neu
  gerechnet (die alten Skripte nutzen `ValData` statt `TrainData` und noch
  freies `KLa`). Für Modell 3 ist Bootstrap nicht praktikabel: 300 Fits ×
  etwa 1 h.
* **`YPhX` und `YXO_eff`** haben relative Unsicherheiten von 27 % bzw. 23 %
  und lassen sich durch **keinen** Glucose-Feedverlauf verbessern. Phosphat
  limitiert in keinem zulässigen Versuch, DOT ist quasistationär. Beide
  bräuchten eine andere Stellgröße (Rührerdrehzahl) oder phosphatlimitierte
  Bedingungen.
* **Ethanol** ist in keinem Datensatz gemessen. Die Overflow-Kinetik lässt
  sich deshalb nicht identifizieren, nur belegen.
* Nur das **A-Kriterium** wurde in der OVP verwendet. Es minimiert die Summe
  der relativen Varianzen und verbessert daher bevorzugt die ohnehin gut
  bestimmten Richtungen — die Konditionszahl steigt leicht (562 → 613). Ein
  modifiziertes E-Kriterium würde stattdessen die schlechteste Richtung
  angehen.