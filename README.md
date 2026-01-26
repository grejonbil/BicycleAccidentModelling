# GEO511 – Master Thesis Programming

Dieses Repository enthält den Code zur Modellierung und Vorhersage von Unfall-Lokationen sowie die zugehörige Datenverarbeitung für die Masterarbeit (GEO511).

## Inhalt

- `R/Location_Prediction_Clean.Rmd`  
  Zentrales Analyse‑Notebook: Datenaufbereitung, Feature‑Engineering, Modelltraining (GBM/RF/XGB), Generierung von Zufallspunkten, Szenario‑Vorhersagen und Raster‑Visualisierung.
- `Data/`, `GIS/`, `JMP/`, `Literature/`  
  Daten, GIS‑Material, Zusatzanalysen und Literatur (lokal, nicht versioniert oder nur teilweise enthalten).
- `random_grid.py`, `random_grid.ipynb`, `dissolve.ipynb`  
  Hilfsskripte/Notebooks für ergänzende Verarbeitung.

## Voraussetzungen

In R benötigt (mindestens):  
`sf`, `raster`, `rasterVis`, `ggplot2`, `dplyr`, `caret`, `gbm`, `randomForest`, `xgboost`, `future`, `future.apply`, `boot`, `gridExtra`

Hinweis: Die exakten Paketversionen sind nicht fixiert. Falls Reproduzierbarkeit nötig ist, bitte ein `renv`/`packrat` Setup ergänzen.

## Datenstruktur (erwartet)

Die R‑Markdown erwartet eine definierte Ordnerstruktur, insbesondere:

- Eingabedaten unter `Data/` und/oder `GIS/`
- Zielausgabe unter `Prediction_Grid/Random/` (wird bei Bedarf erzeugt)

Pfadlogik und Dateinamen sind im Notebook festgelegt. Wenn du lokal abweichende Strukturen nutzt, passe die Pfade in `R/Location_Prediction_Clean.Rmd` an.

## Workflow (Kurzüberblick)

1. Daten laden und bereinigen  
2. Feature‑Berechnung (z. B. Distanzen zu Infrastruktur‑Objekten)  
3. Modelltraining (GBM, RF, XGB)  
4. Generierung von Zufallspunkten in Rasterzellen  
5. Szenario‑Vorhersagen + Rasterisierung/Plots

## Raster‑ und Punkt‑Parameter

Im Notebook wird über Kombinationen aus Rasterzellgröße und Anzahl modellierter Punkte iteriert:

- Rasterzellgrößen (Meter): `4, 5, 6, 8, 10, 15, 20, 30, 40, 50, 60, 70, 80`
- Anzahl Punkte: `15'000, 30'000, 45'000, 100'000, 200'000`

Die generierten Punkt‑Dateien werden so benannt:
`point_network_random_<grid>m_<points>.geojson`

Beispiel: `point_network_random_8m_100000.geojson`

## Ausführen

Öffne `R/Location_Prediction_Clean.Rmd` in RStudio und führe die relevanten Chunks aus.  
Viele Schritte sind `eval=FALSE`, um Laufzeit zu sparen. Entferne das Flag bewusst bei Bedarf.

## Ergebnisse

Output‑Dateien und Karten werden im Notebook konfiguriert (z. B. `Prediction_Grid/Random/` sowie Plot‑Ausgabeordner).

## Hinweise

- Die Verarbeitung kann je nach Punktanzahl und Zellgröße sehr rechenintensiv sein.  
- Achte auf ausreichend RAM und Laufzeit, insbesondere bei 200k Punkten.  
- Für parallele Verarbeitung wird `future` verwendet.

---
