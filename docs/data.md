# Data

All spatial data use the Swiss coordinate system CH1903+/LV95 (EPSG:2056). The data are not redistributed in this repository. `R/prepare_data.R` downloads the open data listed below into `<data directory>/raw/` and builds all model inputs.

## Sources

| Dataset | Source | Version used | Download |
|----|----|----|----|
| Police-reported road traffic accidents with bicycle involvement, 2016–2025 (Accident Assessment Protocol of the Federal Roads Office, FEDRO) | City of Zurich, Open Data Zurich (`sid_dav_strassenverkehrsunfallorte`) | 2026 | automatic |
| Signposted speed limits (road network) | City of Zurich geoportal WFS `Signalisierte_Geschwindigkeiten` | 2026 | automatic |
| Road width and surface | swisstopo, swissTLM3D (`tlm_strassen_strasse`) | 2026-02 | automatic |
| Traffic zones (20/30 km/h zones, pedestrian zones) | City of Zurich geoportal WFS `Verkehrszonen` | 2026 | automatic |
| One-way streets | City of Zurich geoportal WFS `Erlaubte_Fahrtrichtungen` | 2026 | automatic |
| Tram tracks | City of Zurich geoportal WFS `VBZ_Infrastruktur_OGD` | 2026 | automatic |
| Public on-street parking | City of Zurich geoportal WFS `oeffentlich_zugaengliche_Parkplaetze_DAV` | 2026 | automatic |
| Bicycle network (bicycle lanes and paths) | City of Zurich, `Fuss- und Velowegnetz` (`tbl_routennetz`) | 2022 (latest) | automatic |
| Bus lines, public transport stops (ZVV) | Canton of Zurich OGD WFS (`ogd-0108`, `ogd-0140`) | 2026 | automatic |
| Traffic signals, stop signs, marked pedestrian crossings | OpenStreetMap (BBBike extract Zürich) | 2026 | automatic |
| Traffic model GVM-ZH 2018 (average weekday traffic per road segment) | Office of Transport of the Canton of Zurich (Amt für Mobilität); restricted, available on request | 2018 (latest) | manual: `raw/GVM18.shp` (attribute `Exposition`) |
| City boundary, lake, forest, buildings, green areas | Canton and City of Zurich open geodata | 2026 | manual: `raw/context/{city_boundary,lake,forest,buildings,green_areas}.gpkg` |
| Complex intersections of the original study (2023) | derived, provided in this repository | 2023 | copy `data/complex_intersections_2023.gpkg` to `raw/` |
| Pseudo-absence points of the original study (2023) | derived, provided in this repository | 2023 | copy `data/pseudo_absence_points_2023.gpkg` to `raw/` |

## Derived layers

- **Junctions** are intersections of differently named road segments (merged within 5 m). Right-of-way junctions are junctions at which only roads with speed limits of at most 30 km/h meet; junctions between main and side roads connect 30 km/h and 50 km/h roads.
- **Complex intersections** are locations where bus lines, pedestrian crossings, car roads, bicycle routes and tram tracks meet within 50 m. The 47 locations identified in the original study are kept where all five layers are still present within 50 m (45 locations).
- **Pseudo-absence points** are the 10,000 points of the original study along the cycling network (at least 5 m from accidents and 10 m from each other), without those within 5 m of an accident recorded in 2016–2025.
- **Road buffers** of 25 m and 50 m around the road network, without buildings and green areas, define the area of the prediction grids.

## Variables

| Variable (code) | Name in the manuscript | Description | Feature type |
|----|----|----|----|
| `PedestrianCrossing` | PedestrianCrossing | Distance to the nearest marked pedestrian crossing [m] | distance |
| `StopSign` | StopSign | Distance to the nearest stop sign [m] | distance |
| `TramTrack` | TramTrack | Distance to the nearest tram track [m] | distance |
| `ComplexIntersection` | ComplexIntersection | Distance to the nearest complex intersection [m] | distance |
| `ROW` | ROW | Distance to the nearest right-of-way junction [m] | distance |
| `PPS` | PPS | Distance to the nearest public parking space [m] | distance |
| `BusLane` | BusLane | Distance to the nearest bus line [m] | distance |
| `PT_Stops` | PT_Stops | Distance to the nearest public transport stop [m] | distance |
| `RedLight` | RedLight | Distance to the nearest traffic signal [m] (candidate variable, not in the final model) | distance |
| `Junction` | Junction | Distance to the nearest junction of a main and a side road [m] (candidate variable, not in the final model) | distance |
| `SpedLmt` | SpeedLimit | Signposted speed limit of the nearest road segment (20, 30, 50, 60 km/h) | attribute |
| `RodWdth` | RoadWidth | Width class of the nearest road segment (1, 2, 3, 4, 6, 8, 10 m) | attribute |
| `Exposition` | Exposure | Average weekday traffic of the nearest segment (GVM-ZH 2018) [vehicles/day] | attribute |
| `zonentyp_t` | TrafficArea | Speed regime of the traffic area (20, 30; 50 outside 20/30 zones) | attribute |
| `bkTyp_b` | BicycleNetwork | Type of the nearest bicycle-network segment: 0 = no dedicated facility, 1 = bicycle lane, 2 = bicycle path | attribute |
| `direction` | Direction | 1 = within 10 m of a one-way street, 2 = otherwise | attribute |
| `accidnt` | – | 1 = accident, 0 = pseudo-absence point | response |

Points on roads with a speed limit of 0 km/h (pedestrian zones) or of 80 km/h and more are excluded from modelling.

## Directory structure

Paths are relative to the data directory (default `Data/` next to `R/`, or `BAM_DATA_DIR`). Everything except `raw/GVM18.shp`, `raw/context/`, `raw/complex_intersections_2023.gpkg` and `raw/pseudo_absence_points_2023.gpkg` is created by `prepare_data.R`.

```         
Data/
  raw/                                  downloaded and manually provided source data
  accidents/
    accidents.gpkg                      bicycle accidents 2016-2025
    pseudo_absence_points.gpkg          pseudo-absence points
    model_table.gpkg                    accidents + PAPs with all features
  infrastructure/
    road_network.gpkg                   SpedLmt, RodWdth, RodSrfc
    traffic_zones.gpkg  one_way_streets.gpkg  tram_tracks.gpkg  public_parking.gpkg
    bicycle_network.gpkg  bus_lines.gpkg  public_transport_stops.gpkg
    traffic_signals.gpkg  stop_signs.gpkg  pedestrian_crossings.gpkg
    traffic_model_gvm2018.gpkg  junctions.gpkg  right_of_way_junctions.gpkg  complex_intersections.gpkg
  context/
    city_boundary.gpkg  lake.gpkg  forest.gpkg  road_buffer_25m.gpkg  road_buffer_50m.gpkg
  prediction_grid/
    grid_<points>_points_<buffer>m_buffer_<cell>m_cells.gpkg
```

The prediction grid used for the scenario analyses contains 200,000 points within the 25 m road buffer (`grid_200000_points_25m_buffer_5m_cells.gpkg`).
