# =============================================================================
# prepare_data.R
# Builds all model inputs from the original open data sources:
#   1. bicycle accidents 2016-2025 (City of Zurich)
#   2. infrastructure layers (City/Canton of Zurich OGD, swissTLM3D, OpenStreetMap)
#   3. derived layers: junctions, right-of-way junctions, complex intersections
#   4. road buffers, pseudo-absence points, modelling table and prediction grids
#
# Usage (from the R/ directory):  Rscript prepare_data.R
#
# Downloads are cached in <BAM_DATA_DIR>/raw. Three inputs are not downloaded:
#   - GVM-ZH 2018 traffic model (restricted; Amt fuer Mobilitaet, Canton of Zurich):
#       raw/GVM18.shp with attribute "Exposition"
#   - city boundary, lake, forest, buildings and green areas (Canton/City of Zurich):
#       raw/context/{city_boundary,lake,forest,buildings,green_areas}.gpkg
#   - complex intersections and pseudo-absence points of the original study (2023):
#       raw/complex_intersections_2023.gpkg, raw/pseudo_absence_points_2023.gpkg
#       (provided in this repository under data/)
# See docs/data.md for details.
# =============================================================================

suppressPackageStartupMessages({ library(sf); library(dplyr); library(igraph) })
sf_use_s2(FALSE)
options(timeout = 3600)

data_dir <- Sys.getenv("BAM_DATA_DIR", file.path("..", "Data"))
raw_dir  <- file.path(data_dir, "raw")
for (d in c(raw_dir, file.path(data_dir, c("accidents", "infrastructure", "context", "prediction_grid"))))
  dir.create(d, recursive = TRUE, showWarnings = FALSE)
YEARS <- 2016:2025
crs <- 2056
log_msg <- function(...) message(format(Sys.time(), "%H:%M:%S"), "  ", ...)
out <- function(x, ...) { f <- file.path(data_dir, ...); st_write(x, f, delete_dsn = TRUE, quiet = TRUE); log_msg(f, ": ", nrow(x)) }

# --- download helpers ----------------------------------------------------------
fetch <- function(url, file) {
  f <- file.path(raw_dir, file)
  if (!file.exists(f)) { log_msg("download ", file); download.file(url, f, mode = "wb", quiet = TRUE) }
  f
}
# City of Zurich geoportal (QGIS server WFS, GeoJSON in WGS84)
city_wfs <- function(service, typename, file) {
  st_transform(st_read(fetch(paste0("https://www.ogd.stadt-zuerich.ch/wfs/geoportal/", service,
    "?SERVICE=WFS&VERSION=1.1.0&REQUEST=GetFeature&TYPENAME=", typename, "&outputFormat=GeoJSON"), file), quiet = TRUE), crs)
}
# Canton of Zurich OGD WFS
canton_wfs <- function(layer, file, geom_type = NULL) {
  f <- file.path(raw_dir, file)
  if (!file.exists(f))
    sf::gdal_utils("vectortranslate", "WFS:https://maps.zh.ch/wfs/OGDZHWFS", f,
                   options = c(layer, "-spat", "2676000", "1241000", "2690000", "1255000", "-spat_srs", "EPSG:2056",
                               "-t_srs", "EPSG:2056", if (!is.null(geom_type)) c("-nlt", geom_type)))
  st_read(f, quiet = TRUE)
}

# --- context ---------------------------------------------------------------------
ctx <- function(n) st_transform(st_read(file.path(raw_dir, "context", paste0(n, ".gpkg")), quiet = TRUE), crs)
city <- st_union(ctx("city_boundary"))
for (n in c("city_boundary", "lake", "forest")) out(ctx(n), "context", paste0(n, ".gpkg"))
in_city <- function(x) x[lengths(st_intersects(x, city)) > 0, ]

# =============================================================================
# 1. Accidents
acc_csv <- fetch("https://data.stadt-zuerich.ch/dataset/sid_dav_strassenverkehrsunfallorte/download/RoadTrafficAccidentLocations.csv",
                 "RoadTrafficAccidentLocations.csv")
acc <- read.csv(acc_csv, encoding = "UTF-8") %>%
  filter(AccidentInvolvingBicycle == "true", AccidentYear %in% YEARS) %>%
  transmute(id = AccidentUID, year = AccidentYear, month = AccidentMonth, severity = AccidentSeverityCategory_en,
            type = AccidentType_en, E = AccidentLocation_CHLV95_E, N = AccidentLocation_CHLV95_N) %>%
  st_as_sf(coords = c("E", "N"), crs = crs)
out(acc, "accidents", "accidents.gpkg")

# =============================================================================
# 2. Infrastructure layers
speed <- city_wfs("Signalisierte_Geschwindigkeiten", "view_geoserver_tempo_ist", "signalisierte_geschwindigkeiten.geojson")
speed_map <- c(T0 = 0, T20 = 20, T30 = 30, T30N0 = 30, T50 = 50, T50N30 = 50, T60 = 60, T80 = 80, T100 = 100, T120 = 120)
speed$SpedLmt <- unname(speed_map[speed$temporegime_technical])

# road width and surface from swissTLM3D (release 2026-02)
tlm_file <- file.path(raw_dir, "tlm_strassen_zurich.gpkg")
if (!file.exists(tlm_file)) {
  z <- fetch("https://data.geo.admin.ch/ch.swisstopo.swisstlm3d/swisstlm3d_2026-02/swisstlm3d_2026-02_2056_5728.gpkg.zip", "swisstlm3d.gpkg.zip")
  sf::gdal_utils("vectortranslate", paste0("/vsizip/", z, "/SWISSTLM3D_2026_LV95_LN02.gpkg"), tlm_file,
                 options = c("tlm_strassen_strasse", "-spat", "2675000", "1239000", "2692000", "1257000", "-dim", "XY"))
}
width_map <- c("10m Strasse" = 10, "8m Strasse" = 8, "6m Strasse" = 6, "4m Strasse" = 4, "3m Strasse" = 3,
               "2m Weg" = 2, "1m Weg" = 1, "Autobahn" = 10, "Autostrasse" = 10, "Einfahrt" = 6, "Ausfahrt" = 6)
tlm <- st_read(tlm_file, quiet = TRUE) %>% filter(objektart %in% names(width_map)) %>% st_zm()
tlm$width <- unname(width_map[tlm$objektart]); tlm$surface <- as.integer(tlm$belagsart != "Hartbelag")

# road network: speed-limit segments with the modal TLM width/surface along each segment
roads <- speed %>% filter(!is.na(SpedLmt)) %>% st_cast("MULTILINESTRING") %>% st_cast("LINESTRING")
smp <- st_line_sample(st_geometry(roads), n = 5)
smp_pts <- st_cast(smp, "POINT"); seg_id <- rep(seq_len(nrow(roads)), each = 5)   # 5 points per segment
nn <- st_nearest_feature(smp_pts, tlm)
mode_of <- function(v) { u <- unique(v); u[which.max(tabulate(match(v, u)))] }
roads$RodWdth <- tapply(tlm$width[nn], seg_id, mode_of)[as.character(seq_len(nrow(roads)))]
roads$RodSrfc <- tapply(tlm$surface[nn], seg_id, mode_of)[as.character(seq_len(nrow(roads)))]
roads <- roads %>% transmute(name = lokalisationsname, SpedLmt, RodWdth = as.numeric(RodWdth), RodSrfc = as.numeric(RodSrfc))
out(roads, "infrastructure", "road_network.gpkg")

zones <- city_wfs("Verkehrszonen", "view_geoserver_zonen_ist", "verkehrszonen.geojson")
zones$zonentyp_t <- as.numeric(sub("T", "", zones$zonentyp_technical))
out(zones %>% select(zonenname, zonentyp, zonentyp_t), "infrastructure", "traffic_zones.gpkg")

one_way <- city_wfs("Erlaubte_Fahrtrichtungen", "view_geoserver_einbahn_ist", "erlaubte_fahrtrichtungen.geojson")
out(one_way %>% select(name = lokalisationsname, einbahnregime), "infrastructure", "one_way_streets.gpkg")

tram <- city_wfs("VBZ_Infrastruktur_OGD", "strecke_schienen", "vbz_schienen.geojson")
out(tram %>% select(objectid), "infrastructure", "tram_tracks.gpkg")

parking <- city_wfs("oeffentlich_zugaengliche_Parkplaetze_DAV", "oeff_strassenparkierung_dav_p", "parkplaetze.geojson") %>% st_cast("POINT")
out(parking %>% select(art, kategorie), "infrastructure", "public_parking.gpkg")

routes <- city_wfs("Fuss__und_Velowegnetz", "tbl_routennetz", "fuss_und_velowegnetz.geojson")
bicycle <- routes %>% filter(velo == 1) %>%
  mutate(bkTyp_b = case_when(veloweg == 1 ~ 2, velostreifen %in% c("1", "BOTH", "FT", "TF") ~ 1, TRUE ~ 0),
         bikeTyp = c("no dedicated facility", "bicycle lane", "bicycle path")[bkTyp_b + 1]) %>% select(bikeTyp, bkTyp_b)
out(bicycle, "infrastructure", "bicycle_network.gpkg")

pt_lines <- canton_wfs("ogd-0108_giszhpub_zvv_linien_l", "zvv_linien.gpkg", geom_type = "MULTILINESTRING")
bus <- in_city(pt_lines %>% filter(betriebszweig %in% c("02", "03", "04", "05", "06", "07", "08", "10", "31", "51", "61")))
out(bus %>% select(liniennummer, betriebszweig_txt), "infrastructure", "bus_lines.gpkg")
pt_stops <- in_city(canton_wfs("ogd-0140_giszhpub_zvv_haltestellen_p", "zvv_haltestellen.gpkg"))
out(pt_stops %>% select(chstname, vtyp), "infrastructure", "public_transport_stops.gpkg")

# OpenStreetMap: traffic signals, stop signs, marked pedestrian crossings
osm <- st_read(fetch("https://download.bbbike.org/osm/bbbike/Zuerich/Zuerich.osm.pbf", "Zuerich.osm.pbf"), layer = "points",
               query = "SELECT osm_id, highway, other_tags FROM points WHERE highway IN ('traffic_signals','stop','crossing')", quiet = TRUE) %>%
  st_transform(crs) %>% in_city()
out(osm %>% filter(highway == "traffic_signals"), "infrastructure", "traffic_signals.gpkg")
out(osm %>% filter(highway == "stop"), "infrastructure", "stop_signs.gpkg")
marked <- grepl('"crossing"=>"(zebra|marked|traffic_signals|uncontrolled)"|"crossing:markings"=>"(zebra|yes|lines|dashes|surface|dots)"', osm$other_tags)
crossings <- osm %>% filter(highway == "crossing", marked)
out(crossings, "infrastructure", "pedestrian_crossings.gpkg")

# GVM-ZH 2018 traffic model (restricted input, see header)
exposure <- st_read(file.path(raw_dir, "GVM18.shp"), quiet = TRUE) %>% st_transform(crs) %>% select(Exposition)
out(exposure, "infrastructure", "traffic_model_gvm2018.gpkg")

# =============================================================================
# 3. Derived layers
# Junctions: intersection points of road segments with different street names
# (points closer than 5 m are merged).
#   right-of-way junctions: two roads with <= 30 km/h meet ("30/30")
#   junctions (main/side road): a <= 30 km/h road meets a 50 km/h road ("30/50")
drive <- roads %>% filter(SpedLmt > 0)

# Step 1: the intersection point(s) of road i and j, classified by their speed limits.
intersection_type <- function(i, j) {
  lo <- min(drive$SpedLmt[c(i, j)]); hi <- max(drive$SpedLmt[c(i, j)])
  if (hi <= 30) "30/30" else if (lo <= 30 && hi == 50) "30/50" else "other"
}
road_pair_points <- function(i, j) {
  g <- st_intersection(st_geometry(drive)[i], st_geometry(drive)[j])
  if (!length(g)) return(NULL)
  st_sf(type = intersection_type(i, j), geometry = st_cast(st_cast(g, "MULTIPOINT"), "POINT"))
}
# all crossings of road i with later, differently-named roads it touches
road_crossings <- function(i) {
  js <- road_pairs[[i]]; js <- js[js > i & drive$name[js] != drive$name[i]]
  if (!length(js)) return(NULL)
  do.call(rbind, lapply(js, road_pair_points, i = i))
}
road_pairs <- st_intersects(drive, drive)
jx <- do.call(rbind, lapply(seq_len(nrow(drive)), road_crossings))
st_crs(jx) <- crs

# Step 2: merge intersection points within 5 m into one junction, keeping all types present.
jx$group <- components(graph_from_adj_list(st_is_within_distance(jx, jx, 5), mode = "all"))$membership
jx <- jx %>% group_by(group) %>% summarise(types = paste(sort(unique(type)), collapse = "+"), .groups = "drop") %>% st_centroid()

# Step 3: split into right-of-way junctions and (main/side-road) junctions, and write out.
row_junctions <- jx %>% filter(grepl("30/30", types))
junctions <- jx %>% filter(grepl("30/50", types))
out(row_junctions, "infrastructure", "right_of_way_junctions.gpkg")
out(junctions, "infrastructure", "junctions.gpkg")

# Complex intersections: locations where bus lines, pedestrian crossings, car roads,
# cycle routes and tram lines intersect within a 50 m radius. The 47 locations
# identified in the original study (raw/complex_intersections_2023.gpkg) are kept
# if they still meet all criteria with the current data.
complex_2023 <- st_read(file.path(raw_dir, "complex_intersections_2023.gpkg"), quiet = TRUE) %>% st_set_crs(crs)
near <- function(layer) lengths(st_is_within_distance(complex_2023, layer, 50)) > 0
complex <- complex_2023[near(bus) & near(crossings) & near(bicycle) & near(tram) & near(roads), ]
out(complex, "infrastructure", "complex_intersections.gpkg")

# =============================================================================
# 4. Road buffers, pseudo-absence points, modelling table, prediction grids
source("features.R")    # compute_point_features() and the layer list

buildings <- ctx("buildings"); green <- ctx("green_areas")
not_road <- st_union(st_buffer(st_union(rbind(st_geometry(buildings) %>% st_sf(), st_geometry(green) %>% st_sf())), 0))
road_buffer <- function(w) {
  b <- st_intersection(st_union(st_buffer(st_geometry(roads %>% filter(SpedLmt > 0, SpedLmt < 80)), w)), city)
  st_sf(geometry = st_difference(b, not_road))
}
for (w in c(25, 50)) if (!file.exists(file.path(data_dir, "context", sprintf("road_buffer_%dm.gpkg", w))))
  out(road_buffer(w), "context", sprintf("road_buffer_%dm.gpkg", w))

# Pseudo-absence points: the 10,000 PAPs of the original study (2023), placed along the
# cycling network (raw/pseudo_absence_points_2023.gpkg), without those now within 5 m
# of a recorded accident. Length-proportional random PAPs are evaluated as a
# sensitivity analysis in model_evaluation.R (stage "pseudoabs").
pap <- st_read(file.path(raw_dir, "pseudo_absence_points_2023.gpkg"), quiet = TRUE) %>% st_set_crs(crs)
st_geometry(pap) <- "geometry"
pap <- pap[lengths(st_is_within_distance(pap, acc, 5)) == 0, ]; pap$id <- paste0("PAP", seq_len(nrow(pap)))
out(pap, "accidents", "pseudo_absence_points.gpkg")

model_table <- rbind(acc %>% transmute(id, year, accidnt = 1), pap %>% transmute(id, year = NA_integer_, accidnt = 0))
model_table <- compute_point_features(model_table)
out(model_table, "accidents", "model_table.gpkg")

# Prediction grids: one random point in each of n randomly selected cells of a
# square grid (cell size s) covering the road buffer (width w), cf. Table 2
grid_points <- function(n, w, s) {
  buf <- st_read(file.path(data_dir, "context", sprintf("road_buffer_%dm.gpkg", w)), quiet = TRUE)
  pts <- st_sample(buf, size = round(n * 1.6)); xy <- st_coordinates(pts)
  cell <- paste(floor(xy[, 1] / s), floor(xy[, 2] / s)); pts <- pts[!duplicated(cell)]
  st_sf(geometry = pts[sample(length(pts), min(n, length(pts)))])
}
set.seed(2025)
for (cfg in list(c(200000, 25, 5), c(200000, 50, 50), c(15000, 50, 5))) {
  f <- sprintf("grid_%d_points_%dm_buffer_%dm_cells.gpkg", cfg[1], cfg[2], cfg[3])
  if (file.exists(file.path(data_dir, "prediction_grid", f))) next
  out(compute_point_features(grid_points(cfg[1], cfg[2], cfg[3])), "prediction_grid", f)
}
log_msg("data preparation finished")
