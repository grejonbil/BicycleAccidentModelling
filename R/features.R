# =============================================================================
# features.R
# Infrastructure layers and feature computation (minimum distance intersection
# calculation, MDIC) for arbitrary points. Sourced by prepare_data.R and setup.R.
# =============================================================================

suppressPackageStartupMessages({ library(sf); library(dplyr) })
sf_use_s2(FALSE)
if (!exists("data_dir")) data_dir <- Sys.getenv("BAM_DATA_DIR", file.path("..", "Data"))

infra <- function(n) st_read(file.path(data_dir, "infrastructure", paste0(n, ".gpkg")), quiet = TRUE)

RoadNetwork         <- infra("road_network")
TrafficArea         <- infra("traffic_zones")
Exposure            <- infra("traffic_model_gvm2018")
BicycleNetwork      <- infra("bicycle_network")
one_way             <- infra("one_way_streets")

# layers to which the minimum Euclidean distance is computed (MDIC strategy 1)
cols_dist  <- c("RedLight", "PedestrianCrossing", "StopSign", "Junction", "ROW", "TramTrack",
                "ComplexIntersection", "BusLane", "PPS", "PT_Stops")
sf_objects <- lapply(c("traffic_signals", "pedestrian_crossings", "stop_signs", "junctions", "right_of_way_junctions",
                       "tram_tracks", "complex_intersections", "bus_lines", "public_parking", "public_transport_stops"), infra)

# Features of a set of points:
#   * distance [m] to the nearest feature of each layer in `sf_objects`
#   * zonentyp_t: speed regime of the traffic zone (20, 30; 50 outside zones; 0 pedestrian zone)
#   * SpedLmt, RodWdth, RodSrfc: speed limit, width class and surface of the nearest road segment
#   * Exposition: average weekday traffic of the nearest GVM-ZH 2018 segment
#   * bkTyp_b: nearest bicycle-network segment (0 = no dedicated facility, 1 = lane, 2 = path)
#   * direction: 1 = within 10 m of a one-way street, 2 = otherwise
compute_point_features <- function(points) {
  g <- st_geometry(points)
  for (k in seq_along(sf_objects)) {
    y <- st_geometry(sf_objects[[k]])
    points[[cols_dist[k]]] <- as.numeric(st_distance(g, y[st_nearest_feature(g, y)], by_element = TRUE))
  }
  idx <- st_intersects(g, TrafficArea)
  points$zonentyp_t <- vapply(idx, function(i) if (length(i)) as.numeric(TrafficArea$zonentyp_t[i[1]]) else 50, 0)
  ir <- st_nearest_feature(g, RoadNetwork)
  points$SpedLmt <- RoadNetwork$SpedLmt[ir]
  points$RodWdth <- RoadNetwork$RodWdth[ir]
  points$RodSrfc <- RoadNetwork$RodSrfc[ir]
  points$Exposition <- Exposure$Exposition[st_nearest_feature(g, Exposure)]
  points$bkTyp_b <- BicycleNetwork$bkTyp_b[st_nearest_feature(g, BicycleNetwork)]
  io <- st_nearest_feature(g, one_way)
  points$direction <- ifelse(as.numeric(st_distance(g, st_geometry(one_way)[io], by_element = TRUE)) <= 10, 1, 2)
  points
}
