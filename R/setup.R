# =============================================================================
# setup.R
# Shared configuration for model training and evaluation:
#   - infrastructure layers and feature computation (features.R)
#   - modelling table, model formula and hyperparameters
#   - prediction grid and infrastructure scenarios (Table 1)
#
# Run prepare_data.R first. The data directory can be changed with the
# environment variable BAM_DATA_DIR (default: ../Data).
# =============================================================================

data_dir <- Sys.getenv("BAM_DATA_DIR", file.path("..", "Data"))
if (!dir.exists(file.path(data_dir, "accidents")))
  stop("Prepared data not found in ", normalizePath(data_dir, mustWork = FALSE), ". Run prepare_data.R or set BAM_DATA_DIR.")
source("features.R")
suppressPackageStartupMessages(library(caret))

ctx_layer <- function(n) st_read(file.path(data_dir, "context", paste0(n, ".gpkg")), quiet = TRUE)
accidents <- st_read(file.path(data_dir, "accidents", "accidents.gpkg"), quiet = TRUE)

# --- Modelling table: accidents 2016-2025 and pseudo-absence points --------------
# Pedestrian zones (0 km/h) and roads with >= 80 km/h are excluded.
model_sf <- st_read(file.path(data_dir, "accidents", "model_table.gpkg"), quiet = TRUE) %>%
  filter(SpedLmt < 80, SpedLmt != 0)
model_data <- model_sf %>% st_drop_geometry() %>% as.data.frame() %>%
  mutate(accidntFactorClass = factor(accidnt, levels = c(0, 1), labels = c("Class0", "Class1")),
         y = as.integer(accidnt == 1))
xy <- st_coordinates(model_sf); model_data$X <- xy[, 1]; model_data$Y <- xy[, 2]

RoadNetwork <- RoadNetwork %>% filter(SpedLmt < 80, SpedLmt != 0)

# --- Model specification -----------------------------------------------------------
predictors <- c("PedestrianCrossing", "StopSign", "TramTrack", "ComplexIntersection", "RodWdth",
                "Exposition", "zonentyp_t", "ROW", "PPS", "bkTyp_b", "BusLane", "SpedLmt",
                "PT_Stops", "direction")
interactions <- c("RodWdth:SpedLmt", "SpedLmt:Exposition")
make_formula <- function(terms) as.formula(paste("accidntFactorClass ~", paste(terms, collapse = " + ")))
model_formula <- make_formula(c(predictors, interactions))

paramGrid_gbm <- expand.grid(interaction.depth = 18, n.trees = 300, shrinkage = 0.05, n.minobsinnode = 20)
paramGrid_rf  <- expand.grid(mtry = 5)
paramGrid_xgb <- expand.grid(eta = 0.1, max_depth = 8, gamma = 5, colsample_bytree = 0.5,
                             min_child_weight = 10, nrounds = 250, subsample = 0.5)

feature_labels <- c(PedestrianCrossing = "PedestrianCrossing", StopSign = "StopSign", TramTrack = "TramTrack",
                    ComplexIntersection = "ComplexIntersection", RodWdth = "RoadWidth", Exposition = "Exposure",
                    zonentyp_t = "TrafficArea", ROW = "ROW", PPS = "PPS", bkTyp_b = "BicycleNetwork",
                    BusLane = "BusLane", SpedLmt = "SpeedLimit", PT_Stops = "PT_Stops", direction = "Direction",
                    RedLight = "RedLight", Junction = "Junction")

# --- Prediction grids (Table 2); the first is used for all scenario analyses -------
grid_configs <- data.frame(points = c(200000, 200000, 15000), buffer = c(25, 50, 50), cell = c(5, 50, 5))
grid_file <- function(points, buffer, cell)
  file.path(data_dir, "prediction_grid", sprintf("grid_%d_points_%dm_buffer_%dm_cells.gpkg", points, buffer, cell))
load_grid <- function(points = 200000, buffer = 25, cell = 5)
  st_read(grid_file(points, buffer, cell), quiet = TRUE) %>% subset(SpedLmt < 80 & SpedLmt != 0)

# --- Infrastructure scenarios (Table 1) --------------------------------------------
# bkTyp_b codes: 0 = no dedicated bicycle facility (Type 1), 1 = bicycle lane (Type 2),
#                2 = bicycle path (Type 3).
cap_speed_30 <- function(v) ifelse(v >= 30, 30, v)          # existing 20 km/h zones are retained
scenarios <- list(
  Baseline  = function(d) d,
  Optimal   = function(d) { d$SpedLmt <- cap_speed_30(d$SpedLmt); d$RodWdth <- d$RodWdth * 1.1
                            d$Exposition <- 2000; d$bkTyp_b <- 2; d },
  Congested = function(d) { d$SpedLmt <- 50; d$RodWdth <- d$RodWdth * 0.8
                            d$Exposition <- pmax(d$Exposition, 8000); d$bkTyp_b <- 0; d },
  Realistic = function(d) { d$SpedLmt <- cap_speed_30(d$SpedLmt); d },
  Control1  = function(d) { d$SpedLmt <- cap_speed_30(d$SpedLmt); d$RodWdth <- d$RodWdth * 1.3
                            d$Exposition <- 8000; d$bkTyp_b <- 0; d },
  Control2  = function(d) { d$SpedLmt <- 50; d$RodWdth <- d$RodWdth * 0.8
                            d$Exposition <- 2000; d$bkTyp_b <- 1; d }
)
scenario_labels <- c(Baseline = "Baseline", Optimal = "Optimal", Congested = "Congested",
                     Realistic = "Realistic", Control1 = "Control 1", Control2 = "Control 2")
