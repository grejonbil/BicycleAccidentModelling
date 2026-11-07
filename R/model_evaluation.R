# =============================================================================
# model_evaluation.R
# Validation, sensitivity, interpretation and scenario analyses of the bicycle
# accident location models (Sections 3.2-4.3 of the manuscript).
#
# Usage (from the R/ directory):
#   Rscript model_evaluation.R                 # run all stages
#   Rscript model_evaluation.R cv metrics      # run selected stages
#
# Stages (in dependency order):
#   cv          out-of-fold predictions: random and spatial block cross-validation
#   metrics     Table 4, calibration curves (Figure S1), classification thresholds
#   repeats     repeat-accident locations and de-duplicated model
#   features    Table S3 (feature-set sensitivity)
#   scenarios   scenario predictions (Figures 4-7), Table S2, Table 6
#   interpret   SHAP, permutation importance, partial dependence (Figures 8 and 9), Table S5
#   temporal    temporal-mismatch check (accidents 2022-2025 and 2016-2019 only)
#   pseudoabs   Table S1 (pseudo-absence sensitivity)
#   bootstrap   Tables 5 and S6 (bootstrap uncertainty of scenario differences)
#   maps        summaries of the scenario maps (Section 4.2)
#
# Outputs: ../results/tables, ../results/figure_data, ../results/key_results.csv
# The figures are drawn by figures.Rmd from these outputs.
# Intermediate objects are cached in ../results/cache (not versioned).
# Environment variables: BAM_DATA_DIR, BAM_GRID_FILE, BAM_RESULTS_DIR, BAM_CORES, BAM_BOOT
# =============================================================================

source("setup.R")
suppressPackageStartupMessages({
  library(tidyr); library(gbm); library(randomForest); library(xgboost); library(pROC)
  library(pdp); library(kernelshap); library(igraph); library(parallel)
})

stages_all <- c("cv", "metrics", "repeats", "features", "scenarios", "interpret",
                "temporal", "pseudoabs", "bootstrap", "maps")
args   <- commandArgs(trailingOnly = TRUE)
stages <- if (length(args)) args else stages_all
run    <- function(s) s %in% stages

NCORE  <- as.integer(Sys.getenv("BAM_CORES", max(1, detectCores() - 1)))
B_BOOT <- as.integer(Sys.getenv("BAM_BOOT", "100"))
res_dir <- Sys.getenv("BAM_RESULTS_DIR", file.path("..", "results"))
tab_dir <- file.path(res_dir, "tables"); fig_data_dir <- file.path(res_dir, "figure_data"); cache_dir <- file.path(res_dir, "cache")
for (d in c(tab_dir, fig_data_dir, cache_dir)) dir.create(d, recursive = TRUE, showWarnings = FALSE)
cache <- function(n) file.path(cache_dir, paste0(n, ".rds"))
save_cache <- function(x, n) { tmp <- paste0(cache(n), ".tmp"); saveRDS(x, tmp); file.rename(tmp, cache(n)) }   # atomic write
log_msg <- function(...) message(format(Sys.time(), "%H:%M:%S"), "  ", ...)

# key results (single numbers quoted in the text) are collected in one file
key_file <- file.path(res_dir, "key_results.csv")
key <- function(name, value, description) {
  old <- if (file.exists(key_file)) read.csv(key_file, stringsAsFactors = FALSE, colClasses = "character")
         else data.frame(name = character(), value = character(), description = character())
  write.csv(rbind(old[old$name != name, ], data.frame(name = name, value = as.character(value), description = description)),
            key_file, row.names = FALSE)
}
fmt <- function(x, d = 2) formatC(x, format = "f", digits = d)

# --- model helpers ---------------------------------------------------------------
ctrl_none <- trainControl(method = "none", classProbs = TRUE)
grid_of <- function(m) switch(m, gbm = paramGrid_gbm, rf = paramGrid_rf, xgbTree = paramGrid_xgb)
fit_model <- function(dat, method = "gbm", formula = model_formula) {
  a <- list(formula, data = dat, method = method, trControl = ctrl_none, tuneGrid = grid_of(method))
  if (method == "rf") a$ntree <- 500
  if (method == "gbm") a$verbose <- FALSE
  do.call(train, a)
}
prob <- function(fit, newdata) predict(fit, newdata = newdata, type = "prob")[, "Class1"]
brier   <- function(p, y) mean((p - y)^2)
logloss <- function(p, y, eps = 1e-15) { p <- pmin(pmax(p, eps), 1 - eps); -mean(y * log(p) + (1 - y) * log(1 - p)) }
auc_of  <- function(p, y) as.numeric(auc(roc(y, p, quiet = TRUE, direction = "<")))
scores  <- function(p, y) c(AUC = auc_of(p, y), Brier = brier(p, y), LogLoss = logloss(p, y))
model_names <- c(gbm = "GB", rf = "RF", xgbTree = "XGB")

# Spatial blocks: square cells of `size` m, randomly assigned to k folds
block_folds <- function(df, size, k = 10) {
  b <- paste(floor(df$X / size), floor(df$Y / size)); u <- unique(b)
  sample(rep_len(1:k, length(u)))[match(b, u)]
}
cv_predict <- function(df, folds, method = "gbm", formula = model_formula) {
  res <- mclapply(sort(unique(folds)), function(f) {
    fit <- fit_model(df[folds != f, ], method, formula); list(f = f, p = prob(fit, df[folds == f, ])) }, mc.cores = NCORE)
  out <- rep(NA_real_, nrow(df)); for (r in res) out[folds == r$f] <- r$p; out
}

set.seed(2026)
folds <- list(random = sample(rep_len(1:10, nrow(model_data))),
              b500 = block_folds(model_data, 500), b1000 = block_folds(model_data, 1000),
              b2000 = block_folds(model_data, 2000))
key("n_presences", sum(model_data$y), "accident points in the modelling dataset")
key("n_pseudo_absences", sum(model_data$y == 0), "pseudo-absence points in the modelling dataset")
key("ratio_presence_absence", fmt(sum(model_data$y == 0) / sum(model_data$y), 1), "effective presence:absence ratio 1:x")

# =============================================================================
# cv: out-of-fold predictions
if (run("cv")) {
  oof <- if (file.exists(cache("oof"))) readRDS(cache("oof")) else list()
  jobs <- c(paste(rep(names(model_names), each = 2), c("random", "b1000")), "gbm b500", "gbm b2000")
  for (j in setdiff(jobs, names(oof))) {
    m <- strsplit(j, " ")[[1]]; oof[[j]] <- cv_predict(model_data, folds[[m[2]]], m[1])
    save_cache(oof, "oof"); log_msg("cv ", j)
  }
}
oof <- if (file.exists(cache("oof"))) readRDS(cache("oof")) else NULL
need_oof <- function() if (is.null(oof)) stop("Run stage 'cv' first.")

# =============================================================================
# metrics: Table 4, calibration curves (Figure S1), thresholds
if (run("metrics")) {
  need_oof(); y <- model_data$y
  t4 <- do.call(rbind, lapply(names(model_names), function(m) {
    r <- scores(oof[[paste(m, "random")]], y); s <- scores(oof[[paste(m, "b1000")]], y)
    ci <- ci.auc(roc(y, oof[[paste(m, "b1000")]], quiet = TRUE, direction = "<"))
    data.frame(model = model_names[[m]], auc_random = r[["AUC"]], brier_random = r[["Brier"]], logloss_random = r[["LogLoss"]],
               auc_spatial = s[["AUC"]], auc_spatial_ci_low = ci[1], auc_spatial_ci_high = ci[3],
               brier_spatial = s[["Brier"]], logloss_spatial = s[["LogLoss"]]) }))
  write.csv(t4, file.path(tab_dir, "table_4_model_performance.csv"), row.names = FALSE)
  for (b in c("b500", "b1000", "b2000")) key(paste0("auc_gb_spatial_", sub("b", "", b), "m"), fmt(scores(oof[[paste("gbm", b)]], y)[["AUC"]]),
                                             paste("GB ROC-AUC, spatial block CV,", sub("b", "", b), "m blocks"))
  key("auc_gb_random", fmt(t4$auc_random[1]), "GB ROC-AUC, random 10-fold CV")
  rocs <- lapply(names(model_names), function(m) roc(y, oof[[paste(m, "b1000")]], quiet = TRUE, direction = "<"))
  key("delong_p_gb_rf", signif(roc.test(rocs[[1]], rocs[[2]], method = "delong")$p.value, 2), "DeLong test GB vs RF (spatial CV)")
  key("delong_p_gb_xgb", signif(roc.test(rocs[[1]], rocs[[3]], method = "delong")$p.value, 2), "DeLong test GB vs XGB (spatial CV)")

  cal <- do.call(rbind, lapply(names(model_names), function(m)
    data.frame(model = model_names[[m]], p = oof[[paste(m, "b1000")]], y = y) %>% mutate(decile = ntile(p, 10)) %>%
      group_by(model, decile) %>% summarise(mean_predicted = mean(p), observed_share = mean(y), .groups = "drop")))
  write.csv(cal, file.path(tab_dir, "calibration_curves.csv"), row.names = FALSE)

  r_gb <- rocs[[1]]; p_gb <- oof[["gbm b1000"]]
  youden <- coords(r_gb, "best", best.method = "youden", ret = "threshold")$threshold[1]
  ths <- seq(0.05, 0.95, 0.01)
  f1 <- sapply(ths, function(t) { pr <- p_gb >= t; tp <- sum(pr & y == 1); 2 * tp / (2 * tp + sum(pr & y == 0) + sum(!pr & y == 1)) })
  at80 <- coords(r_gb, 0.8, input = "threshold", ret = c("specificity", "sensitivity"))
  save_cache(list(youden = youden), "thresholds")
  key("threshold_youden", fmt(youden), "Youden-optimal threshold (GB, spatial CV)")
  key("threshold_f1", fmt(ths[which.max(f1)]), "F1-maximising threshold (GB, spatial CV)")
  key("specificity_at_0.8", fmt(at80$specificity), "specificity at threshold 0.8")
  key("sensitivity_at_0.8", fmt(at80$sensitivity), "sensitivity at threshold 0.8")
}

# =============================================================================
# repeats: locations with several accidents (within 5 m)
if (run("repeats")) {
  need_oof()
  acc <- model_sf[model_data$y == 1, ]
  cl <- components(graph_from_adj_list(st_is_within_distance(acc, acc, dist = 5), mode = "all"))$membership
  n_at <- table(cl); repeat_ids <- as.integer(names(n_at)[n_at >= 2])
  key("repeat_accidents", sum(n_at[n_at >= 2]), "accidents at locations with >= 2 accidents within 5 m")
  key("repeat_accidents_pct", fmt(100 * sum(n_at[n_at >= 2]) / sum(n_at), 1), "share of accidents at repeat locations [%]")
  key("repeat_locations", length(repeat_ids), "locations with >= 2 accidents within 5 m")
  key("repeat_max", max(n_at), "maximum number of accidents at one location")
  pres <- which(model_data$y == 1)
  dedup <- model_data[c(which(model_data$y == 0), pres[!duplicated(cl)]), ]
  set.seed(7)
  key("auc_gb_spatial_deduplicated", fmt(scores(cv_predict(dedup, block_folds(dedup, 1000)), dedup$y)[["AUC"]]),
      "GB ROC-AUC (1000 m blocks) with one presence point per location")
  p_pres <- oof[["gbm b1000"]][pres]; is_rep <- cl %in% repeat_ids
  key("prob_repeat_locations", fmt(mean(p_pres[is_rep])), "mean out-of-fold probability at repeat-accident locations")
  key("prob_single_locations", fmt(mean(p_pres[!is_rep])), "mean out-of-fold probability at single-accident locations")
  key("wilcoxon_p_repeat", signif(wilcox.test(p_pres[is_rep], p_pres[!is_rep])$p.value, 2), "Wilcoxon test repeat vs single")
}

# =============================================================================
# features: Table S3
if (run("features")) {
  need_oof(); base_auc <- scores(oof[["gbm b1000"]], model_data$y)[["AUC"]]
  variants <- c(setNames(lapply(predictors, function(v) c(setdiff(predictors, v), interactions[!grepl(v, interactions)])),
                         paste("without", feature_labels[predictors])),
                list("14 features, no explicit interactions" = predictors,
                     "16 features (+RedLight, +Junction)" = c(predictors, "RedLight", "Junction", interactions)))
  tA3 <- do.call(rbind, lapply(names(variants), function(n) {
    a <- scores(cv_predict(model_data, folds$b1000, "gbm", make_formula(variants[[n]])), model_data$y)[["AUC"]]; log_msg(n)
    data.frame(variant = n, n_features = sum(!grepl(":", variants[[n]])), auc_spatial = a, delta_auc = a - base_auc) }))
  write.csv(tA3, file.path(tab_dir, "table_S3_feature_sensitivity.csv"), row.names = FALSE)
}

# =============================================================================
# prediction grid and reference model (needed by the remaining stages)
if (any(sapply(c("scenarios", "interpret", "temporal", "pseudoabs", "bootstrap", "maps"), run))) {
  grid_sf <- load_grid(); grid_df <- st_drop_geometry(grid_sf)
  scen_X <- lapply(scenarios, function(f) f(grid_df))
  principal_axes <- grid_df$Exposition >= quantile(grid_df$Exposition, 0.9)
  if (!file.exists(cache("fit_reference"))) { set.seed(11); save_cache(fit_model(model_data), "fit_reference") }
  fit_ref <- readRDS(cache("fit_reference"))
}

# =============================================================================
# scenarios: scenario predictions (Figures 4-7), Table S2, Table 6
if (run("scenarios")) {
  youden <- readRDS(cache("thresholds"))$youden
  sp <- bind_rows(lapply(names(scen_X), function(s) data.frame(scenario = scenario_labels[[s]], p = prob(fit_ref, scen_X[[s]]))))
  sp$scenario <- factor(sp$scenario, levels = scenario_labels)
  save_cache(sp, "scenario_predictions")
  saveRDS(list(predictions = sp, youden = youden), file.path(fig_data_dir, "scenario_predictions.rds"))   # Figures 4-7
  t5 <- sp %>% group_by(scenario) %>% summarise(mean = mean(p), median = median(p), p90 = quantile(p, .9),
          pct_ge_0.5 = 100 * mean(p >= .5), pct_ge_youden = 100 * mean(p >= youden),
          pct_ge_0.8 = 100 * mean(p >= .8), pct_ge_0.9 = 100 * mean(p >= .9), .groups = "drop")
  write.csv(t5, file.path(tab_dir, "table_S2_scenario_thresholds.csv"), row.names = FALSE)
  t7 <- sp %>% group_by(scenario) %>% summarise(p05 = quantile(p, .05), p95 = quantile(p, .95), .groups = "drop")
  write.csv(t7, file.path(tab_dir, "table_6_scenario_probability_ranges.csv"), row.names = FALSE)
}

# =============================================================================
# interpret: SHAP and partial dependence (Figures 8 and 9), Table S5, importance measures
if (run("interpret")) {
  need_oof()
  X <- model_data[, predictors]
  set.seed(5); sub <- sample(nrow(X), 600); bg <- X[sample(nrow(X), 120), ]
  ks <- kernelshap(fit_ref, X = X[sub, ], bg_X = bg, pred_fun = function(object, newdata) prob(object, newdata), verbose = FALSE)
  shap <- ks$S; colnames(shap) <- feature_labels[colnames(shap)]
  Xs <- X[sub, ]; colnames(Xs) <- feature_labels[colnames(Xs)]
  imp <- sort(colMeans(abs(shap)), decreasing = TRUE)
  write.csv(data.frame(feature = names(imp), mean_abs_shap = imp,
                       spearman_value_vs_shap = sapply(names(imp), function(v) cor(Xs[[v]], shap[, v], method = "spearman"))),
            file.path(tab_dir, "shap_importance.csv"), row.names = FALSE)
  key("rank_speed_limit_shap", which(names(imp) == "SpeedLimit"), "rank of speed limit by mean |SHAP| (of 14)")
  key("rank_exposure_shap", which(names(imp) == "Exposure"), "rank of exposure by mean |SHAP| (of 14)")

  perm <- do.call(rbind, mclapply(1:5, function(f) {
    tr <- model_data[folds$b1000 != f, ]; te <- model_data[folds$b1000 == f, ]
    fit <- fit_model(tr); a0 <- auc_of(prob(fit, te), te$y); set.seed(f)
    sapply(predictors, function(v) { t2 <- te; t2[[v]] <- sample(t2[[v]]); a0 - auc_of(prob(fit, t2), t2$y) }) }, mc.cores = NCORE))
  write.csv(data.frame(feature = feature_labels[predictors], delta_auc = colMeans(perm)),
            file.path(tab_dir, "permutation_importance.csv"), row.names = FALSE)

  train_df <- model_data[, c(predictors, "accidntFactorClass")]
  pd1 <- partial(fit_ref, pred.var = c("SpedLmt", "Exposition"), prob = TRUE, which.class = "Class1", train = train_df,
                 pred.grid = expand.grid(SpedLmt = c(20, 30, 50, 60), Exposition = seq(0, 30000, by = 1500)))
  pd2 <- partial(fit_ref, pred.var = c("RodWdth", "SpedLmt"), prob = TRUE, which.class = "Class1", train = train_df,
                 pred.grid = expand.grid(RodWdth = c(1, 2, 3, 4, 6, 8, 10), SpedLmt = c(20, 30, 50, 60)))
  write.csv(pd1, file.path(tab_dir, "partial_dependence_speed_exposure.csv"), row.names = FALSE)
  write.csv(pd2, file.path(tab_dir, "partial_dependence_roadwidth_speed.csv"), row.names = FALSE)
  saveRDS(list(shap = shap, X = Xs), file.path(fig_data_dir, "shap_values.rds"))   # Figure 8

  # Table S5: co-occurring adverse characteristics
  t8 <- model_data %>%
    mutate(p = oof[["gbm b1000"]],
           n_adverse = pmin((ComplexIntersection <= 50) + (TramTrack <= 10) + (SpedLmt >= 50) +
                            (Exposition >= 8000) + (bkTyp_b == 0), 4)) %>%
    group_by(n_adverse) %>% summarise(n_points = n(), share_accident_points = mean(y), mean_predicted = mean(p),
                                      q25 = quantile(p, .25), q75 = quantile(p, .75), .groups = "drop")
  write.csv(t8, file.path(tab_dir, "table_S5_infrastructure_profiles.csv"), row.names = FALSE)
}

# =============================================================================
# temporal: accidents 2022-2025 only (closest to the infrastructure reference date)
# vs. a size-matched random subsample of all accidents 2016-2025
if (run("temporal")) {
  recent <- model_data[model_data$y == 0 | model_data$year >= 2022, ]; n_rec <- sum(recent$y)
  set.seed(3); matched <- model_data[c(which(model_data$y == 0), sample(which(model_data$y == 1), n_rec)), ]
  key("auc_gb_2022_2025", fmt(scores(cv_predict(recent, block_folds(recent, 1000)), recent$y)[["AUC"]]),
      sprintf("GB ROC-AUC (1000 m blocks), accidents 2022-2025 only (n = %d)", n_rec))
  key("auc_gb_size_matched", fmt(scores(cv_predict(matched, block_folds(matched, 1000)), matched$y)[["AUC"]]),
      "GB ROC-AUC (1000 m blocks), size-matched random subsample 2016-2025")
  early <- model_data[model_data$y == 0 | model_data$year <= 2019, ]
  key("auc_gb_2016_2019", fmt(scores(cv_predict(early, block_folds(early, 1000)), early$y)[["AUC"]]),
      sprintf("GB ROC-AUC (1000 m blocks), accidents 2016-2019 only (n = %d)", sum(early$y)))
  ref_surface <- prob(fit_ref, grid_df)
  key("spearman_2022_2025_surface", fmt(cor(prob(fit_model(recent), grid_df), ref_surface, method = "spearman")),
      "Spearman correlation of prediction surfaces (2022-2025 vs. full model)")
  key("spearman_2016_2019_surface", fmt(cor(prob(fit_model(early), grid_df), ref_surface, method = "spearman")),
      "Spearman correlation of prediction surfaces (2016-2019 vs. full model)")
}

# =============================================================================
# pseudoabs: Table S1 - regenerate PAPs on the bicycle network
if (run("pseudoabs")) {
  acc <- model_sf[model_data$y == 1, ]
  if (!file.exists(cache("pap_pool"))) {
    set.seed(42)
    rn <- st_cast(st_cast(st_geometry(BicycleNetwork), "MULTILINESTRING"), "LINESTRING")
    counts <- tabulate(sample(seq_along(rn), 60000, replace = TRUE, prob = as.numeric(st_length(rn))), nbins = length(rn))
    pts <- st_line_sample(rn[counts > 0], n = counts[counts > 0], type = "random")
    pool <- compute_point_features(st_sf(geometry = st_cast(pts[!st_is_empty(pts)], "POINT")))
    pool$d_acc <- as.numeric(st_distance(st_geometry(pool), st_geometry(acc)[st_nearest_feature(pool, acc)], by_element = TRUE))
    save_cache(pool, "pap_pool"); log_msg("PAP candidate pool")
  }
  pool <- readRDS(cache("pap_pool")); pool <- pool[pool$SpedLmt < 80 & pool$SpedLmt != 0, ]
  thin <- function(p, d) {   # greedy minimum-distance thinning (grid hash)
    xy <- st_coordinates(p); cx <- floor(xy[, 1] / d); cy <- floor(xy[, 2] / d)
    h <- new.env(hash = TRUE); keep <- logical(nrow(xy))
    for (i in seq_len(nrow(xy))) {
      ok <- TRUE
      for (dx in -1:1) { for (dy in -1:1) { js <- h[[paste(cx[i] + dx, cy[i] + dy)]]
        if (!is.null(js) && any((xy[js, 1] - xy[i, 1])^2 + (xy[js, 2] - xy[i, 2])^2 < d^2)) { ok <- FALSE; break } }
        if (!ok) break }
      if (ok) { keep[i] <- TRUE; k <- paste(cx[i], cy[i]); h[[k]] <- c(h[[k]], i) }
    }
    p[keep, ] }
  pres <- model_data[model_data$y == 1, c(predictors, "X", "Y", "y", "accidntFactorClass")]
  ref_surface <- prob(fit_ref, grid_df)
  scenario_summary <- function(f) {   # mean probabilities and differences to Baseline on the principal axes, all scenarios
    pb <- lapply(scen_X, function(d) prob(f, d))
    out <- data.frame(mean_baseline = mean(pb$Baseline), mean_optimal = mean(pb$Optimal), mean_congested = mean(pb$Congested),
                      axes_optimal_minus_baseline = mean((pb$Optimal - pb$Baseline)[principal_axes]))
    for (s in c("Realistic", "Control1", "Control2")) out[[paste0("mean_", tolower(s))]] <- mean(pb[[s]])
    for (s in c("Congested", "Realistic", "Control1", "Control2"))
      out[[paste0("axes_", tolower(s), "_minus_baseline")]] <- mean((pb[[s]] - pb$Baseline)[principal_axes])
    out }
  cfg <- rbind(data.frame(ratio = c(1, 2, 3, 5, 10), d_acc = 5, d_pap = 10),
               data.frame(ratio = 3, d_acc = c(10, 25), d_pap = 10), data.frame(ratio = 3, d_acc = 5, d_pap = c(25, 50)))
  tA4 <- do.call(rbind, lapply(seq_len(nrow(cfg)), function(i) {
    cc <- cfg[i, ]; n_target <- round(cc$ratio * nrow(pres)); set.seed(500 + i)
    cand <- pool[pool$d_acc >= cc$d_acc, ]; cand <- thin(cand[sample(nrow(cand)), ], cc$d_pap)
    cand <- cand[seq_len(min(n_target, nrow(cand))), ]
    pd <- st_drop_geometry(cand); xy <- st_coordinates(cand); pd$X <- xy[, 1]; pd$Y <- xy[, 2]
    pd$y <- 0L; pd$accidntFactorClass <- factor("Class0", levels = c("Class0", "Class1"))
    dat <- bind_rows(pres, pd[, names(pres)])
    s <- scores(cv_predict(dat, block_folds(dat, 1000)), dat$y); fit_i <- fit_model(dat); surf <- prob(fit_i, grid_df); log_msg("PAP config ", i)
    data.frame(cc, n_pap = nrow(cand), auc_spatial = s[["AUC"]], brier = s[["Brier"]],
               spearman_vs_reference = cor(surf, ref_surface, method = "spearman"), pct_ge_0.8 = 100 * mean(surf >= 0.8),
               scenario_summary(fit_i)) }))
  ref_row <- data.frame(ratio = NA, d_acc = 5, d_pap = NA, n_pap = sum(model_data$y == 0), auc_spatial = NA, brier = NA,
                        spearman_vs_reference = 1, pct_ge_0.8 = 100 * mean(ref_surface >= 0.8), scenario_summary(fit_ref))
  tA4 <- rbind(cbind(pap = "original (reference)", ref_row), cbind(pap = "length-proportional", tA4))
  write.csv(tA4, file.path(tab_dir, "table_S1_pseudo_absence_sensitivity.csv"), row.names = FALSE)
}

# =============================================================================
# bootstrap: Tables 5 and S6
comparisons <- list(c("Realistic", "Baseline"), c("Optimal", "Baseline"), c("Optimal", "Realistic"),
                    c("Congested", "Optimal"), c("Congested", "Baseline"), c("Control1", "Baseline"), c("Control2", "Baseline"))
if (run("bootstrap")) {
  parts <- if (file.exists(cache("bootstrap"))) readRDS(cache("bootstrap")) else list()
  set.seed(99); cell_sample <- sample(nrow(grid_df), 20000)   # per-cell intervals on a random cell sample
  todo <- setdiff(seq_len(B_BOOT), as.integer(names(parts)))
  for (chunk in split(todo, ceiling(seq_along(todo) / NCORE))) {
    res <- mclapply(chunk, function(b) {
      set.seed(1000 + b)
      idx <- c(sample(which(model_data$y == 1), replace = TRUE), sample(which(model_data$y == 0), replace = TRUE))
      fb <- fit_model(model_data[idx, ]); pb <- lapply(scen_X, function(d) prob(fb, d))
      list(axes = sapply(comparisons, function(cc) mean((pb[[cc[1]]] - pb[[cc[2]]])[principal_axes])),
           p05  = sapply(comparisons, function(cc) quantile(pb[[cc[1]]] - pb[[cc[2]]], 0.05)),
           cells = sapply(comparisons, function(cc) (pb[[cc[1]]] - pb[[cc[2]]])[cell_sample]))
    }, mc.cores = NCORE)
    names(res) <- chunk; parts <- c(parts, res); save_cache(parts, "bootstrap"); log_msg("bootstrap ", max(chunk))
  }
  AX <- sapply(parts, `[[`, "axes"); P5 <- sapply(parts, `[[`, "p05")
  t6 <- do.call(rbind, lapply(seq_along(comparisons), function(k) {
    C <- sapply(parts, function(r) r$cells[, k]); lo <- apply(C, 1, quantile, .025); hi <- apply(C, 1, quantile, .975)
    data.frame(comparison = paste(scenario_labels[comparisons[[k]]], collapse = " – "),
               mean_diff_principal_axes_pp = 100 * median(AX[k, ]),
               ci_low = 100 * quantile(AX[k, ], .025), ci_high = 100 * quantile(AX[k, ], .975),
               p05_cell_diff_pp = 100 * median(P5[k, ]),
               p05_ci_low = 100 * quantile(P5[k, ], .025), p05_ci_high = 100 * quantile(P5[k, ], .975),
               pct_cells_ci_excludes_zero = 100 * mean(lo > 0 | hi < 0), mean_cell_ci_width_pp = 100 * mean(hi - lo)) }))
  write.csv(t6, file.path(tab_dir, "table_S6_scenario_uncertainty.csv"), row.names = FALSE)
  key("bootstrap_replicates", length(parts), "bootstrap replicates")
}

# =============================================================================
# maps: summaries of the scenario maps (Section 4.2)
if (run("maps")) {
  shown <- c("Baseline", "Optimal", "Congested", "Realistic")
  P <- setNames(lapply(shown, function(s) prob(fit_ref, scen_X[[s]])), shown)
  pairs <- list(c("Realistic", "Baseline"), c("Optimal", "Realistic"), c("Optimal", "Baseline"), c("Congested", "Optimal"))
  diff_summary <- do.call(rbind, lapply(pairs, function(pp) { d <- 100 * (P[[pp[1]]] - P[[pp[2]]])
    data.frame(comparison = paste(pp[1:2], collapse = " – "), min = min(d), p05 = quantile(d, .05), median = median(d),
               p95 = quantile(d, .95), max = max(d)) }))
  write.csv(diff_summary, file.path(tab_dir, "scenario_difference_summary.csv"), row.names = FALSE)
  write.csv(data.frame(scenario = shown, points_ge_0.8 = sapply(shown, function(n) sum(P[[n]] >= 0.8))),
            file.path(tab_dir, "scenario_hotspot_counts.csv"), row.names = FALSE)
}
log_msg("done: ", paste(stages, collapse = ", "))
