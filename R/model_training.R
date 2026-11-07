# =============================================================================
# model_training.R
# Hold-out evaluation of GB, RF and XGB (Table 3) and RMSE of the predicted
# probability surfaces for different prediction-grid configurations (Table 2).
#
# Usage (from the R/ directory):  Rscript model_training.R
# Outputs: ../results/tables/table_2_rmse.csv, table_3_holdout_classification.csv
# =============================================================================

source("setup.R")
suppressPackageStartupMessages({ library(gbm); library(randomForest); library(xgboost); library(pROC); library(boot); library(parallel) })
res_dir <- Sys.getenv("BAM_RESULTS_DIR", file.path("..", "results"))
tab_dir <- file.path(res_dir, "tables"); dir.create(tab_dir, recursive = TRUE, showWarnings = FALSE)
log_msg <- function(...) message(format(Sys.time(), "%H:%M:%S"), "  ", ...)
model_names <- c(gbm = "GB", rf = "RF", xgbTree = "XGB")
grid_of <- function(m) switch(m, gbm = paramGrid_gbm, rf = paramGrid_rf, xgbTree = paramGrid_xgb)
train_model <- function(dat, method, ctrl) {
  a <- list(model_formula, data = dat, method = method, trControl = ctrl, tuneGrid = grid_of(method), metric = "ROC")
  if (method == "rf") a$ntree <- 500
  if (method == "gbm") a$verbose <- FALSE
  do.call(train, a)
}
ctrl_cv <- trainControl(method = "cv", number = 10, summaryFunction = twoClassSummary, classProbs = TRUE)

# --- Hold-out evaluation (90/10 split, 10-fold CV on the training data) -------------
set.seed(123)
idx <- createDataPartition(model_data$accidntFactorClass, p = 0.9, list = FALSE)
train_set <- model_data[idx, ]; test_set <- model_data[-idx, ]
t3 <- do.call(rbind, lapply(names(model_names), function(m) {
  set.seed(123); fit <- train_model(train_set, m, ctrl_cv)
  p <- predict(fit, newdata = test_set, type = "prob")[, "Class1"]
  pred <- factor(ifelse(p >= 0.5, "Class1", "Class0"), levels = c("Class0", "Class1"))
  cm <- confusionMatrix(pred, test_set$accidntFactorClass, positive = "Class1")
  log_msg("hold-out ", m)
  data.frame(model = model_names[[m]], cv_auc_train = max(fit$results$ROC),
             auc_test = as.numeric(auc(roc(test_set$y, p, quiet = TRUE, direction = "<"))),
             recall = cm$byClass[["Recall"]], precision = cm$byClass[["Precision"]], f1 = cm$byClass[["F1"]],
             kappa = cm$overall[["Kappa"]], accuracy = cm$overall[["Accuracy"]],
             tn = cm$table["Class0", "Class0"], fp = cm$table["Class1", "Class0"],
             fn = cm$table["Class0", "Class1"], tp = cm$table["Class1", "Class1"])
}))
write.csv(t3, file.path(tab_dir, "table_3_holdout_classification.csv"), row.names = FALSE)

# --- RMSE of predicted probability surfaces (Table 2) --------------------------------
# Final models are trained on the full dataset. Predicted probabilities are averaged
# per 25 x 25 m cell and compared with observed accident occurrence (1 = at least one
# accident in the cell); RMSE is bootstrapped over cells (R = 1000).
ctrl_none <- trainControl(method = "none", classProbs = TRUE)
set.seed(11); fits <- lapply(setNames(names(model_names), names(model_names)), function(m) train_model(model_data, m, ctrl_none))
acc_xy <- st_coordinates(accidents)
rmse_boot <- function(pred, obs) {
  b <- boot(data.frame(pred, obs), function(d, i) sqrt(mean((d$pred[i] - d$obs[i])^2)), R = 1000)
  c(rmse = b$t0, ci_low = unname(quantile(b$t, .025)), ci_high = unname(quantile(b$t, .975)))
}
t2 <- do.call(rbind, lapply(seq_len(nrow(grid_configs)), function(k) {
  cfg <- grid_configs[k, ]; g <- load_grid(cfg$points, cfg$buffer, cfg$cell); gxy <- st_coordinates(g)
  cell <- paste(floor(gxy[, 1] / 25), floor(gxy[, 2] / 25))
  acc_cells <- unique(paste(floor(acc_xy[, 1] / 25), floor(acc_xy[, 2] / 25)))
  obs <- as.integer(unique(cell) %in% acc_cells)
  out <- lapply(names(model_names), function(m) {
    p_points <- predict(fits[[m]], newdata = st_drop_geometry(g), type = "prob")[, "Class1"]
    p <- tapply(p_points, cell, mean)[unique(cell)]
    set.seed(1); r <- rmse_boot(as.numeric(p), obs)
    setNames(r, paste(model_names[[m]], names(r), sep = "_")) })
  log_msg("RMSE grid ", k)
  data.frame(points = cfg$points, buffer_m = cfg$buffer, cell_m = cfg$cell, n_points = nrow(g), n_cells = length(obs), t(unlist(out)))
}))
write.csv(t2, file.path(tab_dir, "table_2_rmse.csv"), row.names = FALSE)
log_msg("done")
