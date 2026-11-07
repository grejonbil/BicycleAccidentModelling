# =============================================================================
# install_packages.R
# Installs the R packages required by the scripts in R/. The versions listed were used to produce
# the results in results/ (R 4.3.3).
#
# Note: caret (7.0.1) is not compatible with xgboost >= 2.0; xgboost 1.7.x is
# therefore installed from the CRAN archive.
# =============================================================================

versions <- c(
  sf = "1.1.3", dplyr = "1.2.1", tidyr = "1.3.2", caret = "7.0.1", gbm = "2.3.1",
  randomForest = "4.7.1.2", pROC = "1.19.1", ggplot2 = "4.0.3", pdp = "0.10.0",
  kernelshap = "0.9.1", shapviz = "0.10.4", patchwork = "1.3.2", igraph = "2.3.4",
  scico = "1.5.0", rmarkdown = "2.32", knitr = "1.52", data.table = "1.18.6.1", raster = "3.6.32", boot = "1.3.30",
  future = "1.80.0", future.apply = "1.21.0", ggthemes = "7.0.0",
  viridis = "0.6.5", gridExtra = "2.3.1"
)

missing <- setdiff(names(versions), rownames(installed.packages()))
if (length(missing)) install.packages(missing)

if (!"xgboost" %in% rownames(installed.packages()) || packageVersion("xgboost") >= "2.0.0") {
  install.packages("https://cran.r-project.org/src/contrib/Archive/xgboost/xgboost_1.7.7.1.tar.gz",
                   repos = NULL, type = "source")
}

# Report installed versions
for (p in c(names(versions), "xgboost")) message(sprintf("%-14s %s", p, as.character(packageVersion(p))))
