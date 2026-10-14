# Inputs are rebuilt by running the upstream step.
suppressMessages({
  library(caret); library(randomForest); library(gbm)
})
# locate this script and run from its own directory, so the paths below work
# (R identifiers cannot start with "_", hence the fa0 / f0 names)
fa0 <- commandArgs(trailingOnly = FALSE)
f0  <- sub("^--file=", "", fa0[grep("^--file=", fa0)])
if (length(f0)) setwd(dirname(normalizePath(f0[1])))
expr <- read.csv("../../02_derived_data/expression_matrix_all_12487genes_106samples.csv", row.names = 1, check.names = FALSE)
d65  <- read.csv("../../02_derived_data/ml_split_train_seed12345.csv", row.names = 1, check.names = FALSE)
gl   <- setdiff(colnames(d65), "Group")
expr169 <- expr[gl, ]
grp <- read.csv("../../02_derived_data/group_labels_106samples.csv", stringsAsFactors = FALSE)
g <- factor(grp$Group)
case_lv <- "low BMD"; ctrl_lv <- "high BMD"

# （）
REF_RF  <- c("NAP1L2","ZEB2","B3GALT2","ASNS","PLCB3","TUBB2B","FAM169A","MST1","SLC4A5","AKT1")
REF_GBM <- c("NAP1L2","ZEB2","MST1","PLCB3","ASNS","ADAM11","ABCC10","IL17RA","GNGT1","GSN")
HUBS <- c("ZEB2","NAP1L2","ASNS","MST1","PLCB3")

run_one <- function(seed, label) {
  set.seed(seed)
  idx <- as.integer(caret::createDataPartition(y = g, p = 0.7, list = FALSE))
  dt <- as.data.frame(t(expr169[, idx]), check.names = FALSE)
  y  <- factor(ifelse(g[idx] == case_lv, "DN", "Control"), levels = c("Control", "DN"))
  cat("\n==============", label, "(seed =", seed, ", n =", nrow(dt), ") ==============\n")

  # ---- RF：randomForest(ntree=1000) + MeanDecreaseGini ----
  set.seed(410 + 1002)
  rf <- randomForest::randomForest(x = dt, y = y, ntree = 1000)
  imp_rf <- randomForest::importance(rf)[, "MeanDecreaseGini"]
  top_rf <- head(names(sort(imp_rf, decreasing = TRUE)), 10)

  # ---- GBM：caret gbm +  + varImp(scale=F)----
  set.seed(410 + 1004)
  fitControl <- caret::trainControl(method = "repeatedcv", number = 5, repeats = 1,
                                    verboseIter = FALSE, allowParallel = FALSE)
  gbmGrid <- expand.grid(interaction.depth = c(1, 3), n.trees = c(50, 100),
                         shrinkage = 0.1, n.minobsinnode = 10)
  fit_gbm <- caret::train(x = dt, y = y, method = "gbm", trControl = fitControl,
                          tuneGrid = gbmGrid, verbose = FALSE)
  imp_gbm <- caret::varImp(fit_gbm, scale = FALSE)$importance
  top_gbm <- head(rownames(imp_gbm)[order(-imp_gbm[[1]])], 10)

  cat("RF  top10:", paste(top_rf, collapse = ", "), "\n")
  cat("    hub", length(intersect(top_rf, HUBS)), "/5    RF ？", setequal(top_rf, REF_RF), "\n")
  cat("GBM top10:", paste(top_gbm, collapse = ", "), "\n")
  cat("    hub", length(intersect(top_gbm, HUBS)), "/5    GBM ？", setequal(top_gbm, REF_GBM), "\n")
  cat("    RF ∩ GBM =", paste(sort(intersect(top_rf, top_gbm)), collapse = ", "), "\n")
  invisible(list(rf = top_rf, gbm = top_gbm))
}

cat("：\n  RF :", paste(REF_RF, collapse = ", "), "\n  GBM:", paste(REF_GBM, collapse = ", "), "\n")
r410   <- run_one(410,   "6.1 ")
r12345 <- run_one(12345, "6.5 ")
