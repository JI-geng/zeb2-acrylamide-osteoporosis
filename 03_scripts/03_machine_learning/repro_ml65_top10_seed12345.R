# Reproduction of the 10-learner feature-selection panel (manuscript Step 4).
#
# LEARNER PANEL -- this list must stay identical to
#   (a) the 10 files in 02_derived_data/per_learner_rankings/
#       C5.0, DecisionTree, GBM, GLM, KNN, LASSO, NNET, RF, SVM_RFE, XGBoost
#   (b) the "Learners (10)" row of 00_parameters/WORKFLOW_PARAMETERS.md
#       LASSO, RF, SVM-RFE, XGBoost, DT, GBM, GLM, NNET, KNN, C5.0
#
# Do NOT add learners here. The consensus count reported in the manuscript
# (ZEB2 9/10; ASNS, B3GALT2, MST1, NAP1L2 7/10; PLCB3 5/10) is defined against a
# denominator of exactly 10; any extra learner silently changes that denominator
# and makes the deposited script contradict the manuscript.
#
# KNOWN SUBSTITUTION (stated in the manuscript's limitations): the platform's
# "SVM_RFE" learner wraps recursive feature elimination around an SVM. caret has
# no single-method equivalent, so svmRadial is used below as the closest
# open-source surrogate. C5.0 additionally requires the C50 package; if it is
# not installed that learner is skipped and the script warns, because the
# consensus denominator is then 9, not 10.

# Inputs are rebuilt by running the upstream step (02_derived_data is not
# distributed with this repository -- see README section 1).
suppressMessages({
  library(caret); library(randomForest); library(xgboost); library(gbm)
})
suppressMessages(try(library(C50), silent = TRUE))

# locate this script and run from its own directory, so the paths below work
# (R identifiers cannot start with "_", hence the fa0 / f0 names)
fa0 <- commandArgs(trailingOnly = FALSE)
f0  <- sub("^--file=", "", fa0[grep("^--file=", fa0)])
if (length(f0)) setwd(dirname(normalizePath(f0[1])))
d <- read.csv("../../02_derived_data/ml_split_train_seed12345.csv", row.names = 1, check.names = FALSE)
d$Group <- factor(d$Group, levels = c("Control", "DN"))
cat("train set:", nrow(d), "x", ncol(d) - 1, "\n")
cat("classes:", paste(names(table(d$Group)), table(d$Group), sep = "=", collapse = "  "), "\n\n")

HUBS <- c("ZEB2", "NAP1L2", "ASNS", "MST1", "PLCB3")
TOP <- 10

# display name -> caret method  (exactly 10 learners)
algos <- list(
  LASSO      = "glmnet",
  RF         = "rf",
  `SVM-RFE`  = "svmRadial",
  XGBoost    = "xgbTree",
  DT         = "rpart",
  GBM        = "gbm",
  GLM        = "glm",
  NNET       = "nnet",
  KNN        = "knn",
  `C5.0`     = "C5.0"
)

ctrl <- trainControl(method = "repeatedcv", number = 5, repeats = 1,
                     classProbs = TRUE, savePredictions = "final",
                     summaryFunction = twoClassSummary)
ctrl$summaryFunction <- twoClassSummary

res <- list()
for (nm in names(algos)) {
  impl <- algos[[nm]]
  set.seed(12345)
  fit <- tryCatch(
    caret::train(Group ~ ., data = d, method = impl, trControl = ctrl,
                 metric = "ROC", tuneLength = 3, verbose = FALSE),
    error = function(e) {
      cat(sprintf("  [%s] failed: %s\n", nm, conditionMessage(e))); NULL
    }
  )
  if (is.null(fit)) next
  vi <- tryCatch(caret::varImp(fit, scale = TRUE), error = function(e) NULL)
  if (is.null(vi)) { cat(sprintf("  [%s] has no variable importance\n", nm)); next }
  imp <- vi$importance
  sc <- if (ncol(imp) == 1) imp[[1]] else apply(imp, 1, max)
  names(sc) <- rownames(imp)
  sc <- sc[is.finite(sc)]
  ord <- names(sort(sc, decreasing = TRUE))
  top <- head(ord, TOP)
  res[[nm]] <- top
  hit <- intersect(top, HUBS)
  cat(sprintf("%-10s top10 = %s\n           hub %d/5   missing: %s\n", nm,
              paste(top, collapse = ", "), length(hit),
              paste(setdiff(HUBS, top), collapse = ",")))
}

if (length(res) != 10) {
  cat("\nWARNING: only", length(res), "of 10 learners produced a ranking.",
      "The consensus denominator below is", length(res),
      "-- it is NOT comparable with the manuscript's /10 counts.\n")
}

cat("\n================ three-ensemble intersection ================\n")
# the panel reported in the manuscript: intersection of the three ensemble learners
tri <- Reduce(intersect, list(res[["RF"]], res[["XGBoost"]], res[["GBM"]]))
cat("RF intersect XGBoost intersect GBM =", length(tri), ":", paste(sort(tri), collapse = ", "), "\n")
cat("identical to the 5 hub genes?", identical(sort(tri), sort(HUBS)), "\n\n")

cat("consensus count (out of", length(res), "learners):\n")
allg <- unique(unlist(res))
cnt <- sapply(allg, function(g) sum(sapply(res, function(t) g %in% t)))
cnt <- sort(cnt, decreasing = TRUE)
for (i in seq_len(min(15, length(cnt)))) {
  g <- names(cnt)[i]
  cat(sprintf("  %-12s %d/%d %s\n", g, cnt[i], length(res),
              ifelse(g %in% HUBS, " <= hub", "")))
}
cat("\nthe 5 hub genes:\n")
for (g in HUBS) cat(sprintf("  %-8s %d/%d\n", g,
                            ifelse(g %in% names(cnt), cnt[g], 0), length(res)))

write.csv(do.call(rbind, lapply(names(res), function(n)
  data.frame(Algorithm = n, Top10 = paste(res[[n]], collapse = ";"),
             stringsAsFactors = FALSE))),
  "rerun_top10.csv", row.names = FALSE)
cat("\nwritten: rerun_top10.csv\n")
