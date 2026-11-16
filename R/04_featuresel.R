# =============================================================================
#  R/04_featuresel.R  --  feature selection (shared)
#
#  Three learners over ONE shared design matrix.  Replaces 05ML.R plus the
#  standalone lasso.R / XGBoost.R copies in UC-3 (same inputs, three blocks).
#  The published 9-gene intersection is RF (n=19) n LASSO n XGBoost.
# =============================================================================

get_candidate_matrix <- function() {
  expr  <- read_expr(file.path(cohort_dir(TRAIN_COHORT), "expr_normalized.tsv"))
  deg   <- trimws(readLines(file.path(DIR_TABLES, "deg", "DEG_main.tsv"), warn = FALSE))[-1]
  pool  <- gene_set_intersection_pool()
  genes <- intersect(rownames(expr), intersect(deg, pool))
  grp   <- read_group(file.path(cohort_dir(TRAIN_COHORT), "group.tsv"))
  common <- intersect(colnames(expr), grp$Sample)
  x <- t(expr[genes, common, drop = FALSE])
  y <- factor(grp$Type[match(common, grp$Sample)], levels = c("control", "UC"))
  list(x = x, y = y)
}

gene_set_intersection_pool <- function() {
  f <- file.path(DIR_TABLES, "gene_sets", "ARG_DEG_intersection.txt")
  if (file.exists(f)) trimws(readLines(f, warn = FALSE)) else character(0)
}

run_lasso <- function(dat) {
  set.seed(SEED$random_forest)
  suppressPackageStartupMessages(library(glmnet))
  fit <- glmnet::cv.glmnet(dat$x, dat$y, family = "binomial",
                           alpha = FEATSEL$lasso$alpha,
                           nfolds = FEATSEL$lasso$nfolds,
                           type.measure = "class")
  cf <- stats::coef(fit, s = FEATSEL$lasso$lambda_rule)
  g  <- rownames(cf)[which(cf != 0)]; g <- g[!g %in% "(Intercept)"]
  list(fit = fit, genes = unique(g),
       coef = data.frame(gene = rownames(cf), coefficient = as.numeric(cf)))
}

run_random_forest <- function(dat) {
  set.seed(SEED$random_forest)
  suppressPackageStartupMessages(library(randomForest))
  rf <- randomForest::randomForest(dat$x, dat$y, ntree = FEATSEL$random_forest$ntree,
                                   importance = TRUE)
  imp <- randomForest::importance(rf)
  g <- rownames(imp)[imp[, "MeanDecreaseGini"] > THRESH$rf_gini]
  list(fit = rf, genes = unique(g),
       importance = data.frame(gene = rownames(imp), gin = imp[, "MeanDecreaseGini"]))
}

run_xgboost <- function(dat) {
  set.seed(SEED$xgboost)
  suppressPackageStartupMessages(library(xgboost))
  y <- as.integer(dat$y) - 1L
  dtrain <- xgboost::xgb.DMatrix(dat$x, label = y)
  pars <- list(objective = "binary:logistic", eval_metric = "auc",
               max_depth = FEATSEL$xgboost$max_depth, eta = FEATSEL$xgboost$eta,
               subsample = FEATSEL$xgboost$subsample,
               colsample_bytree = FEATSEL$xgboost$colsample_bytree)
  fit <- xgboost::xgb.train(pars, dtrain, nrounds = FEATSEL$xgboost$nrounds,
                            nfold = FEATSEL$xgboost$nfold)
  imp <- xgboost::xgb.importance(feature_names = colnames(dat$x), model = fit)
  g <- head(imp$Feature[order(-imp$Gain)], FEATSEL$xgboost_top_n)
  list(fit = fit, genes = unique(g), importance = as.data.frame(imp))
}
