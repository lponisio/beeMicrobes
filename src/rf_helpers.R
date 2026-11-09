# Helper functions for bee_rf.R
# needs: dplyr, tibble, ranger, pROC, randomForestSRC (loaded in main script)

# merge categories with < min_n bees into "other"
# (rare levels can't be learned and break site-held-out CV)
lump <- function(vals, min_n = 30) {
  vals <- as.character(vals)
  factor(if_else(vals %in% names(which(table(vals) < min_n)), "other", vals))
}

# out-of-fold predictions holding out whole groups (site folds or species);
# each bee is predicted by a forest that never saw its group to avoid within group correlations.
# factor resp -> class probabilities (one column per class); numeric resp -> vector
cv_pred <- function(resp, preds, groups, wts = NULL) {
  oof <- if (is.factor(resp)) matrix(NA, length(resp), nlevels(resp), dimnames = list(NULL, levels(resp))) else rep(NA, length(resp))
  for (grp in unique(groups)) {
    test <- groups == grp                                     # bees in the held-out group
    fit <- ranger(x = preds[!test, , drop = FALSE], y = resp[!test], case.weights = wts[!test],
                  probability = is.factor(resp), num.trees = 500,
                  respect.unordered.factors = "order")
    pred_te <- predict(fit, preds[test, , drop = FALSE])$predictions
    if (is.factor(resp)) oof[test, colnames(pred_te)] <- pred_te else oof[test] <- pred_te
  }
  oof
}

# weight = 1 / class size, so detected and not-detected bees count equally
# (stops forests for rare pathogens from always predicting "not detected")
bal_wts <- function(resp) as.numeric(1 / table(resp)[as.character(resp)])

# group-held-out cross-validated AUC for a 0/1 outcome, using balanced weights
cv_auc <- function(resp, preds, groups) {
  resp <- factor(resp)
  as.numeric(auc(resp, cv_pred(resp, preds, groups, wts = bal_wts(resp))[, "1"], quiet = TRUE))
}

# ranger's named importance vector -> two-column tibble
imp_tbl <- function(fit) tibble(variable = names(fit$variable.importance), importance = fit$variable.importance)

# fit an rfsrc forest for a 0/1 outcome, then refit on n_sub subsamples
# to get confidence intervals on each variable's permutation importance.
# rfq = TRUE: quantile classifier, handles rare pathogens (imbalanced classes)
path_imp <- function(resp, preds, n_sub = 50) {
  rf_dat  <- as.data.frame(bind_cols(y = factor(resp), preds))
  fit     <- rfsrc(y ~ ., data = rf_dat, ntree = 500, rfq = TRUE, perf.type = "gmean",
                   importance = "permute")
  sub_fit <- subsample(fit, B = n_sub, verbose = FALSE)        # slow step
  imp_ci  <- extract.subsample(sub_fit)$var.jk.sel.Z |>
    as_tibble(rownames = "variable") |>
    arrange(desc(mean))
  list(fit = fit, sub_fit = sub_fit, imp_ci = imp_ci)
}

# PIRF returns 3 fits (mtry = 10% / sqrt / log2 of ASVs); keep the one with lowest OOB error
best_mtry <- function(pirf_out) pirf_out[[which.min(sapply(pirf_out, `[[`, "oob.err"))]]
