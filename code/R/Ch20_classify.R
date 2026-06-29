# =====================================================================
# Chapter 20 - A Taste of Machine Learning (R)
# Classify tissue (leaf / root / stem) from 20 gene-expression features.
#
# Runs out of the box: it uses tidymodels if you have it, otherwise it falls
# back to base R + nnet + class (both ship WITH R), so no big install is needed.
# If nnet/class are somehow missing they are fetched into a writable user library.
#
# Data: expression.csv (150 samples x 20 genes + tissue label).
# Verified: test accuracy ~0.93, 5-fold CV ~0.95, markers = gene01..gene06.
# =====================================================================

# ---------- 0. helper: install a package into a writable personal library ----
ensure_packages <- function(pkgs, repo = "https://cloud.r-project.org") {
  user_lib <- Sys.getenv("R_LIBS_USER")
  if (!nzchar(user_lib)) user_lib <- file.path(path.expand("~"), "R-libs")
  if (!dir.exists(user_lib)) dir.create(user_lib, recursive = TRUE, showWarnings = FALSE)
  .libPaths(c(user_lib, .libPaths()))
  for (p in pkgs) if (!requireNamespace(p, quietly = TRUE)) {
    message("Installing '", p, "' into ", user_lib, " ...")
    try(install.packages(p, lib = user_lib, repos = repo), silent = TRUE)
  }
}

set.seed(0)
df <- read.csv("expression.csv")
df$tissue <- factor(df$tissue)
gcols <- grep("^gene", names(df), value = TRUE)
df <- df[, c(gcols, "tissue")]
cat("samples:", nrow(df), " genes:", length(gcols), "\n"); print(table(df$tissue))

USE_TM <- requireNamespace("tidymodels", quietly = TRUE)

# =====================================================================
# PATH A — tidymodels (canonical), used only if it is already installed
# =====================================================================
if (USE_TM) {
  suppressMessages(library(tidymodels))
  split <- initial_split(df, prop = 0.70, strata = tissue)
  train <- training(split); test <- testing(split)
  cat("\n[tidymodels]  train:", nrow(train), " test:", nrow(test), "\n")
  rec  <- recipe(tissue ~ ., data = train) |> step_normalize(all_numeric_predictors())
  spec <- multinom_reg(penalty = 0) |> set_engine("nnet")
  wf   <- workflow() |> add_recipe(rec) |> add_model(spec)
  fit  <- fit(wf, data = train)
  pred <- predict(fit, test) |> bind_cols(test["tissue"])
  cat("test accuracy:", round(accuracy(pred, tissue, .pred_class)$.estimate, 3), "\n")
  print(conf_mat(pred, tissue, .pred_class))
  folds <- vfold_cv(df, v = 5, strata = tissue)
  cv <- fit_resamples(wf, resamples = folds, metrics = metric_set(accuracy))
  cat("5-fold CV accuracy:", round(collect_metrics(cv)$mean[1], 3), "\n")

# =====================================================================
# PATH B — base R fallback (no tidymodels needed)
# =====================================================================
} else {
  cat("\n(tidymodels not installed -> base-R fallback with nnet + class)\n")
  ensure_packages(c("nnet", "class"))      # both normally ship with R
  suppressMessages({ library(nnet); library(class) })

  X <- as.matrix(df[, gcols]); y <- df$tissue
  n <- nrow(X)

  # --- standardise helper: fit centre/scale on a training subset ---
  standardize <- function(tr) {
    mu <- colMeans(X[tr, ]); sdv <- apply(X[tr, ], 2, sd); sdv[sdv == 0] <- 1
    scale(X, center = mu, scale = sdv)            # returns all rows, train-fit
  }
  # --- stratified index helpers ---
  strat_test <- function(frac = 0.30, seed = 0) {
    set.seed(seed)
    unlist(lapply(split(seq_len(n), y), function(ix) sample(ix, round(length(ix) * frac))))
  }
  strat_folds <- function(k = 5, seed = 0) {
    set.seed(seed); f <- vector("list", k)
    for (lv in levels(y)) { ix <- sample(which(y == lv))
      for (i in seq_along(ix)) f[[(i - 1) %% k + 1]] <- c(f[[(i - 1) %% k + 1]], ix[i]) }
    f
  }

  # ---------- 1. train/test split ----------
  te <- sort(strat_test(0.30, 0)); tr <- setdiff(seq_len(n), te)
  Z  <- standardize(tr)
  trn <- data.frame(Z[tr, , drop = FALSE], tissue = y[tr])
  tst <- data.frame(Z[te, , drop = FALSE], tissue = y[te])
  cat("\ntrain:", length(tr), " test:", length(te), "\n")

  # ---------- 2. multinomial logistic regression ----------
  m    <- multinom(tissue ~ ., data = trn, trace = FALSE)
  pred <- predict(m, tst)
  cat("test accuracy:", round(mean(pred == tst$tissue), 3), "\n")
  cat("confusion matrix (rows = true, cols = predicted):\n")
  print(table(true = tst$tissue, pred = pred))

  # ---------- 3. helpers that return predictions, for cross-validation ----------
  predictors <- list(
    logreg = function(Xtr, ytr, Xte)
      predict(multinom(ytr ~ ., data = data.frame(Xtr, ytr = ytr), trace = FALSE),
              data.frame(Xte)),
    knn15  = function(Xtr, ytr, Xte) knn(Xtr, Xte, cl = ytr, k = 15)
  )
  if (requireNamespace("ranger", quietly = TRUE))
    predictors$forest <- function(Xtr, ytr, Xte)
      predict(ranger::ranger(ytr ~ ., data = data.frame(Xtr, ytr = ytr),
                             num.trees = 500), data.frame(Xte))$predictions

  cv_accuracy <- function(fun, k = 5) {
    folds <- strat_folds(k, 0); acc <- numeric(k)
    for (i in seq_len(k)) {
      te <- folds[[i]]; tr <- setdiff(seq_len(n), te); Z <- standardize(tr)
      acc[i] <- mean(fun(Z[tr, , drop = FALSE], y[tr], Z[te, , drop = FALSE]) == y[te])
    }
    acc
  }

  # ---------- 4. 5-fold CV for the logistic model ----------
  cv <- cv_accuracy(predictors$logreg)
  cat(sprintf("\n5-fold CV accuracy: %.3f +/- %.3f\n", mean(cv), sd(cv)))

  # ---------- 5. compare a few algorithms ----------
  cat("algorithm comparison (5-fold CV accuracy):\n")
  for (nm in names(predictors))
    cat("  ", nm, ":", round(mean(cv_accuracy(predictors[[nm]])), 3), "\n")
}

# =====================================================================
# 6. Which genes matter? one-way ANOVA F per gene (base R, both paths)
# =====================================================================
Fscores <- sapply(gcols, function(g) summary(aov(df[[g]] ~ df$tissue))[[1]][["F value"]][1])
cat("\ntop 6 genes by ANOVA F (the real markers):\n")
print(round(sort(Fscores, decreasing = TRUE)[1:6], 1))

cat("\nDone.\n")
