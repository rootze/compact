# ------------------------------------------------------------
# Package setup with strict version pinning for SHAP pipeline
# ------------------------------------------------------------

# --- STRICTLY pinned (do NOT auto-install) ---
pinned_pkgs <- list(
  xgboost = "1.7.11.1",
  SHAPforxgboost = "0.1.3"
)

# --- Regular CRAN packages ---
required_pkgs <- c("ggplot2", "data.table", "ggrastr")

# Enforce versions (fail fast, no auto-install)
for (pkg in names(pinned_pkgs)) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    stop(sprintf(
      "Required package '%s' not installed.\nPlease install version %s in this conda env.",
      pkg, pinned_pkgs[[pkg]]
    ))
  }
  installed <- as.character(packageVersion(pkg))
  if (!identical(installed, pinned_pkgs[[pkg]])) {
    stop(sprintf(
      "Package version mismatch for '%s': required %s but found %s",
      pkg, pinned_pkgs[[pkg]], installed
    ))
  }
  library(pkg, character.only = TRUE)
}

# Auto-install only the safe packages
invisible(lapply(required_pkgs, function(pkg) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg)
  }
  library(pkg, character.only = TRUE)
}))


#' Identify Key Driver Genes Using SHAP and XGBoost
#'
#' This function trains an XGBoost classifier to distinguish either:
#' \itemize{
#'   \item a specified case group against all other cells in a metadata column
#'         (\code{set_case} vs "rest"), or
#'   \item two specific groups within that column (\code{ident.1} vs \code{ident.2}, pairwise mode),
#' }
#' using gene expression data from a Seurat object.
#'
#' SHAP (SHapley Additive exPlanations) values are computed to quantify the
#' contribution of each gene to the model predictions. The top SHAP-ranked genes
#' are returned as key driver candidates.
#'
#' Three model training modes are supported via the \code{mode} parameter:
#' \itemize{
#'   \item \code{"full"}: Train on the entire dataset (no test set, no cross-validation).
#'   \item \code{"split"}: Randomly split data into training and testing subsets.
#'   \item \code{"cv"}: Perform k-fold cross-validation.
#' }
#' In \code{"split"} mode, test set AUC is computed and stored. In \code{"cv"} mode,
#' per-fold AUCs and mean AUC are computed. In \code{"full"} mode, AUC is not
#' computed to avoid overestimation.
#'
#' SHAP results and metadata are stored in \code{seurat_obj@misc$shap}.
#' Optionally, results are written to disk under \code{out_dir}, with all files
#' documented in an \code{info.txt}.
#'
#' Exactly one of the following must be supplied:
#' \itemize{
#'   \item \code{set_case}: for a "case vs rest" comparison, or
#'   \item \code{ident.1} and \code{ident.2}: for a pairwise comparison between two groups.
#' }
#'
#' @param seurat_obj A \code{Seurat} object with normalized gene expression and metadata.
#' @param conditions Column name in \code{seurat_obj@meta.data} indicating the condition labels.
#' @param set_case Character. Case label to compare against all other values in \code{conditions}
#'   (case vs rest mode). Must be \code{NULL} if \code{ident.1} and \code{ident.2} are used.
#' @param ident.1 Character. First group label in \code{conditions} for pairwise comparison
#'   (treated as the positive class, label 1). Must be used together with \code{ident.2} and
#'   \code{set_case} must be \code{NULL}.
#' @param ident.2 Character. Second group label in \code{conditions} for pairwise comparison
#'   (treated as the negative class, label 0). Must be used together with \code{ident.1} and
#'   \code{set_case} must be \code{NULL}.
#' @param top_n Integer. Number of top SHAP-ranked genes to return (default: 50).
#' @param max_depth Integer. Maximum depth of each XGBoost tree (default: 4).
#' @param eta Numeric. Learning rate (shrinkage) for XGBoost (default: 0.1).
#' @param nrounds Integer. Number of boosting rounds for XGBoost (default: 50).
#' @param nthread Integer. Number of CPU threads to use for XGBoost (default: 2).
#' @param variable_genes Optional character vector of genes to use. If \code{NULL},
#'   uses \code{VariableFeatures(seurat_obj)}.
#' @param out_dir Optional directory to save results. If \code{NULL}, results are not saved to disk.
#' @param mode Character. Training mode: \code{"full"}, \code{"split"}, or \code{"cv"}
#'   (default: \code{"full"}).
#' @param train_fraction Numeric. Fraction of cells used for training in \code{"split"} mode
#'   (default: 0.8).
#' @param nfold Integer. Number of folds for cross-validation (default: 5).
#' @param seed Integer. Random seed for reproducibility (default: 123).
#' @param expr_layer Character. Expression layer/slot to use when building the feature matrix.
#'   In Seurat v4 this is passed as \code{slot}; in Seurat v5 this is passed as \code{layer}.
#'   Common values are \code{"data"}, \code{"counts"}, or \code{"scale.data"} (default: \code{"data"}).
#'
#' @return A modified \code{Seurat} object with SHAP results stored in
#'   \code{seurat_obj@misc$shap}, including:
#' \itemize{
#'   \item \code{model}: Trained XGBoost model (only in \code{"full"} and \code{"split"} modes; not stored in \code{"cv"}).
#'   \item \code{shap_result}: Raw SHAP result object (only in \code{"full"} and \code{"split"} modes; not stored in \code{"cv"}).
#'   \item \code{shap_long}: Long-format SHAP values (per-cell, per-gene).
#'   \item \code{shap_summary}: Mean absolute SHAP values per gene.
#'   \item \code{key_drivers}: Top \code{top_n} SHAP-ranked genes.
#'   \item \code{variable_genes}: Genes used in the model.
#'   \item \code{test_auc}: AUC on held-out test set (only in \code{"split"} mode; \code{NA} otherwise).
#'   \item \code{auc_per_fold}: Vector of AUC scores per fold (only in \code{"cv"} mode).
#'   \item \code{mean_auc}: Mean AUC across folds (only in \code{"cv"} mode; \code{NA} if AUC not computed).
#'   \item \code{mode}: Training mode used.
#'   \item \code{comparison}: A list describing the comparison type and labels, with fields:
#'     \code{type} (either \code{"case_vs_rest"} or \code{"pairwise"}) and
#'     \code{set_case} or \code{ident.1}/\code{ident.2} accordingly.
#'   \item \code{params}: A list of key settings used (e.g., \code{max_depth}, \code{eta}, \code{nrounds}, \code{nthread}, \code{expr_layer}).
#'   \item \code{package_versions}: A named list of package versions used (e.g., \code{xgboost}, \code{SHAPforxgboost}).
#' }
#'
#' This implementation enforces exact package versions for reproducibility.
#' It will stop with an error unless \code{xgboost == "1.7.11.1"} and
#' \code{SHAPforxgboost == "0.1.3"} are installed in the active R library.
#'
#' @section Output Files (if \code{out_dir} is provided):
#' \itemize{
#'   \item \code{model.txt}: XGBoost model summary (not in \code{"cv"} mode).
#'   \item \code{model.rds}: Serialized XGBoost model object (not in \code{"cv"} mode).
#'   \item \code{shap_result.rds}: Output of \code{SHAPforxgboost::shap.values()} (not in \code{"cv"} mode).
#'   \item \code{shap_summary.txt}: Mean absolute SHAP value per gene.
#'   \item \code{variable_genes.txt}: List of genes used for model training.
#'   \item \code{key_drivers.txt}: Top SHAP-ranked genes.
#'   \item \code{shap_long.rds}: Full SHAP values in long format.
#'   \item \code{auc.txt}: AUC from held-out test set (only in \code{"split"} mode, if computed).
#'   \item \code{auc_per_fold.txt}: AUC per fold (only in \code{"cv"} mode, if computed).
#'   \item \code{info.txt}: Summary of saved outputs.
#' }
#'
#' @export
#'
#' @examples
#' # Case vs rest: AD vs all other diagnoses
#' seurat_obj <- FindShapKeyDriver(
#'   seurat_obj,
#'   conditions = "Diagnosis",
#'   set_case = "AD",
#'   top_n = 30
#' )
#'
#' # Case vs rest in split mode with custom output directory
#' seurat_obj <- FindShapKeyDriver(
#'   seurat_obj,
#'   conditions = "Diagnosis",
#'   set_case = "AD",
#'   mode = "split",
#'   out_dir = "results/shap/"
#' )
#'
#' # Case vs rest in 5-fold CV mode
#' seurat_obj <- FindShapKeyDriver(
#'   seurat_obj,
#'   conditions = "Diagnosis",
#'   set_case = "AD",
#'   mode = "cv",
#'   nfold = 5
#' )
#'
#' # Pairwise comparison: AD vs PiD only
#' seurat_obj <- FindShapKeyDriver(
#'   seurat_obj,
#'   conditions = "Diagnosis",
#'   ident.1 = "AD",
#'   ident.2 = "PiD",
#'   mode = "cv",
#'   nfold = 5
#' )

FindShapKeyDriver <- function(
  seurat_obj, conditions,
  set_case = NULL,
  ident.1 = NULL,
  ident.2 = NULL,
  top_n = 50,
  max_depth = 4, eta = 0.1, nrounds = 50,
  nthread = 2,
  variable_genes = NULL,
  out_dir = NULL,
  mode = c("full", "split", "cv"),
  train_fraction = 0.8,
  nfold = 5,
  seed = 123,
  expr_layer = "data"   # works as slot in v4, layer in v5
) {


  require(Seurat)
  require(xgboost)
  require(SHAPforxgboost)
  require(data.table)

  # =========================
  # HARD version requirements
  # =========================
  .require_pkg_version_exact <- function(pkg, required) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
      stop("Missing package: ", pkg, call. = FALSE)
    }
    installed <- as.character(utils::packageVersion(pkg))
    required  <- as.character(required)
    if (!identical(installed, required)) {
      stop(sprintf(
        "Package version mismatch for '%s': required %s but found %s.\nFix by installing the required version in this env.",
        pkg, required, installed
      ), call. = FALSE)
    }
    invisible(TRUE)
  }

  .require_pkg_version_exact("xgboost",       "1.7.11.1")
  .require_pkg_version_exact("SHAPforxgboost","0.1.3")

  # =========================
  # NOTE Input checks
  # =========================
  if (is.null(set_case) && (is.null(ident.1) || is.null(ident.2))) {
    stop("You must provide either `set_case` (case vs rest) OR both `ident.1` and `ident.2` (pairwise).", call. = FALSE)
  }
  if (!is.null(set_case) && (!is.null(ident.1) || !is.null(ident.2))) {
    stop("Please use either `set_case` OR (`ident.1` and `ident.2`), not both.", call. = FALSE)
  }

  has_pROC <- requireNamespace("pROC", quietly = TRUE)
  if (!has_pROC) warning("pROC not installed: AUC will not be computed.")

  mode <- match.arg(mode, choices = c("full", "split", "cv"))
  if (!is.numeric(train_fraction) || train_fraction <= 0 || train_fraction >= 1) {
    stop("`train_fraction` must be between 0 and 1 (exclusive).", call. = FALSE)
  }
  # if (!is.numeric(nfold) || nfold < 2) {
  #   stop("`nfold` must be >= 2.", call. = FALSE)
  # }

  # if (missing(set_case)) stop("`set_case` is required but was not provided.")

  # if (!(set_case %in% seurat_obj[[conditions]][, 1])) {
  #   stop(sprintf("`set_case = '%s'` not found in `%s` column of Seurat object.", set_case, conditions))
  # }

  stopifnot(inherits(seurat_obj, "Seurat"))


  if (missing(conditions) || !nzchar(conditions)) {
    stop("`conditions` is required.", call. = FALSE)
  }
  if (!conditions %in% colnames(seurat_obj[[]])) {
    stop("`conditions` column not found in seurat_obj@meta.data: ", conditions, call. = FALSE)
  }


  # =========================
  # NOTE Helper: Seurat v4/v5 GetAssayData compatibility
  # =========================
  .msg <- function(...) message(...)

  .get_expr_mat <- function(obj, features, layer = "data") {
    fn <- Seurat::GetAssayData
    use_layer <- "layer" %in% names(formals(fn))  # SeuratObject v5+

    mat <- if (use_layer) {
      Seurat::GetAssayData(obj, layer = layer)
    } else {
      Seurat::GetAssayData(obj, slot  = layer)
    }

    mat <- mat[features, , drop = FALSE]
    # xgboost wants numeric matrix; t() of sparse can be awkward => densify after subset
    X <- t(as.matrix(mat))  # cells x genes
    storage.mode(X) <- "double"
    X
  }


  # =========================
  # NOTE Other Helpers
  # =========================

  # .compute_shap_long <- function(model, X_eval) {
  #   shap_result <- SHAPforxgboost::shap.values(model, X_train = X_eval)
  #   SHAPforxgboost::shap.prep(shap_contrib = shap_result$shap_score, X_train = X_eval)
  # }
  .compute_shap <- function(model, X_eval) {
    shap_result <- SHAPforxgboost::shap.values(model, X_train = X_eval)
    shap_long <- SHAPforxgboost::shap.prep(
      shap_contrib = shap_result$shap_score,
      X_train      = X_eval
    )
    list(shap_result = shap_result, shap_long = shap_long)
  }


  # .summarize_shap <- function(shap_long_dt) {
  #   dt <- data.table::as.data.table(shap_long_dt)
  #   summ <- dt[, list(mean_abs_shap = mean(abs(value))), by = variable]
  #   summ[order(-mean_abs_shap)]
  # }
  ## NOTE ERROR
  # Error in h(simpleError(msg, call)) :
  # error in evaluating the argument 'x' in selecting a method for function 'mean': object 'value' not found

  # NOTE Notes Code below
  # Same root cause — value not found means the data.table column reference is failing inside .summarize_shap() when called from within the package namespace.
  # The fix from earlier needs to be applied more carefully. In data.table, even list() can fail if data.table's special evaluation isn't triggered. The safest fix for package code is to avoid non-standard evaluation entirely:
  .summarize_shap <- function(shap_long_dt) {
    dt <- data.table::as.data.table(shap_long_dt)
    genes <- unique(dt$variable)
    mean_abs <- vapply(genes, function(g) {
      mean(abs(dt$value[dt$variable == g]))
    }, numeric(1))
    summ <- data.table::data.table(variable = genes, mean_abs_shap = mean_abs)
    summ[order(-summ$mean_abs_shap)]
  }
  # NOTE This avoids data.table's [ non-standard evaluation entirely, using base R column access ($) instead — which is fully namespace-safe inside a package.

  

  .comparison_descriptor <- function() {
    if (!is.null(set_case)) {
      list(type = "case_vs_rest", set_case = set_case)
    } else {
      list(type = "pairwise", ident.1 = ident.1, ident.2 = ident.2)
    }
  }


  # =========================
  # NOTE Feature set
  # =========================
  if (is.null(variable_genes)) {
    .msg("No gene list supplied — using VariableFeatures from Seurat object.")
    variable_genes <- Seurat::VariableFeatures(seurat_obj)
  } else {
    .msg("Using user-supplied gene list.")
  }

  # gene existence check
  all_genes <- rownames(seurat_obj)
  missing_genes <- setdiff(variable_genes, all_genes)
  if (length(missing_genes) > 0) {
    stop(sprintf(
      "%d supplied genes not found in Seurat object (showing up to 10): %s",
      length(missing_genes), paste(head(missing_genes, 10), collapse = ", ")
    ), call. = FALSE)
  }

  # should add the check right after you define or confirm variable_genes, and
  # before SHAP computation or top_n usage.
  if (length(variable_genes) < top_n) {
    .msg(sprintf("Only %d genes provided. Setting top_n = %d.", length(variable_genes), length(variable_genes)))
    top_n <- length(variable_genes)
  }
  .msg(sprintf("%d genes selected for SHAP analysis.", length(variable_genes)))


  # =========================
  # Build X, y
  # =========================
  # X <- t(GetAssayData(seurat_obj, slot = "data")[variable_genes, ])
  # NOTE Error Message
  #   Error in t.default(GetAssayData(seurat_obj, slot = "data")[variable_genes,  :
  #   argument is not a matrix
  # This error means that the object you're trying to transpose with t() is not a matrix,
  # likely a sparse matrix or a non-matrix object like a data frame.

  # after define variable_genes and top_n
  # X <- t(as.matrix(GetAssayData(seurat_obj, slot = "data")[variable_genes, ])) # replace with the below code for version robustness
  # GetAssayData (version-robust helper): works on v4 + v5
  X <- .get_expr_mat(seurat_obj, variable_genes, layer = expr_layer)
  condition_vals_full <- seurat_obj[[conditions]][, 1]

  # X <- .get_expr_mat(seurat_obj, variable_genes, layer = "data")
  # condition_vals_full <- seurat_obj[[conditions]][, 1]

  # y <- as.numeric(seurat_obj[[conditions]] == set_case)
  # condition_vals <- seurat_obj[[conditions]][, 1]
  # rest_cases <- unique(condition_vals[condition_vals != set_case])
  # message(sprintf("Creating binary label for condition: '%s' vs rest.", set_case))
  # message(sprintf("The 'rest' includes: %s", paste(rest_cases, collapse = ", ")))
  # condition_vals_full already defined above

  # NOTE Replace the above code y <- ... block with pairwise-aware logic

  if (!is.null(set_case)) {
    # case vs rest
    if (!(set_case %in% condition_vals_full)) {
      stop(sprintf("`set_case = '%s'` not found in `%s`.", set_case, conditions), call. = FALSE)
    }
    y <- as.numeric(condition_vals_full == set_case)
    rest_cases <- unique(condition_vals_full[condition_vals_full != set_case])
    .msg(sprintf("Creating binary label for condition: '%s' (1) vs rest (0).", set_case))
    .msg(sprintf("The 'rest' includes: %s", paste(rest_cases, collapse = ", ")))

  } else {
    # ident.1 vs ident.2
    if (!all(c(ident.1, ident.2) %in% condition_vals_full)) {
      stop(sprintf("`ident.1 = '%s'` and/or `ident.2 = '%s'` not found in `%s`.",
                   ident.1, ident.2, conditions), call. = FALSE)
    }
    keep <- which(condition_vals_full %in% c(ident.1, ident.2))
    X <- X[keep, , drop = FALSE]
    y <- as.numeric(condition_vals_full[keep] == ident.1)
    .msg(sprintf("Creating binary label for '%s' (1) vs '%s' (0).", ident.1, ident.2))
  }

  if (length(unique(y)) != 2) {
    stop("After subsetting, y does not contain exactly two classes. Check labels/subsetting.", call. = FALSE)
  }

  # =========================
  # CV feasibility check
  # =========================
  if (mode == "cv") {
    if (nfold < 2) {
      stop("`nfold` must be >= 2.", call. = FALSE)
    }
    if (nrow(X) < nfold) {
      stop(sprintf(
        "Not enough cells for %d-fold CV: nrow(X) = %d.",
        nfold, nrow(X)
      ), call. = FALSE)
    }
  }

  # shared xgb params
  xgb_params <- list(
    objective = "binary:logistic",
    eval_metric = "auc",
    max_depth = max_depth,
    eta       = eta,
    nthread   = nthread
  )

  # auc <- NA       # my code fix here
  auc <- NA_real_

  # =========================
  # NOTE Mode: CV
  # =========================
  if (mode == "cv") {
    set.seed(seed)
    n <- nrow(X)
    fold_ids <- sample(rep(seq_len(nfold), length.out = n))

    all_shap_long <- vector("list", nfold)
    all_aucs <- rep(NA_real_, nfold)

    for (i in seq_len(nfold)) {
      .msg(sprintf("Running fold %d of %d...", i, nfold))
      test_idx  <- which(fold_ids == i)
      train_idx <- setdiff(seq_len(n), test_idx)

      X_train <- X[train_idx, , drop = FALSE]
      y_train <- y[train_idx]
      X_test  <- X[test_idx,  , drop = FALSE]
      y_test  <- y[test_idx]

      dtrain <- xgboost::xgb.DMatrix(data = X_train, label = y_train)
      dtest  <- xgboost::xgb.DMatrix(data = X_test)

      model <- xgboost::xgboost(
        params  = xgb_params,
        data    = dtrain,
        nrounds = nrounds,
        verbose = 0
      )

      # AUC on held-out fold (optional)
      y_pred <- predict(model, dtest)
      if (has_pROC) {
        all_aucs[i] <- as.numeric(pROC::auc(y_test, y_pred))
        .msg(sprintf("Fold %d AUC: %.3f", i, all_aucs[i]))
      }

      # shap_long <- .compute_shap_long(model, X_test)
      shap_out <- .compute_shap(model, X_test)  # or X
      shap_result <- shap_out$shap_result
      shap_long   <- shap_out$shap_long

      # SAFETY CHECK (ADD HERE)
      if (!all(c("variable", "value") %in% colnames(shap_long))) {
        stop("shap_long missing expected columns (variable, value).", call. = FALSE)
      }

      shap_long$fold <- i
      all_shap_long[[i]] <- shap_long
    }

    combined_shap <- data.table::rbindlist(all_shap_long)
    shap_summary <- .summarize_shap(combined_shap)
    key_drivers <- head(shap_summary$variable, top_n)
    auc <- if (has_pROC) mean(all_aucs, na.rm = TRUE) else NA_real_

    seurat_obj@misc$shap <- list(
      shap_long      = combined_shap,
      shap_summary   = shap_summary,
      key_drivers    = key_drivers,
      variable_genes = variable_genes,
      auc_per_fold   = all_aucs,
      mean_auc       = auc,
      mode           = mode,
      comparison     = .comparison_descriptor(),
      params         = list(max_depth = max_depth, eta = eta, nrounds = nrounds, nthread = nthread, expr_layer = expr_layer),
      package_versions = list(
        xgboost = as.character(utils::packageVersion("xgboost")),
        SHAPforxgboost = as.character(utils::packageVersion("SHAPforxgboost"))
      )
    )

    # if (mode == "cv") {
    #   set.seed(seed)
    #   n <- nrow(X)
    #   fold_ids <- sample(rep(1:nfold, length.out = n))
    #   all_shap_long <- list()
    #   # all_aucs <- numeric(nfold)
    #   all_aucs <- rep(NA_real_, nfold)
    #
    #   for (i in 1:nfold) {
    #     message(sprintf("Running fold %d of %d...", i, nfold))
    #     test_idx <- which(fold_ids == i)
    #     train_idx <- setdiff(seq_len(n), test_idx)

        # X_train <- X[train_idx, ]
        # y_train <- y[train_idx]
        # X_test  <- X[test_idx, ]
        # y_test  <- y[test_idx]

        # dtrain <- xgboost::xgb.DMatrix(data = X_train, label = y_train)
        # dtest  <- xgboost::xgb.DMatrix(data = X_test)

        # NOTE replace all calls to xgboost::xgboost() with xgboost::xgb.train()
        # params <- list(objective = "binary:logistic", max_depth = max_depth,
        #                eta = eta, nthread = nthread)
        # model <- xgboost::xgboost(params = params, data = dtrain,
        #                           nrounds = nrounds, verbose = 0)
        # params <- list(
        #   objective = "binary:logistic",
        #   max_depth = max_depth,
        #   eta = eta,
        #   nthread = nthread,
        #   eval_metric = "auc"
        # )
        #
        # model <- xgboost::xgb.train(
        #   params  = params,
        #   data    = dtrain,
        #   nrounds = nrounds,
        #   verbose = 0
        # )
      #
      #   y_pred <- predict(model, dtest)
      #   if (has_pROC) {
      #     auc_fold <- pROC::auc(y_test, y_pred)
      #     all_aucs[i] <- auc_fold
      #     message(sprintf("Fold %d AUC: %.3f", i, auc_fold))
      #   }
      #
      #   shap_result <- SHAPforxgboost::shap.values(model, X_train = X_test)
      #   shap_long <- SHAPforxgboost::shap.prep(shap_contrib = shap_result$shap_score, X_train = X_test)
      #   shap_long$fold <- i
      #   all_shap_long[[i]] <- shap_long
      # }
      #
      # combined_shap <- data.table::rbindlist(all_shap_long)
      # shap_summary <- combined_shap[, .(mean_abs_shap = mean(abs(value))), by = variable]
      # shap_summary <- shap_summary[order(-mean_abs_shap)]
      # key_drivers <- head(shap_summary$variable, top_n)
      # # auc <- mean(all_aucs)
      # # auc <- if (has_pROC) mean(all_aucs) else NA_real_
      # # Above: This already works, but if for any weird reason one fold’s AUC is NA, mean(all_aucs) would also become NA.
      # auc <- if (has_pROC) mean(all_aucs, na.rm = TRUE) else NA_real_
      #
      # seurat_obj@misc$shap <- list(
      #   shap_long = combined_shap,
      #   shap_summary = shap_summary,
      #   key_drivers = key_drivers,
      #   variable_genes = variable_genes,
      #   auc_per_fold = all_aucs,
      #   mean_auc = auc,
      #   mode = mode
      # )
      #
      # comparison <- if (!is.null(set_case)) {
      #   list(type = "case_vs_rest", set_case = set_case)
      # } else {
      #   list(type = "pairwise", ident.1 = ident.1, ident.2 = ident.2)
      # }
      #
      # seurat_obj@misc$shap$comparison <- comparison

  } else {
    # =========================
    # Mode: split / full
    # =========================
    if (mode == "split") {
      set.seed(seed)
      n <- nrow(X)
      train_idx <- sample(seq_len(n), size = floor(train_fraction * n))
      test_idx  <- setdiff(seq_len(n), train_idx)

      X_train <- X[train_idx, , drop = FALSE]
      y_train <- y[train_idx]
      X_test  <- X[test_idx,  , drop = FALSE]
      y_test  <- y[test_idx]

      .msg(sprintf("Split data into training (%d) and testing (%d) cells", length(train_idx), length(test_idx)))

      dtrain <- xgboost::xgb.DMatrix(data = X_train, label = y_train)
      dtest  <- xgboost::xgb.DMatrix(data = X_test)

      model <- xgboost::xgboost(
        params  = xgb_params,
        data    = dtrain,
        nrounds = nrounds,
        verbose = 0
      )

      y_pred <- predict(model, dtest)
      if (has_pROC) {
        auc <- as.numeric(pROC::auc(y_test, y_pred))
        .msg(sprintf("Test AUC: %.3f", auc))
      }

      # shap_long <- .compute_shap_long(model, X_test)
      # shap_result <- SHAPforxgboost::shap.values(model, X_train = X_test)
      shap_out <- .compute_shap(model, X_test)  # or X
      shap_result <- shap_out$shap_result
      shap_long   <- shap_out$shap_long

      # SAFETY CHECK (ADD HERE)
      if (!all(c("variable", "value") %in% colnames(shap_long))) {
        stop("shap_long missing expected columns (variable, value).", call. = FALSE)
      }

    } else {
      # =========================
      # Mode: Full
      # =========================
      dtrain <- xgboost::xgb.DMatrix(data = X, label = y)

      model <- xgboost::xgboost(
        params  = xgb_params,
        data    = dtrain,
        nrounds = nrounds,
        verbose = 0
      )

      # shap_long <- .compute_shap_long(model, X)
      # shap_result <- SHAPforxgboost::shap.values(model, X_train = X)
      shap_out <- .compute_shap(model, X) # Full
      shap_result <- shap_out$shap_result
      shap_long   <- shap_out$shap_long

      # SAFETY CHECK (ADD HERE)
      if (!all(c("variable", "value") %in% colnames(shap_long))) {
        stop("shap_long missing expected columns (variable, value).", call. = FALSE)
      }

    }

    shap_summary <- .summarize_shap(shap_long)
    key_drivers <- head(shap_summary$variable, top_n)

    seurat_obj@misc$shap <- list(
      model         = model,
      shap_result   = shap_result,   # keeps old behavior
      shap_long     = shap_long,
      shap_summary  = shap_summary,
      key_drivers   = key_drivers,
      variable_genes = variable_genes,
      test_auc      = if (mode == "split") auc else NA_real_,
      mode          = mode,
      comparison    = .comparison_descriptor(),
      params         = list(max_depth = max_depth, eta = eta, nrounds = nrounds, nthread = nthread, expr_layer = expr_layer),
      package_versions = list(
        xgboost = as.character(utils::packageVersion("xgboost")),
        SHAPforxgboost = as.character(utils::packageVersion("SHAPforxgboost"))
      )
    )
  }

  # =========================
  # Save outputs
  # =========================
  if (!is.null(out_dir)) {
    shap_dir <- file.path(out_dir, paste0("shap_output_", mode))
    if (!dir.exists(shap_dir)) dir.create(shap_dir, recursive = TRUE)
    .msg(sprintf("Saving SHAP outputs to: %s", shap_dir))

    cmp_line <- {
      cmp <- .comparison_descriptor()
      paste(sprintf("%s=%s", names(cmp), unlist(cmp)), collapse = " | ")
    }

    info_lines <- c(
      sprintf("SHAP output summary (mode = '%s'):", mode),
      sprintf("conditions: %s", conditions),
      # sprintf("comparison: %s", paste(unlist(.comparison_descriptor()), collapse = " | ")),
      sprintf("comparison: %s", cmp_line),
      sprintf("n_cells: %d", nrow(X)),
      sprintf("n_features: %d", length(variable_genes)),
      sprintf("xgboost: %s", as.character(utils::packageVersion("xgboost"))),
      sprintf("SHAPforxgboost: %s", as.character(utils::packageVersion("SHAPforxgboost"))),
      ""
    )

    if (mode %in% c("full", "split")) {
      capture.output(str(seurat_obj@misc$shap$model), file = file.path(shap_dir, "model.txt"))
      saveRDS(seurat_obj@misc$shap$model, file = file.path(shap_dir, "model.rds"))
      saveRDS(seurat_obj@misc$shap$shap_result, file = file.path(shap_dir, "shap_result.rds"))
      info_lines <- c(info_lines,
        "- model.txt: XGBoost model structure",
        "- model.rds: Serialized model",
        "- shap_result.rds: SHAP decomposition object"
      )
    }

    utils::write.table(seurat_obj@misc$shap$shap_summary,
      file = file.path(shap_dir, "shap_summary.txt"),
      sep = "\t", quote = FALSE, row.names = FALSE
    )
    utils::write.table(data.frame(variable_genes = seurat_obj@misc$shap$variable_genes),
      file = file.path(shap_dir, "variable_genes.txt"),
      sep = "\t", quote = FALSE, row.names = FALSE
    )
    utils::write.table(data.frame(key_drivers = seurat_obj@misc$shap$key_drivers),
      file = file.path(shap_dir, "key_drivers.txt"),
      sep = "\t", quote = FALSE, row.names = FALSE
    )
    saveRDS(seurat_obj@misc$shap$shap_long, file = file.path(shap_dir, "shap_long.rds"))

    info_lines <- c(info_lines,
      "- shap_summary.txt: Mean absolute SHAP per gene",
      "- variable_genes.txt: Feature list",
      "- key_drivers.txt: Top SHAP genes",
      "- shap_long.rds: Long SHAP values"
    )

    if (mode == "cv") {
      utils::write.table(data.frame(auc_per_fold = seurat_obj@misc$shap$auc_per_fold),
        file = file.path(shap_dir, "auc_per_fold.txt"),
        sep = "\t", quote = FALSE, row.names = FALSE
      )
      info_lines <- c(info_lines, "- auc_per_fold.txt: AUC per fold")
    }

    if (mode == "split") {
      writeLines(sprintf("%.5f", seurat_obj@misc$shap$test_auc), file.path(shap_dir, "auc.txt"))
      info_lines <- c(info_lines, "- auc.txt: Test AUC")
    }

    writeLines(info_lines, con = file.path(shap_dir, "info.txt"))
    .msg("SHAP results saved. See info.txt for details.")
  }

  .msg(sprintf("SHAP driver gene analysis complete (mode = '%s')", mode))
  # seurat_obj
  return(seurat_obj)

}



#
#' Reload SHAP Output from Disk
#'
#' This function reloads previously saved SHAP analysis results from a specified directory
#' (e.g., created by \code{FindShapKeyDriver()}) and attaches them to a Seurat object.
#'
#' It is useful when SHAP results were not stored in \code{@misc} or when working across sessions.
#'
#' @param seurat_obj A \code{Seurat} object to attach SHAP results to.
#' @param shap_dir Path to the directory containing SHAP outputs (e.g., "results/shap_output_split").
#'
#' @return A \code{Seurat} object with SHAP results loaded into \code{seurat_obj@misc$shap}.
#' @export
#'
#' @examples
#' seurat_obj <- ReloadShapOutput(seurat_obj, shap_dir = "results/shap_output_split")

ReloadShapOutput <- function(seurat_obj, shap_dir, enforce_pins = FALSE) {

  stopifnot(inherits(seurat_obj, "Seurat"))
  if (!dir.exists(shap_dir)) stop(sprintf("Directory not found: %s", shap_dir), call. = FALSE)

  message(sprintf("Reloading SHAP outputs from: %s", shap_dir))
  shap_list <- list()

  # -------------------------
  # Optional strict version checks (matches FindShapKeyDriver pins)
  # -------------------------
  if (isTRUE(enforce_pins)) {
    pinned <- list(xgboost = "1.7.11.1", SHAPforxgboost = "0.1.3")
    for (pkg in names(pinned)) {
      if (!requireNamespace(pkg, quietly = TRUE)) {
        stop(sprintf("Missing package '%s' (required %s).", pkg, pinned[[pkg]]), call. = FALSE)
      }
      inst <- as.character(utils::packageVersion(pkg))
      if (!identical(inst, pinned[[pkg]])) {
        stop(sprintf("Package version mismatch for '%s': required %s but found %s", pkg, pinned[[pkg]], inst),
             call. = FALSE)
      }
    }
  }

  # -------------------------
  # Load files if present
  # -------------------------
  f_summary <- file.path(shap_dir, "shap_summary.txt")
  if (file.exists(f_summary)) {
    shap_list$shap_summary <- data.table::fread(f_summary)
  }

  f_vg <- file.path(shap_dir, "variable_genes.txt")
  if (file.exists(f_vg)) {
    shap_list$variable_genes <- read.table(f_vg, header = TRUE, stringsAsFactors = FALSE)$variable_genes
  }

  f_kd <- file.path(shap_dir, "key_drivers.txt")
  if (file.exists(f_kd)) {
    shap_list$key_drivers <- read.table(f_kd, header = TRUE, stringsAsFactors = FALSE)$key_drivers
  }

  f_long <- file.path(shap_dir, "shap_long.rds")
  if (file.exists(f_long)) {
    shap_list$shap_long <- readRDS(f_long)
  }

  f_shapres <- file.path(shap_dir, "shap_result.rds")
  if (file.exists(f_shapres)) {
    shap_list$shap_result <- readRDS(f_shapres)
  }

  f_model <- file.path(shap_dir, "model.rds")
  if (file.exists(f_model)) {
    shap_list$model <- readRDS(f_model)
  }

  f_auc <- file.path(shap_dir, "auc.txt")
  if (file.exists(f_auc)) {
    shap_list$test_auc <- suppressWarnings(as.numeric(readLines(f_auc)[1]))
  }

  f_aucfold <- file.path(shap_dir, "auc_per_fold.txt")
  if (file.exists(f_aucfold)) {
    auc_df <- data.table::fread(f_aucfold)
    shap_list$auc_per_fold <- auc_df$auc_per_fold
    shap_list$mean_auc <- mean(shap_list$auc_per_fold, na.rm = TRUE)
  }

  # -------------------------
  # Attempt to infer mode from directory name or info.txt
  # Try to extract mode from folder name
  # -------------------------
  mode_guess <- sub("^shap_output_", "", basename(shap_dir))
  if (!mode_guess %in% c("full", "split", "cv")) mode_guess <- NA_character_

  info_file <- file.path(shap_dir, "info.txt")
  if (file.exists(info_file)) {
    info_lines <- readLines(info_file)

    # New header line looks like: "SHAP output summary (mode = 'cv'):"
    header <- grep("^SHAP output summary \\(mode = '", info_lines, value = TRUE)
    if (length(header) > 0) {
      mode_from_info <- sub(".*\\(mode = '([^']+)'\\).*", "\\1", header[1])
      if (mode_from_info %in% c("full", "split", "cv")) mode_guess <- mode_from_info
    }

    # Pull conditions/comparison if you want (optional, harmless if missing)
    cond_line <- grep("^conditions:", info_lines, value = TRUE)
    if (length(cond_line) > 0) shap_list$conditions <- sub("^conditions:\\s*", "", cond_line[1])

    cmp_line <- grep("^comparison:", info_lines, value = TRUE)
    if (length(cmp_line) > 0) shap_list$comparison_string <- sub("^comparison:\\s*", "", cmp_line[1])

    # Pull package versions if present
    xgb_line <- grep("^xgboost:", info_lines, value = TRUE)
    shap_line <- grep("^SHAPforxgboost:", info_lines, value = TRUE)
    if (length(xgb_line) > 0 || length(shap_line) > 0) {
      shap_list$package_versions <- list(
        xgboost = if (length(xgb_line) > 0) sub("^xgboost:\\s*", "", xgb_line[1]) else NA_character_,
        SHAPforxgboost = if (length(shap_line) > 0) sub("^SHAPforxgboost:\\s*", "", shap_line[1]) else NA_character_
      )
    }
  }

  # mode_guess <- gsub("^shap_output_?", "", basename(shap_dir))
  #
  # # If not a valid mode, fallback to reading from info.txt
  # if (!mode_guess %in% c("full", "split", "cv")) {
  #   info_file <- file.path(shap_dir, "info.txt")
  #   if (file.exists(info_file)) {
  #     info_lines <- readLines(info_file)
  #     mode_line <- grep("mode = '", info_lines, value = TRUE)
  #     if (length(mode_line) > 0) {
  #       mode_guess <- sub(".*mode = '([^']+)'.*", "\\1", mode_line[1])
  #     } else {
  #       warning("Mode could not be inferred from info.txt.")
  #     }
  #   } else {
  #     warning("info.txt not found. Mode could not be inferred.")
  #   }
  # }

  shap_list$mode <- mode_guess

  seurat_obj@misc$shap <- shap_list
  message("SHAP outputs successfully loaded into seurat_obj@misc$shap.")
  return(seurat_obj)
}
