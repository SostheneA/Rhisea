#' @title Run HISEA Mixed Stock Analysis
#'
#' @description
#' This is the main wrapper function and core set of utilities for running the HISEA mixed-stock analysis
#' framework, allowing simulation, analysis, or bootstrap estimation of stock composition from mixture samples.
#'
#' Supported operation modes:
#' - **SIMULATION**: Simulate mixtures based on known proportions and evaluate performance of classification and estimators.
#' - **ANALYSIS**: Apply trained classifier to real mixture data to estimate stock proportions.
#' - **BOOTSTRAP**: Resample real mixture to evaluate variability of estimates.
#'
#' Supported classifiers: LDA, QDA, Random Forest, SVM, k-NN, ANN, XGBoost, Naive Bayes, Mclust, MLR.
#' Supported estimators: RAW, Cook, Constrained Cook, EM (Millar), Maximum Likelihood.
#'
#' Includes integrated 10-fold cross-validation and model quality evaluation (accuracy, kappa, F1, etc.).
#'
#' @param type Character. "SIMULATION", "ANALYSIS" or "BOOTSTRAP".
#' @param np Integer. Number of populations (stocks).
#' @param nv Integer. Number of variables.
#' @param seed_val Integer. Random seed for reproducibility.
#' @param nsamps Integer. Number of replicates.
#' @param Nmix Integer. Sample size of the simulated mixture (for SIMULATION only).
#' @param actual Numeric vector. True proportions used in simulation.
#' @param baseline_path Character. File path to the baseline `.std` file.
#' @param mix_path Character. File path to the mixture `.mix` file.
#' @param export_csv Logical. Whether to export summary and confusion matrix to CSV.
#' @param output_dir Character. Output directory.
#' @param verbose Logical. Print progress messages.
#' @param method_class Character. Classification method (e.g., "LDA", "RF", "SVM", etc.).
#' @param stocks_names Character vector. Optional vector of stock names.
#' @param resample_baseline Logical. Resample the baseline for each replicate.
#' @param resampled_baseline_sizes Integer vector. Sizes of resamples per stock.
#' @param phi_method Character. "standard" or "cv" (cross-validation-based confusion matrix).
#' @param mclust_model_names Character vector. Models to test with Mclust.
#' @param mclust_perform_cv Logical. Whether to cross-validate Mclust.
#'
#' @return A list with:
#' \describe{
#'   \item{estimation_summary}{Summary table with mean, SD, and RMSE of estimates.}
#'   \item{classification_model}{Final trained classifier object.}
#'   \item{baseline_classification_quality}{Accuracy, Kappa, and per-class metrics.}
#'   \item{phi_matrix}{Estimated confusion matrix used in corrections.}
#'   \item{mixture_classification_details}{List with predicted pseudo-classes and likelihoods.}
#' }
#'
#' A `.rda` file of results is also saved in `output_dir`.
#'
#' @seealso \code{compute_cook_estimators}, \code{estimate_millar}, \code{estimate_ml},
#' \code{get_cv_metrics_and_phi}, \code{train_model}, \code{predict_model}, \code{.resample_baseline_data_helper}
#'
#' @examples
#' \dontrun{
#' run_hisea_all(type="SIMULATION",
#'              np=3, nv=5,
#'              actual=c(0.2,0.3,0.5),
#'              Nmix=200,
#'              baseline_path="baseline.std",
#'              method_class="RF",
#'              resample_baseline=TRUE,
#'              resampled_baseline_sizes=c(100,100,100),
#'              output_dir="results")
#' }
#'
#' @importFrom stats as.formula predict
#' @importFrom utils write.csv
#' @importFrom data.table data.table .SD
#' @importFrom MASS ginv
#' @importFrom utils write.csv
#' @export




# -------------------------------------------------------------------
# Helper: resample baseline data per stock
# -------------------------------------------------------------------
.resample_baseline_data_helper <- function(original_baseline_list,
                                          resampled_sizes,
                                          stock_names_for_error,
                                          nv_fallback) {
  if (length(original_baseline_list) != length(resampled_sizes)) {
    stop("Length of original_baseline_list must match resampled_sizes.")
  }
  new_list <- vector("list", length(original_baseline_list))
  names(new_list) <- names(original_baseline_list)
  for (j in seq_along(original_baseline_list)) {
    dat_j <- original_baseline_list[[j]]
    n_j   <- resampled_sizes[j]
    stock  <- if (!is.null(stock_names_for_error) &&
                  length(stock_names_for_error)==length(original_baseline_list)) {
      stock_names_for_error[j]
    } else paste0("Stock_", j)
    # zero requested empty matrix
    if (n_j == 0) {
      nc <- if (!is.null(dat_j) && !is.null(ncol(dat_j))) ncol(dat_j) else NA_integer_
      if (is.na(nc) && j>1) {
        for (k in (j-1):1) {
          if (!is.null(original_baseline_list[[k]]) &&
              !is.null(ncol(original_baseline_list[[k]]))) {
            nc <- ncol(original_baseline_list[[k]])
            break
          }
        }
      }
      if (is.na(nc)) nc <- nv_fallback
      new_list[[j]] <- matrix(numeric(0), nrow=0, ncol=nc)
      next
    }
    if (is.null(dat_j) || nrow(dat_j)==0) {
      stop(sprintf("Cannot sample %d from '%s' (baseline empty).", n_j, stock))
    }
    idx <- sample(seq_len(nrow(dat_j)), size=n_j, replace=TRUE)
    new_list[[j]] <- dat_j[idx, , drop=FALSE]
  }
  new_list
}


split_by_class_list <- function(data, labels) {
  fac <- as.factor(labels)
  lev <- levels(fac)
  out <- vector("list", length(lev))
  names(out) <- lev
  for (i in seq_along(lev)) {
    out[[i]] <- as.matrix(data[fac == lev[i], , drop = FALSE])
  }
  out
}

# -------------------------------------------------------------------
# Helper: train one classifier
# -------------------------------------------------------------------
`%||%` <- function(x, y) if (is.null(x)) y else x
train_model <- function(train_data, train_labels,
                        method_class,
                        np, nv,
                        generic_colnames,...) {
  m <- toupper(method_class)
  out <- list()
  args <- list(...)
  if (m=="LDA") {
    out$model <- compute_ldf_coefficients(split_by_class_list(train_data, train_labels))
  } else if (m=="LDA_MASS") {
    out$model <- MASS::lda(train_data, grouping=train_labels, prior=rep(1/np,np))
  } else if (m=="MCLUST") {
    out$model <- mclust::MclustDA(data=train_data, class=train_labels, verbose=FALSE)
  } else if (m=="QDA") {
    out$model <- MASS::qda(train_data, grouping=train_labels, prior=rep(1/np,np))

  } else if (m=="MLR") {
    # Prepare data with consistent column names
    feature_names <- paste0("Feature_", seq_len(ncol(train_data)))
    train_df <- data.frame(as.matrix(train_data))
    colnames(train_df) <- feature_names

    # Create complete dataframe with class
    df <- data.frame(
      Class = factor(train_labels),
      train_df
    )

    # Create formula
    f <- as.formula(paste("Class ~", paste(feature_names, collapse = " + ")))

    # Train the model
    out$model <- nnet::multinom(
      f,
      data = df,
      trace = FALSE,
      MaxNWts = max(2000, (ncol(train_data)+1)*(np-1)*2+500)
    )
    out$feature_names <- feature_names
    out$levels <- levels(factor(train_labels))
  } else if (m=="KNN") {
    k <- args$k %||% 5
    out$model <- list(train_data=train_data, train_labels=train_labels, k=k)
  } else if (m=="SVM") {
    kernel <- args$kernel %||% "radial"
    out$model <- e1071::svm(x=train_data, y=train_labels, kernel=kernel, probability=TRUE)
  } else if (m=="RF") {
    ntree <- args$ntree %||% 500
    # Ensure data is in the correct format
    train_matrix <- as.data.frame(train_data)
    colnames(train_matrix) <- paste0("Feature_", seq_len(ncol(train_data)))
    out$model <- randomForest::randomForest(x=train_matrix,
                                            y=factor(train_labels),
                                            ntree=ntree)
    out$feature_names <- colnames(train_matrix)
  } else if (m=="XGB") {
    nrounds <- args$nrounds %||% 100
    feature_names <- paste0("Feature_", seq_len(ncol(train_data)))
    train_matrix <- as.matrix(train_data)
    colnames(train_matrix) <- feature_names
    labels_num <- as.integer(as.factor(train_labels)) - 1
    dmat <- xgboost::xgb.DMatrix(data=train_matrix, label=labels_num)
    params <- modifyList(list(
      objective = "multi:softprob",
      num_class = np,
      eval_metric = "mlogloss",
      eta = 0.3,
      max_depth = 6
    ), args)  # surcharge avec ce qui est dans ...
    out$model <- xgboost::xgb.train(
      params = params,
      data = dmat,
      nrounds = nrounds,
      verbose = 0
    )
    out$feature_names <- feature_names
    out$np <- np
  } else if (m=="ANN") {
    size <- args$size %||% 5
    df <- data.frame(Class=train_labels, train_data)
    if (nv>0) colnames(df)[-1] <- generic_colnames
    f <- if (nv>0)
      as.formula(paste("Class ~", paste(generic_colnames, collapse=" + ")))
    else
      as.formula("Class ~ 1")
    out$model <- nnet::nnet(f, data=df,
                            size=size, trace=FALSE,
                            MaxNWts=max(5000,((nv+1)*size+(size+1)*np)*3+2000),
                            maxit=200)
  } else if (m=="NB") {
    out$model <- e1071::naiveBayes(x=train_data, y=train_labels)
  } else {
    stop("Unsupported method: ", method_class)
  }
  out
}


# -------------------------------------------------------------------
# Helper: predict from a trained model
# -------------------------------------------------------------------
predict_model <- function(model_info, newdata, method_class,
                          stocks_names, generic_colnames) {
  m <- toupper(method_class)
  np <- length(stocks_names)
  n  <- nrow(newdata)
  pc <- integer(n); prob <- matrix(1/np, nrow=n, ncol=np,
                                   dimnames=list(NULL,stocks_names))

  if (m=="LDA") {
    r <- classify_samples(newdata, model_info$model)
    pc   <- r$class
    prob <- r$likelihood

  } else if (m %in% c("LDA_MASS","QDA")) {
    pr   <- predict(model_info$model, newdata)
    pc   <- as.integer(pr$class)
    prob <- pr$posterior
    colnames(prob) <- stocks_names

  } else if (m=="MCLUST") {
    pr   <- predict(model_info$model, newdata=newdata)
    prob <- pr$z; colnames(prob) <- stocks_names
    pc   <- apply(prob,1,which.max)

  } else if (m=="MLR") {
    # Prepare test data
    newdata_df <- data.frame(as.matrix(newdata))
    colnames(newdata_df) <- model_info$feature_names

    tryCatch({
      # Predict classes
      pc_f <- predict(model_info$model, newdata = newdata_df, type = "class")
      pc <- as.integer(factor(pc_f, levels = stocks_names))

      # Predict probabilities
      pmat <- predict(model_info$model, newdata = newdata_df, type = "probs")

      # Ensure pmat is a matrix
      if (is.null(dim(pmat))) {
        pmat <- matrix(pmat, nrow = nrow(newdata), ncol = np)
      }

      # Check and normalize probabilities
      if (is.matrix(pmat)) {
        if (ncol(pmat) == np) {
          # Normalize if needed
          pmat <- t(apply(pmat, 1, function(x) x/sum(x)))
          colnames(pmat) <- stocks_names
          prob <- pmat
        } else {
          warning("MLR: Unexpected number of columns in probability matrix")
        }
      } else {
        warning("MLR: Could not get probability matrix")
      }
    }, error = function(e) {
      warning("MLR prediction error: ", e$message)
      # Return uniform probabilities in case of error
      prob <- matrix(1/np, nrow=nrow(newdata), ncol=np,
                     dimnames=list(NULL, stocks_names))
      pc <- rep(1, nrow(newdata))
    })

  } else if (m=="KNN") {
    k   <- model_info$model$k
    td  <- model_info$model$train_data
    tl  <- model_info$model$train_labels
    if (nrow(td)<=k) k <- max(1,nrow(td)-1)
    pc_f <- class::knn(train=td, test=newdata, cl=tl, k=k)
    pc   <- as.integer(factor(pc_f, levels=stocks_names))
    if (requireNamespace("kknn", quietly=TRUE)) {
      dftr <- data.frame(Class=tl, td)
      colnames(dftr)[-1] <- generic_colnames
      dfte <- data.frame(newdata)
      colnames(dfte) <- generic_colnames
      km   <- kknn::train.kknn(Class~., data=dftr,
                               kmax=k, kernel="rectangular", scale=FALSE)
      pp   <- predict(km, newdata=dfte, type="prob")
      if (is.matrix(pp) && ncol(pp)==np) {
        colnames(pp) <- stocks_names
        prob <- pp
      }
    }

  } else if (m=="SVM") {
    pr   <- predict(model_info$model, newdata, probability=TRUE)
    pc_f <- pr; pc <- as.integer(factor(pc_f, levels=stocks_names))
    pa   <- attr(pr,"probabilities")
    if (!is.null(pa) && is.matrix(pa)) {
      tmp <- matrix(0, nrow=n, ncol=np, dimnames=list(NULL,stocks_names))
      cn  <- colnames(pa)
      for (nm in intersect(cn,stocks_names)) tmp[,nm] <- pa[,nm]
      tmp <- tmp/rowSums(tmp)
      prob <- tmp
    }

  } else if (m=="RF") {
    newdata_df <- as.data.frame(newdata)
    colnames(newdata_df) <- model_info$feature_names

    # Predict classes
    pc_f <- predict(model_info$model, newdata_df, type="response")
    pc <- as.integer(factor(pc_f, levels=stocks_names))

    # Predict probabilities
    pmat <- predict(model_info$model, newdata_df, type="prob")
    if(is.matrix(pmat) && ncol(pmat) == np) {
      colnames(pmat) <- stocks_names
      prob <- pmat
    }

  } else if (m=="XGB") {
    # Prepare data with same feature names
    newdata_matrix <- as.matrix(newdata)
    colnames(newdata_matrix) <- model_info$feature_names

    # Create DMatrix for prediction
    dmat <- xgboost::xgb.DMatrix(newdata_matrix)

    # Predict probabilities
    pred_prob <- predict(model_info$model, dmat)

    # Reshape probabilities into matrix
    prob <- matrix(pred_prob, nrow=nrow(newdata), ncol=np, byrow=TRUE)
    colnames(prob) <- stocks_names

    # Determine predicted classes
    pc <- max.col(prob)

    # Normalize probabilities if needed
    prob <- t(apply(prob, 1, function(x) {
      x / sum(x)
    }))

  } else if (m=="ANN") {
    df    <- data.frame(newdata)
    if (ncol(newdata)>0) colnames(df) <- generic_colnames
    pc_f  <- predict(model_info$model, df, type="class")
    pc    <- as.integer(factor(pc_f, levels=stocks_names))
    pr    <- predict(model_info$model, df, type="raw")
    if (is.matrix(pr) && ncol(pr)==np) {
      colnames(pr) <- stocks_names; prob <- pr
    } else if (np==2 && is.vector(pr)) {
      prob[, stocks_names[1]] <- 1-pr; prob[, stocks_names[2]] <- pr
    }

  } else if (m=="NB") {
    pc_f <- predict(model_info$model, newdata=newdata, type="class")
    pc   <- as.integer(factor(pc_f, levels=stocks_names))
    pr   <- predict(model_info$model, newdata=newdata, type="raw")
    if (is.matrix(pr) && ncol(pr)==np) {
      colnames(pr) <- stocks_names; prob <- pr
    }
  }

  list(class=pc, likelihood=prob)
}

# -------------------------------------------------------------------
# Helper: 10-fold CV metrics + Phi
# -------------------------------------------------------------------
get_cv_metrics_and_phi <- function(data, labels, method_class,
                                   np, nv, stocks_names,
                                   generic_colnames,
                                   folds = 10,...) {
  if (!requireNamespace("caret", quietly=TRUE))
    stop("Package 'caret' is required for cross-validation.")

  set.seed(123)
  fold_idx <- caret::createFolds(labels, k=folds, list=TRUE)

  # Collect all predictions
  preds <- integer(length(labels))
  for (f in seq_along(fold_idx)) {
    test_i  <- fold_idx[[f]]
    train_i <- setdiff(seq_along(labels), test_i)

    tr_data <- data[train_i,,drop=FALSE]
    tr_lab  <- labels[train_i]
    ts_data <- data[test_i,,drop=FALSE]

    mi  <- train_model(tr_data, tr_lab,
                       method_class,
                       np, nv,
                       generic_colnames,...)
    pr  <- predict_model(mi, ts_data, method_class,
                         stocks_names, generic_colnames)
    preds[test_i] <- pr$class
  }

  # Build confusion
  pred_fac <- factor(preds, levels=seq_len(np), labels=stocks_names)
  cm       <- table(Predicted=pred_fac, Actual=labels)

  # Use caret for detailed metrics
  cmr <- suppressWarnings(
    caret::confusionMatrix(data=pred_fac,
                           reference=labels,
                           mode="everything")
  )

  overall_acc   <- as.numeric(cmr$overall["Accuracy"])
  overall_kappa <- as.numeric(cmr$overall["Kappa"])

  byC <- cmr$byClass

  # Handle multiclass (matrix) vs binary (vector)
  if (is.matrix(byC)) {
    sens    <- as.numeric(byC[,"Sensitivity"])
    spec    <- as.numeric(byC[,"Specificity"])
    prec    <- as.numeric(byC[,"Precision"])
    recall  <- as.numeric(byC[,"Recall"])
    f1      <- as.numeric(byC[,"F1"])
    balAcc  <- as.numeric(byC[,"Balanced Accuracy"])
  } else {
    # Binary: byC is a named vector of length n_metrics
    sens    <- rep(as.numeric(byC["Sensitivity"]), np)
    spec    <- rep(as.numeric(byC["Specificity"]), np)
    prec    <- rep(as.numeric(byC["Precision"]), np)
    recall  <- rep(as.numeric(byC["Recall"]), np)
    f1      <- rep(as.numeric(byC["F1"]), np)
    balAcc  <- rep(as.numeric(byC["Balanced Accuracy"]), np)
  }

  # Build Phi = P(assigned=i | true=j)
  phi_m <- as.matrix(prop.table(cm, margin = 2))

  return(list(
    confusion_matrix            = cm,
    accuracy                    = overall_acc,
    kappa                       = overall_kappa,
    sensitivity_by_class        = sens,
    specificity_by_class        = spec,
    precision_by_class          = prec,
    recall_by_class             = recall,
    f1_by_class                 = f1,
    balanced_accuracy_by_class  = balAcc,
    phi_matrix                  = phi_m
  ))
}

#' Process HISEA Input Data
#'
#' @param input Input data (file path or data frame)
#' @param type Type of input ("std" or "mix")
#' @param nv Number of variables
#' @param stock_col Stock column name
#' @param var_cols_mix Variables for mixture data
#' @param var_cols_std Variables for standard data
#' @param out_path Output path
#' @param stocks_names Stock names
#' @return Processed data
#' @keywords internal
process_hisea_input <- function(input, type = c("std", "mix"),
                                nv = NULL,
                                stock_col = NULL,
                                var_cols_mix = NULL,
                                var_cols_std = NULL) {
  type <- match.arg(type)

  # Case 1 read from existing file
  if (is.character(input) && file.exists(input)) {
    if (is.null(nv)) stop("nv is required to read a .std or .mix file")
    if (type == "std") {
      return(read_baseline(filepath = input, nv = nv))
    } else {
      return(read_mixture(filepath = input, nv = nv))
    }
  }

  # Case 2 input is a data.frame / data.table
  if (inherits(input, "data.frame") || inherits(input, "data.table")) {
    input <- data.table::data.table(input)

    if (type == "std") {
      # Check stock column
      if (is.null(stock_col)) stop("stock_col is required for .std")
      if (!(stock_col %in% names(input)))
        stop("The stock column '", stock_col, "' does not exist in the provided data.frame.")

      # Check variable columns
      if (!all(var_cols_std %in% names(input)))
        stop("At least one of the specified variable columns does not exist in the provided data.frame.")

      # Get stock levels in the original order of appearance
      stock_order <- unique(input[[stock_col]])

      # Create list of matrices in original order
      std_list <- lapply(stock_order, function(s) {
        as.matrix(input[input[[stock_col]] == s, ..var_cols_std])
      })

      # Remove names to return [[1]], [[2]], etc.
      names(std_list) <- NULL

      # Reorder to match stocks_names
      std_list <- std_list[match(stocks_names, stock_order)]

      return(std_list)

    } else { # type == "mix"
      # Check variable columns
      if (!all(var_cols_mix %in% names(input)))
        stop("At least one of the specified variable columns does not exist in the provided data.frame.")

      # Convert directly to matrix
      mix_matrix <- as.matrix(
        do.call(cbind, lapply(input[, ..var_cols_mix], as.numeric))
      )
      colnames(mix_matrix) <- c("V1", "V2")
      return(mix_matrix)
    }
  }

  stop("Input must be either an existing file path or a data.frame/data.table.")
}

#' Write Standard Data to File
#'
#' @param df Data frame
#' @param var_cols_std Variables to write
#' @param stock_col Stock column
#' @param file_path Output file path
#' @return None
#' @keywords internal
#' @export
write_std_from_dataframe <- function(df, stock_col, var_cols_std, file_path = "hisea.std") {
  # Check that the stock column exists
  if (!(stock_col %in% names(df))) {
    stop("The specified stock column does not exist in the data frame.")
  }

  # Check that all selected variable columns exist
  if (!all(var_cols_std %in% names(df))) {
    stop("One or more specified variable columns are missing from the data frame.")
  }
  df=data.table(df)
  # Split data by stock
  stock_groups <- split(df[, ..var_cols_std], df[[stock_col]])

  # Write data in .std format
  con <- file(file_path, "w")
  for (i in seq_along(stock_groups)) {
    mat <- as.matrix(stock_groups[[i]])
    apply(mat, 1, function(row) writeLines(paste(row, collapse = " "), con))
    if (i < length(stock_groups)) writeLines("NEXT STOCK", con)
  }
  writeLines("End of baseline data", con)
  writeLines("End of file", con)
  close(con)
}

#' Write Mixture Data to File
#'
#' @param df Data frame
#' @param var_cols_mix Variables to write
#' @param file_path Output file path
#' @return None
#' @keywords internal
#' @export
write_mix_from_dataframe <- function(df, var_cols_mix, file_path = "hisea.mix") {
  # Check that all selected variable columns exist
  if (!all(var_cols_mix %in% names(df))) {
    stop("One or more specified variable columns are missing from the data frame.")
  }

  # Convert to data.table for flexibility
  df <- data.table::data.table(df)

  # Extract only the variables of interest
  data_to_write <- as.matrix(df[, ..var_cols_mix])

  # Write in HISEA .mix format
  con <- file(file_path, "w")
  apply(data_to_write, 1, function(row) writeLines(paste(row, collapse = " "), con))
  writeLines("End of mixed sample", con)
  writeLines("End of file", con)
  close(con)
}


#' @title Internal Main Function for HISEA Run
#' @description Does the full analysis/simulation/bootstrap run

# -------------------------------------------------------------------
# Main wrapper: run_hisea_all()
# -------------------------------------------------------------------
#' Run Complete HISEA Stock Composition Analysis
#'
#' This function performs comprehensive stock composition analysis using various
#' statistical methods including LDA, Random Forest, SVM, XGBoost, and others.
#'
#' @param type Character. Type of analysis to perform. Default: "ANALYSIS"
#' @param np Integer. Number of populations/stocks. Default: 2
#' @param nv Integer. Number of variables/loci. Default: 2
#' @param seed_val Integer. Random seed for reproducibility. Default: 123456
#' @param nsamps Integer. Number of bootstrap samples. Default: 1000
#' @param Nmix Integer. Size of mixture sample. Default: 100
#' @param actual Numeric vector. True composition proportions. Default: c(0.5, 0.5)
#' @param var_cols_std Character vector of column names for baseline variables
#' @param var_cols_mix Character vector of column names for mixture variables
#' @param stock_col Character name of stock column in baseline data
#' @param baseline_input Data frame or file path for baseline data
#' @param mix_input Data frame or file path for mixture data
#' @param export_csv Logical. Whether to export results to CSV. Default: TRUE
#' @param output_dir Character. Output directory path. Default: getwd()
#' @param verbose Logical. Whether to print progress messages. Default: TRUE
#' @param method_class Character. Classification method ("LDA", "RF", "SVM", etc.). Default: "LDA"
#' @param stocks_names Character vector. Names of stocks. Default: NULL
#' @param resample_baseline Logical. Whether to resample baseline data. Default: FALSE
#' @param resampled_baseline_sizes Integer vector. Sizes for resampled baseline. Default: NULL
#' @param phi_method Character. Method for phi calculation. Default: "default"
#' @param mclust_model_names Character vector. Model names for mclust. Default: NULL
#' @param mclust_perform_cv Logical. Whether to perform cross-validation for mclust. Default: FALSE
#' @param ... Additional arguments passed to the underlying classifier.
#'   Certain classifiers support hyperparameter tuning via these arguments:
#'
#'   - **Random Forest (RF)**: `ntree`, `mtry`, `nodesize`, `maxnodes`
#'     - Example: `run_hisea_all(..., method_class="RF", ntree=5000, mtry=3)`
#'
#'   - **Support Vector Machine (SVM)**: `cost`, `gamma`, `kernel`
#'     - Example: `run_hisea_all(..., method_class="SVM", cost=10, kernel="radial")`
#'
#'   - **XGBoost (XGB)**: `nrounds`, `max_depth`, `eta`, `subsample`, `colsample_bytree`
#'     - Example: `run_hisea_all(..., method_class="XGB", nrounds=100, max_depth=4)`
#'
#'   - **k-Nearest Neighbors (KNN)**: `k`
#'     - Example: `run_hisea_all(..., method_class="KNN", k=5)`
#'
#'   - **Artificial Neural Network (ANN)**: `size`, `decay`, `maxit`
#'     - Example: `run_hisea_all(..., method_class="ANN", size=10, decay=0.01)`
#'
#'   - **Naive Bayes (NB)**: `laplace`
#'     - Example: `run_hisea_all(..., method_class="NB", laplace=1)`
#'
#'   ⚠ Ensure that the names of arguments match those expected by the classifier.
#'     Unrecognized arguments may be ignored or raise an error.

#'
#' @return List containing:
#' \describe{
#'   \item{estimates}{Matrix of stock composition estimates}
#'   \item{mean}{Mean estimates across bootstrap samples}
#'   \item{sd}{Standard deviations of estimates}
#'   \item{performance}{Performance metrics if applicable}
#' }
#'
#' @export
#' @examples
#' \dontrun{
#' # Basic analysis with LDA
#' result <- run_hisea_all(
#'   np = 2,
#'   nv = 3,
#'   method_class = "LDA",
#'   baseline_path = "my_baseline.txt",
#'   mix_path = "my_mixture.txt"
#' )
#' }
run_hisea_all <- function(type = "ANALYSIS",
                          np, nv,
                          seed_val = 123456,
                          var_cols_std=NULL,
                          var_cols_mix=NULL,
                          stock_col=NULL,
                          nsamps   = 1000,
                          Nmix     = 100,
                          actual   = NULL,
                          baseline_input = NULL,
                          mix_input      = NULL,
                          export_csv    = FALSE,
                          output_dir    = ".",
                          verbose       = FALSE,
                          method_class  = "LDA",
                          stocks_names  = NULL,
                          resample_baseline    = FALSE,
                          resampled_baseline_sizes = NULL,
                          phi_method  = c("standard","cv"),
                          mclust_model_names = NULL,
                          mclust_perform_cv  = TRUE,...) {

  type         <- toupper(type)
  method_class <- toupper(method_class)
  phi_method   <- match.arg(phi_method)
  if (is.null(var_cols_std)) var_cols_std <- paste0("V", 1:nv)

  if (type=="ANALYSIS" && resample_baseline) {
    warning("resample_baseline only for SIMULATION/BOOTSTRAP; disabling.")
    resample_baseline <- FALSE
  }
  if (resample_baseline) {
    if (!type %in% c("SIMULATION","BOOTSTRAP"))
      stop("resample_baseline=TRUE requires type='SIMULATION' or 'BOOTSTRAP'.")
    if (is.null(resampled_baseline_sizes) ||
        length(resampled_baseline_sizes)!=np)
      stop("resampled_baseline_sizes must be length np.")
  }

  # read baseline
  base_list <- process_hisea_input(baseline_input, type="std", nv=nv,
                                   stock_col=stock_col, var_cols_std=var_cols_std)
  if (length(base_list)!=np)
    stop(sprintf("Baseline has %d stocks but np=%d.", length(base_list), np))
  base_data <- do.call(rbind, base_list)
  base_num  <- rep(seq_along(base_list), times=sapply(base_list,nrow))

  if (!is.null(stocks_names) && length(stocks_names)==np)
    stock_names_internal <- stocks_names
  else
    stock_names_internal <- paste0("Stock_", seq_len(np))

  base_fac <- factor(base_num, levels=1:np, labels=stock_names_internal)

  # allocate results
  est_names <- c("RAW","COOK","COOKC","EM","ML")
  all_res   <- array(NA_real_, dim=c(nsamps,np,length(est_names)),
                     dimnames=list(NULL, stock_names_internal, est_names))

  generic_colnames <- make.names(paste0("V",1:nv), unique=TRUE)
  set.seed(seed_val)
  full_mix_boot   <- NULL
  final_model     <- NULL
  final_metrics   <- NULL
  final_Phi       <- NULL
  final_mix_pc    <- NULL
  final_mix_like  <- NULL

  if (verbose) message("Running ",nsamps," x ",type," with ",method_class,
                       " (phi=",phi_method,")")

  for (i in seq_len(nsamps)) {
    if (verbose && (i==1||i==nsamps||i%%max(1,floor(nsamps/10))==0))
      message(" Iter ",i,"/",nsamps)

    # baseline for this iter
    if (resample_baseline) {
      bl <- .resample_baseline_data_helper(base_list,
                                          resampled_baseline_sizes,
                                          stock_names_internal, nv)
      bdata <- do.call(rbind,bl)
      bnum  <- rep(seq_along(bl), times=sapply(bl,nrow))
      bfac  <- factor(bnum, levels=1:np, labels=stock_names_internal)
      if (nrow(bdata)==0 && sum(resampled_baseline_sizes)>0) {
        warning("Iter ",i,": empty resampled baseline; skipping.")
        all_res[i,,] <- NA; next
      }
    } else {
      bdata <- base_data; bfac <- base_fac
    }

    # compute metrics & Phi
    if (phi_method=="cv") {
      cmv    <- get_cv_metrics_and_phi(bdata, bfac,
                                       method_class,
                                       np,nv,
                                       stock_names_internal,
                                       generic_colnames, folds=10,...)
      quality <- cmv; iter_Phi <- as.matrix(cmv$phi_matrix)
      # train final model on full baseline
      mi <- train_model(bdata, bfac,
                        method_class, np, nv,
                        generic_colnames,...)
    } else {
      mi <- train_model(bdata, bfac,
                        method_class, np, nv,
                        generic_colnames,...)
      prb <- predict_model(mi, bdata, method_class,
                           stock_names_internal, generic_colnames)
      cm  <- table(Predicted=factor(prb$class,
                                    levels=1:np,
                                    labels=stock_names_internal),
                   Actual   =bfac)
      acc <- sum(diag(cm))/sum(cm)
      tot <- sum(cm)
      pe  <- sum(rowSums(cm)*colSums(cm))/(tot^2)
      kap <- if (abs(1-pe)>.Machine$double.eps) (acc-pe)/(1-pe) else NA_real_
      quality    <- list(confusion_matrix=cm,
                         accuracy=acc, kappa=kap)
      iter_Phi   <- as.matrix(prop.table(cm, margin=2))
    }

    # invert or pseudo-invert Phi
    iter_Phi <- as.matrix(iter_Phi)
    iter_Phi <- matrix(as.numeric(iter_Phi),
                                 nrow = nrow(iter_Phi),
                                 ncol = ncol(iter_Phi))

    phi_det <- tryCatch({
      det(iter_Phi)
    }, error = function(e) {
      warning("Could not compute det(iter_Phi): ", e$message)
      0
    })

    if (abs(phi_det) < 1e-12) {
      inv_Phi <- MASS::ginv(iter_Phi)
    } else {
      inv_Phi <- solve(iter_Phi)
    }

    # store final run artifacts
    if (i==1 || resample_baseline) {
      final_model   <- mi$model
      final_metrics <- quality
      final_Phi     <- iter_Phi
    }

    # generate mixture
    if (type=="SIMULATION") {
      if (is.null(actual)) stop("'actual' needed for SIMULATION.")
      mix_s <- simulate_mixture(base_list, actual, Nmix)
    } else if (type=="ANALYSIS") {
      mix_s <- process_hisea_input(mix_input, type="mix", nv=nv, var_cols_mix=var_cols_mix)
    } else if (type=="BOOTSTRAP") {
      if (i==1) {
        full_mix_boot <-process_hisea_input(mix_input, type="mix", nv=nv, var_cols_mix=var_cols_mix)
        if (nrow(full_mix_boot)==0)
          stop("Bootstrap: mixture empty.")
      }
      idx   <- sample(seq_len(nrow(full_mix_boot)), replace=TRUE)
      mix_s <- full_mix_boot[idx,,drop=FALSE]
    } else stop("Unknown type: ",type)

    if (nrow(mix_s)==0) {
      warning("Iter ",i,": empty mixture; skipping.")
      all_res[i,,] <- NA; next
    }

    # classify mixture
    pr_m <- predict_model(mi, mix_s, method_class,
                          stock_names_internal, generic_colnames)
    mix_pc   <- pr_m$class
    mix_like <- pr_m$likelihood
    mix_like[is.na(mix_like)] <- 1/np
    mix_like <- t(apply(mix_like,1,function(r) r/sum(r)))

    if (i==1 || resample_baseline) {
      final_mix_pc   <- mix_pc
      final_mix_like <- mix_like
    }

    # HISEA estimators
    raw   <- prop.table(tabulate(mix_pc,nbins=np))
    ck    <- compute_cook_estimators(mix_pc, inv_Phi, np)
    emv   <- estimate_millar(mix_pc, iter_Phi, np, verbose=FALSE, max_iter=200)
    mlp   <- estimate_ml(mix_like, np, verbose=FALSE, max_iter=200)

    all_res[i,,"RAW"]   <- raw
    all_res[i,,"COOK"]  <- ck$cook
    all_res[i,,"COOKC"] <- ck$cook_constrained
    all_res[i,,"EM"]    <- emv
    all_res[i,,"ML"]    <- mlp
  }

  # summary & optional CSV
  actual_sum <- if (type=="SIMULATION") actual else NULL
  summary_l <- create_hisea_summary_report(all_res,
                                           actual_sum,
                                           run_type=type)
  if (export_csv) {
    if (!dir.exists(output_dir)) dir.create(output_dir, recursive=TRUE)
    suf <- paste0("_",method_class,"_",type)
    write.csv(as.data.frame(summary_l$mean_estimates),
              file.path(output_dir,paste0("mean",suf,".csv")),
              row.names=FALSE)
    write.csv(as.data.frame(summary_l$sd_estimates),
              file.path(output_dir,paste0("sd",suf,".csv")),
              row.names=FALSE)
    if (!is.null(summary_l$mse_estimates)) {
      write.csv(as.data.frame(summary_l$mse_estimates),
                file.path(output_dir,paste0("rmse",suf,".csv")),
                row.names=FALSE)
    }
    cmf <- final_metrics$confusion_matrix
    if (!is.null(cmf)) {
      write.csv(as.data.frame.matrix(cmf),
                file.path(output_dir,paste0("confMat",suf,".csv")),
                row.names=TRUE)
    }
  }

  timestamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
  # save full return
  out <- list(estimation_summary             = summary_l,
              classification_model           = final_model,
              baseline_classification_quality = final_metrics,
              phi_matrix                     = final_Phi,
              mixture_classification_details = list(
                pseudo_classes = final_mix_pc,
                likelihoods    = final_mix_like
              ))
  save(out, file=file.path(output_dir,
                           paste0("result_", method_class, "_", type, "_", timestamp,".rda")))
  if (verbose) message("run_hisea_all() done.")
  summary_l
}
