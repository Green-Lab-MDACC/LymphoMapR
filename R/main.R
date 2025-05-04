# prep expression data
data_prep_expr <- function(expr_pred, data_type) {
  if (data_type %in% c("RNAseq_TPM_log2p1", "microarray_lognorm")) {
    expr_train <- as.matrix(assay(train_data, "logTPM"))
    genes_common <- intersect(rownames(expr_train), rownames(expr_pred))
    expr_train <- expr_train[genes_common, ]
    expr_pred <- expr_pred[genes_common, ]
  } else if (data_type %in% c("RNAseq_counts", "NanoString_norm_counts")) {
    expr_train <- as.matrix(assay(train_data, "counts"))
    genes_common <- intersect(rownames(expr_train), rownames(expr_pred))
    expr_train <- expr_train[genes_common, ]
    expr_pred <- expr_pred[genes_common, ]
    
    # convert counts to vst
    expr_train <- counts2vst(expr_train)
    expr_pred <- counts2vst(expr_pred)
    
    genes_common <- intersect(rownames(expr_train), rownames(expr_pred))
    expr_train <- expr_train[genes_common, ]
    expr_pred <- expr_pred[genes_common, ]
  } else {
    stop("data_type is not defined in the function")
  }

  return(list(train = expr_train,
              pred = expr_pred))
}


# all DEGs between archetypes in train data
get_train_deg_archetype <- function(gene_pool) {
  genes_diff <- setdiff(gene_pool, rownames(train_data))
  if (length(genes_diff) > 0) {
    message("Removed ", length(genes_diff), " genes that are not in the training data")
  }
  gene_pool <- intersect(rownames(train_data), gene_pool)
  message(length(gene_pool), " genes are used in total")
  
  # deg
  message("Calculating DEGs between archetypes")
  tb_deg <- get_deg_archetype(gene_pool)
  return(tb_deg)
}


# data prep and feature selection
data_prep <- function(expr_pred,
                      data_type,
                      n_features_per_archetype = 50,
                      min_fract_poszscore_group = 0.5,
                      max_fract_poszscore_rest = 0.5,
                      min_log2FoldChange = 0.58,
                      padj_cutoff = 0.05,
                      verbose = FALSE) {
  # train data
  y_train <- train_data$LymphoMAP
  names(y_train) <- train_data$sample_id
  
  # expr
  expr <- data_prep_expr(expr_pred, data_type)
  expr_train <- expr$train
  expr_pred <- expr$pred
  
  mt <- cbind(expr_train, expr_pred)
  sample_datasets <- structure(c(rep("train", ncol(expr_train)), 
                                 rep("pred", ncol(expr_pred))),
                               names = colnames(mt))
  
  # align pred expr to train expr
  message("Aligning expression distribution to training data for all genes")
  if (verbose) {
    expr_pred <- combat_pred_to_train(expr_train, expr_pred)
  } else {
    invisible(capture.output(
      suppressMessages(
        expr_pred <- combat_pred_to_train(expr_train, expr_pred)
      )
    ))
  }
  
  # deg
  message("Calculating DEGs between archetypes")
  gene_pool <- rownames(expr_train)
  if (verbose) {
    tb_deg <- get_deg_archetype(gene_pool)
  } else {
    invisible(capture.output(
      suppressMessages(
        tb_deg <- get_deg_archetype(gene_pool)
      )
    ))
  }
  
  # feature selection
  message("Feature selection")
  tb <- tb_deg %>%
    filter(log2FoldChange > min_log2FoldChange & padj < padj_cutoff) %>%
    dplyr::filter(fract_poszscore_group > min_fract_poszscore_group & 
                    fract_poszscore_rest < max_fract_poszscore_rest)
  
  tb_feature <- get_top_unique_genes(n_features_per_archetype, tb)
  features <- tb_feature$gene_name
  feature_archetypes <- tb_feature %>%
    select(gene_name, LymphoMAP) %>%
    deframe() %>%
    factor(levels = unique(tb_feature$LymphoMAP))
  
  # data matrix reduced to selected features
  mt_train <- expr_train[features, ]
  mt_pred <- expr_pred[features, ]
  
  # combat on reduced features
  # optional, not harmful
  # required if batch effect still remains in the reduced feature space
  message("Aligning expression distribution to training data for selected features")
  if (verbose) {
    mt_pred <- combat_pred_to_train(mt_train, mt_pred)
  } else {
    invisible(capture.output(
      suppressMessages(
        mt_pred <- combat_pred_to_train(mt_train, mt_pred)
      )
    ))
  }
  
  # transpose to samples by features
  mt_train <- t(mt_train)
  mt_pred <- t(mt_pred)
  
  # scale expression across samples based on mean and sd of train data
  message("Scale expressions based on training data")
  scaled_data <- scale_expr(mt_train, mt_pred)
  mt_train <- scaled_data$train
  mt_pred <- scaled_data$pred

  return(list(X_train = mt_train,
              y_train = y_train,
              X_pred = mt_pred,
              feature = tb_feature,
              feature_archetypes = feature_archetypes))
}


#' run_LymphoMap
#'
#' @description
#' Run LymphoMAP classification
#'
#' @param expr_pred gene by sample expression matrix of prediction data. Row names: gene symbols, column names: sample names. 
#' @param data_type expression data type from "RNAseq_TPM_log2p1", "microarray_lognorm", "RNAseq_counts", "NanoString_norm_counts". Default = "RNAseq_TPM_log2p1".
#' @param LymphoMAP_prob_cutoff LymphoMAP assignment threshold of classification probability. Default = 0.
#' @param n_features_per_archetype number of features per archetype used in model training. Default = 50.
#' @param min_fract_poszscore_group criterion for selected features that have minimum fraction of samples within archetype. Default = 0.5.
#' @param max_fract_poszscore_rest criterion for selected features that have maximum fraction of samples outside archetype. Default = 0.5.
#' @param min_log2FoldChange criterion for selected features that have minimum log2FoldChange of archetype vs rest. Default = 0.58.
#' @param padj_cutoff criterion for selected features that have maximum padj of archetype vs rest. Default = 0.05.
#' @param verbose logical value if showing all messages. Default = FALSE.
#' 
#' @details
#' run_LymphoMap integrates prediction expression matrix with training data, selects archetype-defining features from common genes shared by both datasets, fits a naive Bayes model, and finally assign LymphoMAP labels based on classification probabilities for prediction samples. 
#'
#' @return A named list of class `"LymphoMapObject"` with the following components:
#' \describe{
#'   \item{data}{A named list containing:
#'     \describe{
#'       \item{\code{X_train}}{Training expression matrix (numeric).}
#'       \item{\code{y_train}}{Training labels (named factor).}
#'       \item{\code{X_pred}}{Prediction expression matrix (numeric).}
#'       \item{\code{feature}}{Selected features (tibble).}
#'       \item{\code{feature_archetypes}}{Archetypes of selected features (named factor).}
#'     }
#'   }
#'   \item{\code{nb_model}}{An object of class \code{"naive_bayes"} containing the trained model.}
#'   \item{\code{cv_confmat}}{A confusion matrix object returned by \code{caret::confusionMatrix()}.}
#'   \item{\code{pred}}{Final LymphoMAP label assignments and probabilities. The last three columns contain the classification probability matrix for all archetypes. (tibble).}
#'   \item{\code{y_pred}}{Predicted labels (named factor).}
#'   \item{\code{y_prob}}{Probabilities of predicted labels (named numeric vector).}
#'   \item{\code{parameters}}{A list of parameters used in the function call.}
#' }
#' @export
run_LymphoMap <- function(expr_pred,
                          data_type = c("RNAseq_TPM_log2p1", "microarray_lognorm", "RNAseq_counts", "NanoString_norm_counts"),
                          LymphoMAP_prob_cutoff = 0,
                          n_features_per_archetype = 50,
                          min_fract_poszscore_group = 0.5,
                          max_fract_poszscore_rest = 0.5,
                          min_log2FoldChange = 0.58,
                          padj_cutoff = 0.05,
                          verbose = FALSE) {
  data_type <- match.arg(data_type)
  
  # data prep 
  message("Data prep")
  data <- data_prep(expr_pred,
                    data_type,
                    n_features_per_archetype = n_features_per_archetype,
                    min_fract_poszscore_group = min_fract_poszscore_group,
                    max_fract_poszscore_rest = max_fract_poszscore_rest,
                    min_log2FoldChange = min_log2FoldChange,
                    padj_cutoff = padj_cutoff,
                    verbose = verbose)
  mt_train <- data$X_train
  mt_pred <- data$X_pred
  y_train <- data$y_train
  feature_archetypes <- data$feature_archetypes
  
  # naive bayes
  message("Naive Bayes model training")
  nb_model <- nb_clf(mt_train, y_train)
  
  # cv performance
  message("Performance evaluation via cross validation")
  if (verbose) {
    cv_confmat <- nb_clf_cv(mt_train, y_train)
  } else {
    invisible(capture.output(
      suppressMessages(
        cv_confmat <- nb_clf_cv(mt_train, y_train)
      )
    ))
  }
  
  # pred
  message("Final prediction")
  tb_pred <- nb_clf_predict(nb_model, mt_pred)
  
  # pack res
  parameters <- list(LymphoMAP_prob_cutoff = LymphoMAP_prob_cutoff,
                     n_features_per_archetype = n_features_per_archetype,
                     min_fract_poszscore_group = min_fract_poszscore_group,
                     max_fract_poszscore_rest = max_fract_poszscore_rest,
                     min_log2FoldChange = min_log2FoldChange,
                     padj_cutoff = padj_cutoff)
  
  classes <- c(levels(y_train), "Unassigned")
  tb_pred <- tb_pred %>%
    mutate(LymphoMAP = as.character(LymphoMAP)) %>%
    mutate(LymphoMAP = ifelse(LymphoMAP_prob > LymphoMAP_prob_cutoff, LymphoMAP, "Unassigned")) %>%
    mutate(LymphoMAP = factor(LymphoMAP, levels = classes))

  y_pred <- tb_pred %>%
    select(sample_id, LymphoMAP) %>%
    deframe()
  y_prob <- tb_pred %>%
    select(sample_id, LymphoMAP_prob) %>%
    deframe()
  
  res <- list(nb_model = nb_model,
              cv_confmat = cv_confmat,
              data = data,
              pred = tb_pred,
              y_pred = y_pred,
              y_prob = y_prob,
              parameters = parameters)
  class(res) <- "LymphoMapObject"
  
  message("Done")
  return(res)
}


#' Print method for LymphoMapObject
#'
#' @param res A \code{LymphoMapObject}.
#' @param ... Additional arguments passed to \code{print()}.
#' @export
print.LymphoMapObject <- function(res, ...) {
  cat("LymphoMAP Model and Prediction\n")
  cat("Accuracy from cross-validation:", round(res$cv_confmat$overall["Accuracy"], 4), "\n")
  cat("Final label assignments:\n")
  if (sum(res$y_pred == "Unassigned") == 0) {
    print(table(factor(res$y_pred, levels = c("FMAC", "LN", "TEX"))))
  } else {
    print(table(res$y_pred))
  }
  invisible(res)
}


