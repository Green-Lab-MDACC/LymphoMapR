#' fpkm_to_tpm
#'
#' @description
#' convert FPKM to TPM
#'
#' @param fpkm gene by sample matrix in FPKM. 
#' 
#' @details
#' Convert Fragments Per Kilobase of transcript per Million mapped reads (FPKM) to Transcripts Per Million (TPM).
#'
#' @return gene by sample matrix in TPM. 
#' 
#' @export
fpkm_to_tpm <- function(fpkm) {
  sum_fpkm <- colSums(fpkm)  # Sum of all FPKMs per sample
  tpm <- sweep(fpkm, 2, sum_fpkm, FUN = "/") * 1e6  # Normalize
  return(tpm)
}


# convert counts to vst
counts2vst <- function(cts) {
  # convert cts to integers that is required by DESeq2
  genes <- rownames(cts)
  cts <- apply(cts, 2, round)
  cts <- apply(cts, 2, as.integer)
  rownames(cts) <- genes
  
  # construct a DESeqDataSet
  tb <- tibble(sample_id = colnames(cts)) %>%
    mutate(condition = paste0("cond", c(1:ncol(cts)) %% 2) %>% 
             as.factor())
  coldata <- tb[, -1] %>% as.data.frame()
  rownames(coldata) <- tb[[1]]
  dds <- DESeqDataSetFromMatrix(countData = cts,
                                colData = coldata,
                                design = ~ condition)
  
  # pre-filtering to keep only rows that have a count of at least 10 for a minimal number of samples.
  smallestGroupSize <- 3
  keep <- rowSums(counts(dds) >= 10) >= smallestGroupSize
  dds <- dds[keep,]
  
  # Extracting transformed values
  # Variance stabilizing transformation (The transformed data is on the log2 scale)
  vsd <- varianceStabilizingTransformation(dds, blind = TRUE)
  
  return(assay(vsd))
}


# align pred expr matrix to train expr matrix
combat_pred_to_train <- function(mt_train, mt_pred) {
  stopifnot(length(intersect(colnames(mt_train), colnames(mt_pred))) == 0)
  stopifnot(all.equal(rownames(mt_train), rownames(mt_pred)))
  
  batch <- c(rep("train", ncol(mt_train)), 
             rep("pred", ncol(mt_pred)))
  mt_comb <- cbind(mt_train, mt_pred)
  mt_adj <- ComBat(mt_comb, batch = batch, ref.batch = "train")
  return(mt_adj[, colnames(mt_pred)])
}


# scale train expr matrix and apply mean/scale to pred expr matrix
scale_expr <- function(mt_train, mt_pred) {
  stopifnot(all.equal(colnames(mt_train), colnames(mt_pred)))
  
  # scale train
  mt_train <- scale(mt_train)
  x_train_center <- attr(mt_train, "scaled:center")
  x_train_scale <- attr(mt_train, "scaled:scale")
  
  # scale pred to match train
  stopifnot(all.equal(colnames(mt_pred), names(x_train_center)))
  mt_pred <- scale(mt_pred, center = x_train_center, scale = x_train_scale)
  
  return(list(train = mt_train,
              pred = mt_pred))
}


