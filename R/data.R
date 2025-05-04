#' Training Dataset
#'
#' A `SummarizedExperiment` object of RNA-seq data and LymphoMAP labels, with two assays, column metadata, and row metadata.
#'
#' @docType data
#' @name train_data
#' @usage data(train_data)
#' 
#' @format A `SummarizedExperiment` object with:
#' \describe{
#'   \item{assays}{A list of two assays:
#'     \describe{
#'       \item{\code{"counts"}}{Expression counts matrix.}
#'       \item{\code{"logTPM"}}{Expression matrix of log2(TPM + 1).}
#'     }
#'   }
#'   \item{colData}{A `DataFrame` of sample-level metadata:
#'     \describe{
#'       \item{\code{sample_id}}{Sample IDs.}
#'       \item{\code{LymphoMAP}}{A factor indicating LymphoMAP labels.}
#'     }
#'   }
#'   \item{rowData}{A `DataFrame` of gene-level metadata:
#'     \describe{
#'       \item{\code{gene_name}}{Gene symbols.}
#'       \item{\code{protein_coding_nonMtRp}}{A binary indicator (1 or 0) showing whether a gene is protein-coding and not mitochondrial or ribosomal.}
#'     }
#'   }
#' }
#'
#' @source Internal
NULL
