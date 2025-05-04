#' plot_umap
#'
#' @description
#' Plot UMAP
#'
#' @param mt sample by gene expression matrix. 
#' @param sample_info sample categories to be shown on the scatter plot.
#' @param npcs number of principal components used for UMAP calculation; if NULL, calculate UMAP using the input matrix without PCA. Default = NULL.
#' @param palette palette of sample categories. Default = NULL
#' 
#' @details
#' Calculate and generate a UMAP plot.
#'
#' @return a ggplot object.
#' 
#' @export
plot_umap <- function(mt, sample_info, npcs = NULL, palette = NULL) {
  stopifnot(all.equal(names(sample_info), rownames(mt)))
  
  # pca
  if (!is.null(npcs)) {
    if (!is.numeric(npcs)) {
      stop("`npcs` must be numeric.")
    }
    mt <- scale(mt)
    pca_res <- prcomp(mt, rank = npcs)
    mt <- pca_res$x[, 1:npcs]
  }
  
  # umap
  d <- as.dist(1 - cor(t(mt)))
  set.seed(42)
  umap2d <- umap(d,
                 # metric = "correlation",
                 min_dist = 1,
                 n_neighbors = 20,
                 dens_scale = 0.1,
                 # init = "spca",
                 init_sdev = "range")
  
  colnames(umap2d) <- c("UMAP1", "UMAP2")
  tb <- tibble(sample_id = rownames(umap2d)) %>%
    bind_cols(as_tibble(umap2d)) %>%
    mutate(sample_info = sample_info)
  
  p <- ggscatter(tb, "UMAP1", "UMAP2",
                 color = "sample_info",
                 fill =  "sample_info",
                 palette = palette)
  return(p)
}


#' plot_umap_archetype
#'
#' @description
#' Plot UMAP by archetypes
#'
#' @param mt sample by gene expression matrix. 
#' @param sample_archetypes sample archetypes.
#' @param npcs number of principal components used for UMAP calculation; if NULL, calculate UMAP using the input matrix without PCA. Default = NULL.
#' 
#' @details
#' Calculate and generate a UMAP plot showing archetypes.
#'
#' @return a ggplot object.
#' 
#' @export
plot_umap_archetype <- function(mt, sample_archetypes, npcs = NULL) {
  col_archetype <- structure(
    c("#377eb8", "#4daf4a", "#e41a1c", "#BBBBBB"),
    names = c("FMAC", "LN", "TEX", "Unassigned")
  )
  return(plot_umap(mt, sample_archetypes, npcs = npcs, palette = col_archetype))
}


#' plot_confmat
#'
#' @description
#' Plot heatmap of confusion matrix
#'
#' @param confmat a confusion matrix object returned by caret::confusionMatrix().
#' 
#' @details
#' Plot heatmap of confusion matrix returned by caret::confusionMatrix().
#'
#' @return a Heatmap object returned by ComplexHeatmap::draw().
#' 
#' @export
plot_confmat <- function(confmat) {
  mt <- confmat$table %>%
    as.matrix() %>%
    apply(2, as.integer)
  rownames(mt) <- rownames(confmat$table)
  
  # plot
  col_fun <- colorRamp2(c(0, max(mt)), c("white", "red"))
  ht <- Heatmap(mt, name = "N",  
                col = col_fun,
                row_title = "Prediction",
                column_title = "Reference",
                row_names_side = "left",
                column_names_side = "top",
                rect_gp = gpar(col = "black"),
                cluster_rows = FALSE,
                cluster_columns = FALSE,
                cell_fun = function(j, i, x, y, width, height, fill) {
                  if(mt[i, j] > 0)
                    grid.text(sprintf("%.0f", mt[i, j]), x, y, gp = gpar(fontsize = 10))
                })
  p <- draw(ht,
            newpage = FALSE)
  return(NULL)
}


