# feature selection
# for each gene, get fraction of samples with zscore > 0
get_fract_poszscore_group_rest <- function(mt, y) {
  stopifnot(all.equal(colnames(mt), names(y)))
  
  # zscore
  mt <- t(mt)
  mt <- scale(mt)
  
  # wide to long and get poszscore
  tb <- tibble(sample_id = rownames(mt)) %>%
    bind_cols(as_tibble(mt)) %>%
    gather("gene_name", "zscore", -sample_id) %>%
    mutate(poszscore = (zscore > 0))
  
  # labels
  tb <- tb %>%
    left_join(y %>% enframe("sample_id", "LymphoMAP"),
              by = "sample_id") 
  tb_data <- tb
  
  # get fract for each archetype
  tb_fract <- tibble()
  for (s in unique(tb_data$LymphoMAP)) {
    tb <- tb_data %>%
      mutate(cat = ifelse(LymphoMAP == s,
                          "group",
                          "rest"))
    
    # calc fract
    tb_tot <- tb %>%
      group_by(cat, gene_name) %>%
      summarise(n_tot = n(), .groups = "drop") %>%
      ungroup()
    tb <- tb %>%
      group_by(cat, gene_name) %>%
      summarise(n_poszscore = sum(poszscore), .groups = "drop") %>%
      ungroup() %>%
      left_join(tb_tot, by = c("cat", "gene_name")) %>%
      mutate(fract_poszscore = n_poszscore/n_tot) %>%
      select(-n_poszscore, -n_tot)
    
    # long to wide
    tb <- tb %>%
      spread(cat, fract_poszscore) %>%
      dplyr::rename(fract_poszscore_group = group,
                    fract_poszscore_rest = rest) %>%
      mutate(LymphoMAP = s) %>%
      select(LymphoMAP, everything())
    
    tb_fract <- tb_fract %>%
      bind_rows(tb)
  }
  
  return(tb_fract)
}


# DEGs between archetypes in train data
get_deg_archetype <- function(gene_pool, exclude_noncoding_MtRp = TRUE) {
  stopifnot(length(setdiff(gene_pool, rownames(train_data))) == 0)
  
  # intersect gene_pool with protein_coding_nonMtRp
  if (exclude_noncoding_MtRp) {
    pcg_nonMtRp <- rowData(train_data) %>% 
      as_tibble() %>% 
      dplyr::filter(protein_coding_nonMtRp == 1) %>% 
      .$gene_name
    gene_pool <- intersect(gene_pool, pcg_nonMtRp)
  }
  
  # construct a DESeqDataSet
  df_meta <- colData(train_data) %>%
    as.data.frame() %>%
    dplyr::mutate(condition = LymphoMAP) # condition as the design variable to select group of interest in DEG analysis
  dds <- DESeqDataSetFromMatrix(countData = assay(train_data[gene_pool, ], "counts"),
                                colData = df_meta,
                                design = ~ condition)
  
  # pre-filtering to keep only rows that have a count of at least 10 for a minimal number of samples.
  smallestGroupSize <- 3
  keep <- rowSums(counts(dds) >= 10) >= smallestGroupSize
  dds <- dds[keep,]
  
  # DEGs between archetypes
  archetypes <- levels(dds$LymphoMAP)
  dds_ori <- dds
  
  tb_deg <- tibble()
  for (s in archetypes) {
    message("calculating DEGs of ", s, " vs rest")
    # deg
    dds <- dds_ori
    dds$condition <- ifelse(dds$LymphoMAP == s,
                            s,
                            "rest")
    dds$condition <- factor(dds$condition, levels = c("rest", s))
    dds <- DESeq(dds)
    
    # The function DESeq runs the following functions in order:
    # dds <- estimateSizeFactors(dds)
    # dds <- estimateDispersions(dds)
    # dds <- nbinomWaldTest(dds)
    
    res <- results(dds)
    tb <- tibble(LymphoMAP = s,
                 gene_name = rownames(res)) %>%
      bind_cols(as_tibble(res)) %>%
      arrange(desc(sign(log2FoldChange)*(-log(padj,10))), desc(log2FoldChange)) 
    # tb <- tb %>%
    #   filter(log2FoldChange > 0.58 &
    #            padj < 0.05)
    tb_deg <- tb_deg %>%
      bind_rows(tb)
  }
  
  # get fraction of positive zscore in group and rest 
  mt <- as.matrix(assay(train_data, "logTPM"))
  mt <- mt[gene_pool, ]
  y <- train_data$LymphoMAP
  names(y) <- train_data$sample_id
  tb_fract_poszscore <- get_fract_poszscore_group_rest(mt, y)
  
  tb_deg <- tb_deg %>%
    left_join(tb_fract_poszscore, by = c("LymphoMAP", "gene_name")) 
  
  return(tb_deg)
}


# filter top genes and further remove genes shared by multiple archetypes
get_top_unique_genes <- function(n_top, tb_deg) {
  tb <- tb_deg %>%
    filter(log2FoldChange > 0.58 &
             padj < 0.05) %>%
    arrange(desc(stat)) %>%
    group_by(LymphoMAP) %>%
    do(head(., n = n_top)) %>%
    ungroup()
  
  dup_genes <- tb$gene_name[duplicated(tb$gene_name)]
  tb <- tb %>%
    filter(!(gene_name %in% dup_genes))
  return(tb)
}


