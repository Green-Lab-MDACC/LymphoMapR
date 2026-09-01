# LymphoMapR

**LymphoMapR** provides an R implementation of the LymphoMAP classification framework for large B-cell lymphoma (LBCL), using bulk gene expression data to classify tumors based on their microenvironmental archetypes. These archetypes reflect distinct, functionally co-occurring microenvironmental cell subsets that are associated with disease progression and treatment response.

## Installation

```r
# install.packages("remotes") # if not already installed
remotes::install_github("Green-Lab-MDACC/LymphoMapR")
```

## Example

The main function in **LymphoMapR** is `run_LymphoMap()`, which performs archetype classification based on a gene expression matrix. The input should be a **gene-by-sample matrix** (e.g., log2-transformed TPM values). The function integrates your input with an internal training dataset and applies a naive Bayes model to assign each sample to one of **LymphoMAP archetypes**.

```r
library(LymphoMapR)

# expr_matrix: your gene-by-sample expression matrix (gene symbols as row names, sample names as column names)
result <- run_LymphoMap(expr_matrix)

# View classification results
print(result)

# Predicted labels and probabilities
result$pred
```

Classification accuracy is evaluated using cross-validation of the training data in LymphoMapR. To verify whether the predicted archetypes capture the expected biology, users should perform a post-classification differential expression analysis among the three archetypes. The expected differentially expressed genes (DEGs) are highlighted in Figures S7A and S8A of the LymphoMAP paper (Li et al., Cancer Cell, 2025).

## Citation

If you use **LymphoMapR** in your research, please cite:

> Li X, Singhal K, Deng Q, Chihara D, Russler-Germain D, Harkins RA, Henderson J, Arita K, Kizhakeyil A, Sun R, Lakra P, Hussein U, Foltz JA, Wilson A, Schmidt E, Nizamuddin I, Dinh T, Kesaraju A, Hamilton MP, Allen C, Gandhi MK, Tobin J, Jiang A, Hilton L, Scott DW, Vega F, Flowers CR, Westin JR, Griffith OL, Fehniger TA, Griffith M, Green MR. Large B cell lymphoma microenvironment archetype profiles. Cancer Cell. 2025 Jul 14;43(7):1347-1364.e13.

