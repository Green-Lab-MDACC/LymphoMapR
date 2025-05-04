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

## Citation

If you use **LymphoMapR** in your research, please cite:

> Xubin Li, Kartik Singhal, Qing Deng, Dai Chihara, David A Russler-Germain, Usama Hussein, Jennifer Ann Foltz, Jared Henderson, Ashley Wilson, Evelyn Schmidt, Imran A. Nizamuddin, Tommy Dinh, Ryan Sun, Akhil Kesaraju, Laura K. Hilton, David W. Scott, Francisco Vega, Christopher R. Flowers, Obi Griffith, Todd A Fehniger, Malachi Griffith, Michael R. Green, Comprehensive Characterization and Validation of the Tumor Microenvironment in Patients with Relapsed/Refractory Large B-Cell Lymphoma Identifies Subgroups with Greatest Benefit from CD19 CAR T-Cell Therapy, Blood, Volume 144, Supplement 1, 2024. https://doi.org/10.1182/blood-2024-204166.

