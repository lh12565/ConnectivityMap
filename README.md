# [ConnectivityMap](https://www.science.org/doi/10.1126/science.1132939)
Calculate Connectivity Map (CMap) scores to identify potential drugs ([Science paper](https://www.science.org/doi/10.1126/science.1132939)). 

Using the Connectivity Map database, differential gene lists were compared with the transcriptomic effects of over 1,300 known drugs to identify small molecule drugs capable of inducing a similar transcriptional state.

This repository provides an optimized and modularized R pipeline for identifying potential drugs from a list of differentially expressed genes using the Connectivity Map analysis framework.

## Dependencies
Please ensure the following R packages are installed before running the scripts:
- ConnectivityMap
- memoise
- ggrepel
- dplyr
- magrittr
- glue
- ggplot2

Use the script below to automatically check and install any missing packages from both CRAN and Bioconductor：

```r
# Required packages
packages <- c(
  "memoise","ggrepel","dplyr","magrittr","glue","ggplot2","ConnectivityMap"
)

# Function to install CRAN packages
install_cran <- function(pkg) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg, dependencies = TRUE)
  }
}

# Function to install Bioconductor packages
install_bioc <- function(pkg) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    if (!requireNamespace("BiocManager", quietly = TRUE)) {
      install.packages("BiocManager")
    }
    BiocManager::install(pkg, dependencies = TRUE)
  }
}

# Install all packages
for (pkg in packages) {
  if (pkg %in% c("ConnectivityMap")) {
    install_bioc(pkg)
  } else {
    install_cran(pkg)
  }
}
```

## Project Structure
```
├── cmap_functions.R             # Core function definitions
├── compute_scores.R             # Main analysis script
├── data/
    ├── DEG.RDS                  # Example input file for DEG
    ├── MSigDB_enrich.RData      # Precalculate msgibdb enrichment to compute specificity score
    ├── GPL96.RData              # For gene to probe ID conversion
├── results/                     # Output directory
    ├── out_cmap.RData           # Score file
    ├── ranks.RData              # Filtered and ranked score file
    ├── dotplot_cmap.pdf         # dotplot for each ranked instance
└── README.md
```            

## Score
- **Connectivity score**: with a value between +1 and −1 described the similarity between the target transcriptomic profile (DEG) and an individual compound’s transcriptomic profile, where +1 was most similar and −1 was least similar to the input dataset.
- **Specificity_score**: Specificity is computed as the proportion of MSigDB gene lists that give a higher connectivity score than the query gene list, which we subtracted from 1.0 so that high scores correspond to higher specificity.
- **Reliability score**: The reliability measure is "true" (coded as 1) if the perturbagen is represented by more than one profile in CMap and if the majority of experiments show enrichment in the same direction for the query gene list, and is "false" otherwise (coded as 0).

More details, please see [PavlidisLab/regenerationCMap](https://github.com/PavlidisLab/regenerationCMap) and [Nature paper](https://www.nature.com/articles/s41586-025-09647-y).

## Input File Format
Your input file containing the differentially expressed genes must be a data.frame with the following columns:

Column Name  Description
```
gene_name            Official gene symbol (e.g., TP53, STAT1)
log2FoldChange       Log2 fold change of differential expression
pvalue               Statistical p-value from the significance test
padj                 Adjusted p-value (e.g., FDR) for multiple test correction
```

Example (DEG.RDS):
```
                mean_c     mean_e log2FoldChange   pvalue    padj    gene_name
MAPKAPK5-AS1  5.743735 10.0778887      0.8111324 2.54e-07 0.00190 MAPKAPK5-AS1
ZNF93         6.469094  0.8937912     -2.8555539 2.38e-07 0.00190        ZNF93
PML          59.666754 22.3601812     -1.4159954 5.15e-07 0.00225          PML
IFIT3        45.767544  3.9143528     -3.5474792 7.51e-07 0.00225        IFIT3
BTN2A2       18.058070  6.9460147     -1.3783863 7.05e-07 0.00225       BTN2A2
ADRB2        13.122826  0.8066832     -4.0239324 9.05e-07 0.00226        ADRB2
```

## Usage
1. Run Analysis:

      Execute the compute_scores.R script in R.

      Filtering and sorting, can be appropriately modified according to your own data
3. Output
  - Score file (out_cmap.RData)

      Includes all calculated scores (Connectivity, Specificity, Reliability) for easy filtering and prioritization.
  - Ranked List of Potential Compounds (ranks.RData):

      A file listing all candidate small molecules/compounds, ranked by their association (Scores) with your input gene signature.
  - Visualization Plot (dotplot_cmap.pdf):
  
      A dot plot visualizing the top candidate compounds and their Connectivity Scores.
  
      Typically, strong positive scores suggest the compound may induce a biological state similar to your signature, while strong negative scores suggest it may reverse the state.




