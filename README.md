# [ConnectivityMap](https://www.science.org/doi/10.1126/science.1132939)
Calculate Connectivity Map scores to identify potential drugs. 

Using the Connectivity Map database, differential gene lists were compared with the transcriptomic effects of over 1,300 known drugs to identify small molecule drugs capable of inducing a similar transcriptional state.

This repository provides an optimized and modularized R pipeline for identifying potential therapeutic compounds or small molecules from a list of differentially expressed genes using the Connectivity Map (CMap) analysis framework.

## Project Structure
```
├── cmap_functions.R             # Core function definitions
├── compute_scores.R             # Main analysis script
├── data/
    ├── DEG.RDS                  # Example input file for DEG
    ├── MSigDB_enrich.RData      # Precalculate msgibdb enrichment to compute specificity score
    ├── GPL96.RData              # For gene to probe ID conversion
├── results/                     # Output directory
└── README.md
```            

## Score
- Connectivity score: with a value between +1 and −1 described the similarity between the target transcriptomic profile (DEG) and an individual compound’s transcriptomic profile, where +1 was most similar and −1 was least similar to the input dataset.
- Specificity_score: Specificity is computed as the proportion of MSigDB gene lists that give a higher connectivity score than the query gene list, which we subtracted from 1.0 so that high scores correspond to higher specificity.
- Reliability score: The reliability measure is "true" (coded as 1) if the perturbagen is represented by more than one profile in CMap and if the majority of experiments show enrichment in the same direction for the query gene list, and is "false" otherwise (coded as 0).
More details, please see [PavlidisLab/regenerationCMap](https://github.com/PavlidisLab/regenerationCMap) and [Nature paper](https://www.nature.com/articles/s41586-025-09647-y).

## Input File Format
Your input file containing the differentially expressed genes must be a CSV file with the following columns:

Column Name	Description
```
gene_name	Official gene symbol (e.g., TP53, STAT1)
log2FoldChange	Log2 fold change of differential expression
pvalue	Statistical p-value from the significance test
padj	Adjusted p-value (e.g., FDR) for multiple test correction
```

Example (DEG.RDS):
```
gene_name,log2FoldChange,pvalue,padj
STAT1,3.2,0.00001,0.0005
CDKN1A,2.5,0.00002,0.0006
```

## Usage
1. Run Analysis: Execute the compute_scores.R script in R.
2. Output
  - Ranked List of Potential Compounds (ranks.RData):

      A CSV file listing all candidate small molecules/compounds, ranked by their association (Connectivity Score) with your input gene signature.
  
      Includes all calculated scores (Connectivity, Specificity, Reliability) for easy filtering and prioritization.
  
  - Visualization Plot (dotplot_cmap.pdf):
  
      A bar chart visualizing the top candidate compounds and their Connectivity Scores.
  
      Typically, strong positive scores suggest the compound may induce a biological state similar to your signature, while strong negative scores suggest it may reverse the state.

## Dependencies
Please ensure the following R packages are installed before running the scripts:
```
ConnectivityMap
memoise
ggrepel
dplyr
magrittr
gglue
ggplot2
```

You can install them using:

r
install.packages(c("dplyr", "ggplot2", "readr"))



