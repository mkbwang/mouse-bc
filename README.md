# mouse-bc
GSRA Project Analyzing BZA/CE treatment.

* `EDA.R` explores the microbiome data and plots the relative abundance values on phyla and species level for different groups of mice (ovx vs intact, duavee vs control).
* `DAA.R` carry out differential abundance analysis for microbiome species between different groups of samples (ovx vs intact, duavee vs control).
* `GSEA.R`  carries out differential expression analysis and gene set enrichment analysis based on the microarray data. Limit to ovx mice only.
* `gene_microbiome_correlation.R` focuses on differentially expressed genes that belong to enriched pathways. Then I calculate the spearman correlation between gene expression and microbiome species relative abundance.
