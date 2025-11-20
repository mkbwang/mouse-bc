

rm(list=ls())
# load gene expression data and gene information


library(SummarizedExperiment)
library(dplyr)
library(ComplexHeatmap)

## load microarray
microarray_ovx <- readRDS("data/microarray/processed_microarray.rds")
metadata_microarray <- colData(microarray_ovx) |> as.data.frame()
microarray_feature_data <- rowData(microarray_ovx) |> as.data.frame()

microarray_assay <- assay(microarray_ovx)
lean_mask <- metadata_microarray$Adipo == "l"
obese_mask <- metadata_microarray$Adipo == "ob"

metadata_microarray_lean <- metadata_microarray[lean_mask, ]
metadata_microarray_obese <- metadata_microarray[obese_mask, ]


## load microbiome data
microbiome_data <- read.csv("data/microbiome/microbiome_data.csv", row.names=1)
metadata_microbiome <- microbiome_data[, seq(1,3)]
microbiome_assay <- microbiome_data[, -seq(1,9)]

metadata_lean <- metadata_microarray_lean[intersect(rownames(metadata_microarray_lean),
                                                    rownames(metadata_microbiome)), ]
metadata_obese <- metadata_microarray_obese[intersect(rownames(metadata_microarray_obese),
                                                      rownames(metadata_microbiome)), ]


microarray_assay_lean <- microarray_assay[, rownames(metadata_lean)]
microarray_assay_obese <- microarray_assay[, rownames(metadata_obese)]
microbiome_assay_lean <- microbiome_assay[rownames(metadata_lean), ]
microbiome_assay_obese <- microbiome_assay[rownames(metadata_obese), ]



## load limma and GSEA results and focus on the differentially expressed genes only

limma_lean_result <- read.csv("data/microarray/DE/limma_lean_result.csv")
limma_obese_result <- read.csv("data/microarray/DE/limma_obese_result.csv")
gsea_lean_result <- read.csv("data/microarray/DE/gsea_lean_result.csv")
gsea_obese_result <- read.csv("data/microarray/DE/gsea_obese_result.csv")



## load the matching between genes and pathways

hallmarks <- readLines("data/microarray/mh.all.v2025.1.Mm.symbols.gmt")
hallmarks_list <- list()
gene_hallmark_match_list <- list()
for (j in 1:length(hallmarks)){
    contents <- strsplit(hallmarks[j], split="\t")[[1]]
    pathway_name <- gsub("HALLMARK_", "", contents[1])
    hallmarks_list[[pathway_name]] <- contents[-c(1,2)]
    gene_hallmark_match_list[[pathway_name]] <- data.frame(Gene=contents[-c(1,2)], Pathway=pathway_name)
}
gene_hallmark_match <- do.call(rbind, gene_hallmark_match_list)



# only keep the genes that belong to significant gene sets
## p value adjustment within each pathway separately
gsea_lean_subset <- gsea_lean_result %>% filter(pval.adjust < 0.05)
lean_genes_list <- list()
for (j in 1:nrow(gsea_lean_subset)){
    pathway_name <- gsea_lean_subset$Pathway[j]
    selected_genes_df <- limma_lean_result %>% filter(Gene %in% hallmarks_list[[pathway_name]])
    selected_genes_df$adj.P.Val <- p.adjust(selected_genes_df$P.Value, method="BH")
    lean_genes_list[[pathway_name]] <- selected_genes_df
}
lean_marker_genes <- list()
for (j in 1:length(lean_genes_list)){
    pathway_name <- names(lean_genes_list)[j]
    selected_genes_df <- lean_genes_list[[pathway_name]]
    lean_marker_genes[[pathway_name]] <- selected_genes_df %>% filter(adj.P.Val < 0.05)
}
lean_marker_genes_combined <- do.call(rbind, lean_marker_genes)
limma_lean_markergenes <- limma_lean_result %>% filter(P.Value <= max(lean_marker_genes_combined$P.Value))


gsea_obese_subset <- gsea_obese_result %>% filter(pval.adjust < 0.05)
obese_genes_list <- list()
for (j in 1:nrow(gsea_obese_subset)){

    pathway_name <- gsea_obese_subset$Pathway[j]
    selected_genes_df <- limma_obese_result %>% filter(Gene %in% hallmarks_list[[pathway_name]])
    selected_genes_df$adj.P.Val <- p.adjust(selected_genes_df$P.Value, method="BH")
    obese_genes_list[[pathway_name]] <- selected_genes_df

}
obese_marker_genes <- list()
for (j in 1:length(obese_genes_list)){
    pathway_name <- names(obese_genes_list)[j]
    selected_genes_df <- obese_genes_list[[pathway_name]]
    obese_marker_genes[[pathway_name]] <- selected_genes_df %>% filter(adj.P.Val < 0.05)
}
obese_marker_genes_combined <- do.call(rbind, obese_marker_genes)
limma_obese_markergenes <- limma_obese_result %>% filter(P.Value <= max(obese_marker_genes_combined$P.Value))




microarray_assay_lean_markers <- microarray_assay_lean[limma_lean_markergenes$X, ]
microarray_assay_obese_markers <- microarray_assay_obese[limma_obese_markergenes$X, ]


prevalence_microbiome_lean <- colMeans(microbiome_assay_lean > 0)
microbiome_assay_lean_subset <- microbiome_assay_lean[, prevalence_microbiome_lean >= 0.8]
prevalence_microbiome_obese <- colMeans(microbiome_assay_obese > 0)
microbiome_assay_obese_subset <- microbiome_assay_obese[, prevalence_microbiome_obese >= 0.8]



corr_mat_lean <- matrix(0, nrow=nrow(limma_lean_markergenes), ncol=ncol(microbiome_assay_lean_subset))
pval_mat_lean <- matrix(0, nrow=nrow(limma_lean_markergenes), ncol=ncol(microbiome_assay_lean_subset))
rownames(corr_mat_lean) <- rownames(pval_mat_lean) <- limma_lean_markergenes$Gene
colnames(corr_mat_lean) <- colnames(pval_mat_lean) <- colnames(microbiome_assay_lean_subset)
for (i in 1:nrow(corr_mat_lean)){

    for(j in 1:ncol(corr_mat_lean)){

        test_result <- cor.test(microarray_assay_lean_markers[i, ], microbiome_assay_lean_subset[, j],
                                method="spearman")
        corr_mat_lean[i, j] <- test_result$estimate
        pval_mat_lean[i, j] <- test_result$p.value

    }

}
significance_count_lean <- colSums(pval_mat_lean < 0.05)
gene_hallmark_match_lean <- gene_hallmark_match %>% filter(Gene %in% rownames(corr_mat_lean)) %>%
    arrange(Pathway)
pval_mat_lean_subset <- pval_mat_lean[unique(gene_hallmark_match_lean$Gene), significance_count_lean > 0]
corr_mat_lean_subset <- corr_mat_lean[unique(gene_hallmark_match_lean$Gene), significance_count_lean > 0]
bin_mat_lean_subset <- 1*(corr_mat_lean_subset > 0) * (pval_mat_lean_subset < 0.05) -
    1*(corr_mat_lean_subset < 0) * (pval_mat_lean_subset < 0.05)
write.csv(gene_hallmark_match_lean, "data/microarray/correlation/gene_pathway_lean.csv", row.names=FALSE)
write.csv(t(corr_mat_lean_subset), "data/microarray/correlation/correlation_lean.csv", row.names=T, quote=FALSE)

library(circlize)
col_fun <- colorRamp2(c(-1, 0, 1), c("blue", "white", "red"))
ht_corr_lean <- Heatmap(t(corr_mat_lean_subset), name="Correlation", cluster_rows=FALSE, cluster_columns=FALSE,
              show_column_names = TRUE, show_row_names = TRUE,
              show_heatmap_legend = T,
              border_gp = gpar(col="black", lwd=0.5),
              rect_gp = gpar(col = "gray", lwd = 0.5),
              width = unit(0.8, "npc"),
              height = unit(0.8, "npc"),
              row_names_side = "left",
              column_names_side="top",
              na_col = "#555555",
              col = col_fun)


ht_bin_lean <- Heatmap(t(bin_mat_lean_subset), name="Correlation", cluster_rows=FALSE, cluster_columns=FALSE,
                        show_column_names = TRUE, show_row_names = TRUE,
                        show_heatmap_legend = FALSE,
                        border_gp = gpar(col="black", lwd=0.5),
                       rect_gp = gpar(col = "gray", lwd = 0.5),
                        width = unit(0.8, "npc"),
                        height = unit(0.8, "npc"),
                        row_names_side = "left",
                        column_names_side="top",
                        na_col = "#555555",
                        col = col_fun)


corr_mat_obese <- matrix(0, nrow=nrow(limma_obese_markergenes), ncol=ncol(microbiome_assay_obese_subset))
pval_mat_obese <- matrix(0, nrow=nrow(limma_obese_markergenes), ncol=ncol(microbiome_assay_obese_subset))
rownames(corr_mat_obese) <- rownames(pval_mat_obese) <- limma_obese_markergenes$Gene
colnames(corr_mat_obese) <- colnames(pval_mat_obese) <- colnames(microbiome_assay_obese_subset)
for (i in 1:nrow(corr_mat_obese)){
    for(j in 1:ncol(corr_mat_obese)){
        test_result <- cor.test(microarray_assay_obese_markers[i, ], microbiome_assay_obese_subset[, j],
                                method="spearman")
        corr_mat_obese[i, j] <- test_result$estimate
        pval_mat_obese[i, j] <- test_result$p.value
    }
}
significance_count_obese <- colSums(pval_mat_obese < 0.05)
gene_hallmark_match_obese <- gene_hallmark_match %>% filter(Gene %in% rownames(corr_mat_obese)) %>%
    arrange(Pathway)
pval_mat_obese_subset <- pval_mat_obese[unique(gene_hallmark_match_obese$Gene), significance_count_obese > 0]
corr_mat_obese_subset <- corr_mat_obese[unique(gene_hallmark_match_obese$Gene), significance_count_obese > 0]
bin_mat_obese_subset <- 1*(corr_mat_obese_subset > 0) * (pval_mat_obese_subset < 0.05) -
    1*(corr_mat_obese_subset < 0) * (pval_mat_obese_subset < 0.05)
ht_corr_obese <- Heatmap(t(corr_mat_obese_subset), name="Correlation", cluster_rows=FALSE, cluster_columns=FALSE,
                        show_column_names = TRUE, show_row_names = TRUE,
                        show_heatmap_legend = T,
                        border_gp = gpar(col="black", lwd=0.5),
                        rect_gp = gpar(col = "gray", lwd = 0.5),
                        width = unit(0.8, "npc"),
                        height = unit(0.8, "npc"),
                        row_names_side = "left",
                        column_names_side="top",
                        na_col = "#555555",
                        col = col_fun)
ht_bin_obese <- Heatmap(t(bin_mat_obese_subset), name="Correlation", cluster_rows=FALSE, cluster_columns=FALSE,
                       show_column_names = TRUE, show_row_names = TRUE,
                       show_heatmap_legend = FALSE,
                       border_gp = gpar(col="black", lwd=0.5),
                       rect_gp = gpar(col = "gray", lwd = 0.5),
                       width = unit(0.8, "npc"),
                       height = unit(0.8, "npc"),
                       row_names_side = "left",
                       column_names_side="top",
                       na_col = "#555555",
                       col = col_fun)

write.csv(gene_hallmark_match_obese,
          "data/microarray/correlation/gene_pathway_obese.csv", row.names=FALSE)
write.csv(t(corr_mat_obese_subset),
          "data/microarray/correlation/correlation_obese.csv", row.names=T, quote=FALSE)



