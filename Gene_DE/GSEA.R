
rm(list=ls())
library(SummarizedExperiment)
library(dplyr)

## load microarray
microarray_ovx <- readRDS("data/microarray/processed_microarray.rds")
metadata_microarray <- colData(microarray_ovx) |> as.data.frame()
microarray_feature_data <- rowData(microarray_ovx) |> as.data.frame()
# microarray_feature_data$Gene <- tolower(microarray_feature_data$Gene)
microarray_assay <- assay(microarray_ovx)

lean_mask <- metadata_microarray$Adipo == "l"
obese_mask <- metadata_microarray$Adipo == "ob"

metadata_microarray_lean <- metadata_microarray[lean_mask, ]
metadata_microarray_obese <- metadata_microarray[obese_mask, ]



metadata_lean <- metadata_microarray_lean[rownames(metadata_microarray_lean), ]
metadata_obese <- metadata_microarray_obese[rownames(metadata_microarray_obese), ]



# load the hallmark genes (MSigDB)
hallmarks <- readLines("data/microarray/mh.all.v2025.1.Mm.symbols.gmt")
hallmarks_list <- list()
for (j in 1:length(hallmarks)){
    contents <- strsplit(hallmarks[j], split="\t")[[1]]
    hallmarks_list[[contents[1]]] <- contents[-c(1,2)]
}
combined_hmgenes <- unlist(hallmarks_list) |> unname() |> unique()
existence <- combined_hmgenes %in% microarray_feature_data$Gene
subset_genes <- intersect(combined_hmgenes, microarray_feature_data$Gene)

microarray_feature_hallmark <- microarray_feature_data %>% dplyr::filter(Gene %in% subset_genes)

gene_counts <- table(microarray_feature_hallmark$Gene)



# use limma for differential expression analysis
library(limma)
library(ggplot2)
library(ggrepel)

volcano_plot <- function(lm_result, unit=0.2, pval_cutoff=1e-3, label_col="Gene"){

    lm_result$highlight <- 0
    isannot <- lm_result$P.Value < pval_cutoff
    lm_result$highlight[isannot] <- 1
    lm_result$label <- ""
    lm_result$label[isannot] <- lm_result[isannot, label_col]

    lm_result$highlight <- as.factor(lm_result$highlight)
    lm_result$neg10pval <- -log10(lm_result$P.Value)

    max_y <- ceiling(max(lm_result$neg10pval))

    min_xrange <- floor(min(lm_result$logFC)/unit)*unit
    max_xrange <- ceiling(max(lm_result$logFC)/unit)*unit

    myplot <- ggplot(lm_result, aes(x=logFC, y=neg10pval, color=highlight)) + geom_point() +
        scale_x_continuous(breaks=seq(min_xrange, max_xrange, unit), limits=c(min_xrange, max_xrange)) +
        scale_y_continuous(breaks=seq(0, max_y), limits=c(0, max_y))+
        labs(x="Log2FC", y="-log10 pval")+
        geom_vline(xintercept=0, linetype="dashed") +
        scale_color_manual(values=c("#222222", "#F00A75"))+
        geom_label_repel(aes(label = .data$label),
                         size=2,
                         max.overlaps = 20,
                         box.padding   = 0.35,
                         point.padding = 0.5,
                         segment.color = 'grey50') +
        theme_bw()+
        theme(legend.position = "none",
              axis.title.x = element_blank())

    return(myplot)

}


microarray_lean <- microarray_assay[microarray_feature_hallmark$ID,
                                    rownames(metadata_lean)]
lean_designmat <- cbind(rep(1, nrow(metadata_lean)), 1*(metadata_lean$Treat == "duavee"))
colnames(lean_designmat) <- c("Intercept", "DvsC")
lm_lean <- lmFit(microarray_lean, design=lean_designmat)
lm_lean <- eBayes(lm_lean)
lm_lean_result <- topTable(lm_lean, coef="DvsC", adjust="BH",
                           number=nrow(microarray_feature_hallmark))
lm_lean_result$Gene <- microarray_feature_hallmark[rownames(lm_lean_result), "Gene"]
lean_volcanoplot <- volcano_plot(lm_lean_result)
write.csv(lm_lean_result, "Gene_DE/limma_lean_result.csv",
          quote=FALSE)
ggsave(filename="Gene_DE/limma_lean_result.svg",
       plot=lean_volcanoplot,
       width=5, height=3.5)


microarray_obese <- microarray_assay[microarray_feature_hallmark$ID,
                                     rownames(metadata_obese)]
obese_designmat <- cbind(rep(1, nrow(metadata_obese)), 1*(metadata_obese$Treat == "duavee"))
colnames(obese_designmat) <- c("Intercept", "DvsC")
lm_obese <- lmFit(microarray_obese, design=obese_designmat)
lm_obese <- eBayes(lm_obese)
lm_obese_result <- topTable(lm_obese, coef="DvsC", adjust="BH",
                           number=nrow(microarray_feature_hallmark))
lm_obese_result$Gene <- microarray_feature_hallmark[rownames(lm_obese_result), "Gene"]
obese_volcanoplot <- volcano_plot(lm_obese_result)
write.csv(lm_obese_result, "Gene_DE/limma_obese_result.csv",
          quote=FALSE)
ggsave(filename="Gene_DE/limma_obese_result.svg",
       plot=obese_volcanoplot,
       width=5, height=3.5)


# gene set enrichment analysis
library(fgsea)

lm_lean_result_unique <- lm_lean_result %>% group_by(Gene) %>%
    slice_min(P.Value, with_ties=FALSE)
lean_t_statistic <- lm_lean_result_unique$t
names(lean_t_statistic) <- lm_lean_result_unique$Gene
lean_t_statistic <- sort(lean_t_statistic, decreasing=TRUE)

fgseaRes_lean <- fgsea(
    pathways=hallmarks_list,
    stats=lean_t_statistic,
    minSize=15,
    maxSize=500
)
lean_gsea <- data.frame(Pathway=gsub("HALLMARK_", "", fgseaRes_lean$pathway),
                        P.Value=fgseaRes_lean$pval,
                        pval.adjust=fgseaRes_lean$padj,
                        ES=fgseaRes_lean$ES,
                        NES=fgseaRes_lean$NES,
                        logFC=fgseaRes_lean$NES) %>% arrange(P.Value)
lean_gsea_volcano <- volcano_plot(lm_result=lean_gsea, unit=0.5, label_col="Pathway")
lean_gsea$logFC <- NULL
write.csv(lean_gsea, "Gene_DE/gsea_lean_result.csv",
          row.names=FALSE, quote=FALSE)
ggsave(filename="Gene_DE/lean_gsea_result.svg",
       plot=lean_gsea_volcano,
       width=5, height=3.5)

lm_obese_result_unique <- lm_obese_result %>% group_by(Gene) %>%
    slice_min(P.Value, with_ties=FALSE)
obese_t_statistic <- lm_obese_result_unique$t
names(obese_t_statistic) <- lm_obese_result_unique$Gene
obese_t_statistic <- sort(obese_t_statistic, decreasing=TRUE)

fgseaRes_obese <- fgsea(
    pathways=hallmarks_list,
    stats=obese_t_statistic,
    minSize=15,
    maxSize=500
)
obese_gsea <- data.frame(Pathway=gsub("HALLMARK_", "", fgseaRes_obese$pathway),
                        P.Value=fgseaRes_obese$pval,
                        pval.adjust=fgseaRes_obese$padj,
                        ES=fgseaRes_obese$ES,
                        NES=fgseaRes_obese$NES,
                        logFC=fgseaRes_obese$NES) %>% arrange(P.Value)
obese_gsea_volcano <- volcano_plot(lm_result=obese_gsea, unit=0.5, label_col="Pathway")
obese_gsea$logFC <- NULL
write.csv(obese_gsea, "Gene_DE/gsea_obese_result.csv",
          row.names=FALSE, quote=FALSE)
ggsave(filename="Gene_DE/obese_gsea_result.svg",
       plot=obese_gsea_volcano,
       width=5, height=3.5)

