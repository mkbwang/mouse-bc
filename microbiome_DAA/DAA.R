

library(ADAPT)
library(phyloseq)
library(ggplot2)
library(ggrepel)

rm(list=ls())
mouse_microbiome <- read.csv("data/microbiome/microbiome_data.csv", row.names=1)
mouse_meta <- mouse_microbiome[, seq(1,3)]

mouse_metagenomics <- mouse_microbiome[, -seq(1,3)]
mouse_phyla <- mouse_metagenomics[, 1:6]
mouse_species <- mouse_metagenomics[, 7:118]
physeq_all <- phyloseq(otu_table(as.matrix(mouse_species), taxa_are_rows = FALSE),
                       sample_data(mouse_meta))
min_positive_val <- min(mouse_species[mouse_species > 0])
DAA_ovary <- adapt(physeq_all, cond.var="ovx_or_intact",
                           base.cond="intact", adj.var=c("Rx", "Adipo"),
                   depth.filter=0,
                           censor=min_positive_val,
                           prev.filter=0.05)
ovx_vs_intact <- DAA_ovary@details
sorted_pvals <- sort(ovx_vs_intact$pval)
sorted_effsizes <- sort(ovx_vs_intact$log10foldchange)
pval_mask <- ovx_vs_intact$pval <= sorted_pvals[5]
effsize_mask <- ovx_vs_intact$log10foldchange <= sorted_effsizes[5] |
    ovx_vs_intact$log10foldchange >= tail(sorted_effsizes, 5)[1]
ovx_vs_intact$Labels <- ""
ovx_vs_intact$Labels[pval_mask | effsize_mask] <- ovx_vs_intact$Taxa[pval_mask | effsize_mask]
ovx_vs_intact$labeled <- pval_mask | effsize_mask
ovx_vs_intact$neglog10pval <- -log10(ovx_vs_intact$pval)

write.csv(ovx_vs_intact[, 1:6], "plots/DAA_ovary.csv",
          row.names=FALSE, quote=FALSE)

ovx_vs_intact_plot <- ggplot(ovx_vs_intact, aes(x=log10foldchange, y=neglog10pval))+
    geom_point(alpha=0.8, aes(color=.data$labeled)) +
    xlab("") + ylab("-Log10 p-value") + theme_bw() +
    scale_x_continuous(limits=c(-2, 2), breaks=seq(-2, 2, 0.4))+
    scale_y_continuous(limits=c(0, 6), breaks=seq(0, 6))+
    theme(legend.position="none", axis.title=element_text(size=10),
          axis.text=element_text(size=10)) +
    scale_color_manual(values=c("#616161", "#ff0066")) +
    geom_vline(xintercept=0, linetype="dashed", color = "blue") +
    geom_label_repel(aes(label = Labels),
                     size=2,
                     max.overlaps = 20,
                     box.padding   = 0.35,
                     point.padding = 0.5,
                     segment.color = 'grey50')






DAA_trt <- adapt(physeq_all, cond.var="Rx",
                 base.cond="control",
                 depth.filter=0, adj.var=c("ovx_or_intact", "Adipo"),
                 censor=min_positive_val,
                 prev.filter=0.05)

duavee_vs_control <- DAA_trt@details
sorted_pvals <- sort(duavee_vs_control$pval)
sorted_effsizes <- sort(duavee_vs_control$log10foldchange)
pval_mask <- duavee_vs_control$pval <= sorted_pvals[5]
effsize_mask <- duavee_vs_control$log10foldchange <= sorted_effsizes[5] |
    duavee_vs_control$log10foldchange >= tail(sorted_effsizes, 5)[1]
duavee_vs_control$Labels <- ""
duavee_vs_control$Labels[pval_mask | effsize_mask] <- duavee_vs_control$Taxa[pval_mask | effsize_mask]
duavee_vs_control$labeled <- pval_mask | effsize_mask
duavee_vs_control$neglog10pval <- -log10(duavee_vs_control$pval)


duavee_vs_control_plot <- ggplot(duavee_vs_control, aes(x=log10foldchange, y=neglog10pval))+
    geom_point(alpha=0.8, aes(color=.data$labeled)) +
    xlab("") + ylab("-Log10 p-value") + theme_bw() +
    theme(legend.position="none", axis.title=element_text(size=10),
          axis.text=element_text(size=10)) +
    scale_x_continuous(limits=c(-2, 2), breaks=seq(-2, 2, 0.4))+
    scale_y_continuous(limits=c(0, 4), breaks=seq(0, 4))+
    scale_color_manual(values=c("#616161", "#ff0066")) +
    geom_vline(xintercept=0, linetype="dashed", color = "blue") +
    geom_label_repel(aes(label = Labels),
                     size=2,
                     max.overlaps = 20,
                     box.padding   = 0.35,
                     point.padding = 0.5,
                     segment.color = 'grey50')


write.csv(duavee_vs_control[, 1:6], "plots/DAA_trt.csv",
          row.names=FALSE, quote=FALSE)
