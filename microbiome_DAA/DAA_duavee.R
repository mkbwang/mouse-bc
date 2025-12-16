

library(ADAPT)
library(phyloseq)
library(dplyr)
library(ggplot2)
library(ggrepel)
library(openxlsx)

rm(list=ls())
mouse_microbiome <- read.csv("data/microbiome/microbiome_data.csv", row.names=1)
mouse_meta <- mouse_microbiome[, seq(1,3)]

mouse_covariate_count <- mouse_meta %>% group_by(ovx_or_intact, Adipo, Rx) %>%
    summarise(Count=n())

mouse_metagenomics <- mouse_microbiome[, -seq(1,3)]
mouse_phyla <- mouse_metagenomics[, 1:6]
mouse_species <- mouse_metagenomics[, 7:118]

physeq_all <- phyloseq(otu_table(as.matrix(mouse_species), taxa_are_rows = FALSE),
                       sample_data(mouse_meta))


adapt_analysis <- function(phyobj, prev_cutoff=0.1){
    count_table <- otu_table(phyobj)@.Data
    min_positive_val <- min(count_table[count_table > 0])
    DAA_output <- adapt(phyobj, cond.var="Rx",
                        base.cond="control",
                        depth.filter=0,
                        censor=min_positive_val,
                        prev.filter=prev_cutoff)
    details_df <- DAA_output@details
    details_df <- details_df %>% arrange(pval)
    return(details_df)
}


viz_helper <- function(DAA_details, p_cutoff=0.05, size=10){

    DAA_details$Labels <- ""
    DAA_details$Labels[DAA_details$pval < p_cutoff] <- DAA_details$Taxa[DAA_details$pval < p_cutoff]
    DAA_details$labeled <- DAA_details$pval < p_cutoff
    DAA_details$neglog10pval <- -log10(DAA_details$pval)

    x_min <- min(DAA_details$log10foldchange) |> floor()
    x_max <- max(DAA_details$log10foldchange) |> ceiling()
    y_max <- ceiling(max(DAA_details$neglog10pval))

    DAA_plot <- ggplot(DAA_details, aes(x=log10foldchange, y=neglog10pval))+
        geom_point(alpha=0.8, aes(color=.data$labeled)) +
        xlab("") + ylab("-Log10 p-value") + theme_bw() +
        scale_x_continuous(limits=c(x_min, x_max), breaks=seq(x_min, x_max, 0.5))+
        scale_y_continuous(limits=c(0, y_max), breaks=seq(0, y_max, 1))+
        theme(legend.position="none", axis.title=element_text(size=size),
              axis.text=element_text(size=size)) +
        scale_color_manual(values=c("#616161", "#ff0066")) +
        geom_vline(xintercept=0, linetype="dashed", color = "blue") +
        geom_label_repel(aes(label = Labels),
                         size=2,
                         max.overlaps = 20,
                         box.padding   = 0.35,
                         point.padding = 0.5,
                         segment.color = 'grey50')

    return(DAA_plot)

}



# lean intact
physeq_lean_intact <- subset_samples(physeq_all, Adipo == "L" &
                                         ovx_or_intact == "intact")
details_lean_intact <- adapt_analysis(phyobj=physeq_lean_intact)
plot_lean_intact <- viz_helper(details_lean_intact, p_cutoff=0.1)



# lean ovx
physeq_lean_ovx <- subset_samples(physeq_all, Adipo == "L" &
                                         ovx_or_intact == "ovx")
details_lean_ovx <- adapt_analysis(phyobj=physeq_lean_ovx)
plot_lean_ovx <- viz_helper(details_lean_ovx, p_cutoff=0.1)


# obese intact
physeq_obese_intact <- subset_samples(physeq_all, Adipo == "Ob" &
                                         ovx_or_intact == "intact")
details_obese_intact <- adapt_analysis(phyobj=physeq_obese_intact)
plot_obese_intact <- viz_helper(details_obese_intact, p_cutoff=0.1)


# obese ovx
physeq_obese_ovx <- subset_samples(physeq_all, Adipo == "Ob" &
                                          ovx_or_intact == "ovx")
details_obese_ovx <- adapt_analysis(phyobj=physeq_obese_ovx)
plot_obese_ovx <- viz_helper(details_obese_ovx, p_cutoff=0.1)



results_list <- list(sample_count=mouse_covariate_count,
                     lean_intact=details_lean_intact,
                     lean_ovx=details_lean_ovx,
                     obese_intact=details_obese_intact,
                     obese_ovx=details_obese_ovx)


wb <- createWorkbook()
for(sheet_name in names(results_list)) {
    addWorksheet(wb, sheet_name)
    writeData(wb, sheet = sheet_name, results_list[[sheet_name]])
}
saveWorkbook(wb, "microbiome_DAA/TrtDAA.xlsx", overwrite = TRUE)


