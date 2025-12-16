

library(ADAPT)
library(phyloseq)
library(dplyr)
library(ggplot2)
library(ggrepel)
library(openxlsx)

rm(list=ls())
mouse_microbiome <- read.csv("data/microbiome/microbiome_data.csv", row.names=1)
mouse_meta <- mouse_microbiome[, seq(1,3)]

mouse_covariate_count <- mouse_meta %>% group_by(Adipo, Rx) %>%
    summarise(Count=n())

mouse_metagenomics <- mouse_microbiome[, -seq(1,3)]
mouse_phyla <- mouse_metagenomics[, 1:6]
mouse_species <- mouse_metagenomics[, 7:118]

physeq_all <- phyloseq(otu_table(as.matrix(mouse_species), taxa_are_rows = FALSE),
                       sample_data(mouse_meta))


adapt_analysis <- function(phyobj, prev_cutoff=0.1){
    count_table <- otu_table(phyobj)@.Data
    min_positive_val <- min(count_table[count_table > 0])
    DAA_output <- adapt(phyobj, cond.var="ovx_or_intact",
                        base.cond="intact",
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



# lean control
physeq_lean_control <- subset_samples(physeq_all, Adipo == "L" & Rx == "control")
details_lean_control <- adapt_analysis(phyobj=physeq_lean_control)
plot_lean_control <- viz_helper(details_lean_control, p_cutoff=0.1)


# lean duavee
physeq_lean_duavee <- subset_samples(physeq_all, Adipo == "L" & Rx == "duavee")
details_lean_duavee <- adapt_analysis(phyobj=physeq_lean_duavee)
plot_lean_duavee <- viz_helper(details_lean_duavee, p_cutoff=0.1)


# obese control
physeq_obese_control <- subset_samples(physeq_all, Adipo == "Ob" & Rx == "control")
details_obese_control <- adapt_analysis(phyobj=physeq_obese_control)
plot_obese_control <- viz_helper(details_obese_control, p_cutoff=0.1)



# obese duavee
physeq_obese_duavee <- subset_samples(physeq_all, Adipo == "Ob" & Rx == "duavee")
details_obese_duavee <- adapt_analysis(phyobj=physeq_obese_duavee)
plot_obese_duavee <- viz_helper(details_obese_duavee, p_cutoff=0.1)



results_list <- list(lean_control=details_lean_control,
                     lean_duavee=details_lean_duavee,
                     obese_control=details_obese_control,
                     obese_duavee=details_obese_duavee)


wb <- createWorkbook()
for(sheet_name in names(results_list)) {
    addWorksheet(wb, sheet_name)
    writeData(wb, sheet = sheet_name, results_list[[sheet_name]])
}
saveWorkbook(wb, "microbiome_DAA/ovaryDAA.xlsx", overwrite = TRUE)


