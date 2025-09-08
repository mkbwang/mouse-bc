
library(dplyr)
library(tidyr)
library(reshape2)
library(patchwork)
library(ggplot2)

rm(list=ls())
mouse_microbiome <- read.csv("microbiome_data.csv", row.names=1)
mouse_meta <- mouse_microbiome[, seq(1,3)]

mouse_metagenomics <- mouse_microbiome[, -seq(1,3)]
mouse_phyla <- mouse_metagenomics[, 1:6]
mouse_species <- mouse_metagenomics[, 7:118]



# sample types
table(mouse_meta$ovx_or_intact)
table(mouse_meta$Adipo)
table(mouse_meta$Rx)

sample_type_count <- mouse_meta %>% group_by(ovx_or_intact, Rx) %>%
    summarise(count=n())


# Visualize relative abundances of phyla
mouse_phyla$Ovary <- mouse_meta$ovx_or_intact
mouse_phyla$Arm <- mouse_meta$Rx
mouse_phyla<- mouse_phyla %>% select(Ovary, Arm, everything())
write.csv(mouse_phyla, "data/mouse_phyla.csv", row.names=FALSE)
mouse_phyla_long <- reshape2::melt(mouse_phyla, id.vars=c("Arm", "Ovary"),
                         variable.name="Phyla")
colnames(mouse_phyla_long)[4] <- "relabd"
mouse_phyla_long$relabd <- mouse_phyla_long$relabd*100
mouse_phyla_long$Type <- sprintf("%s_%s",
                                 mouse_phyla_long$Ovary,
                                 mouse_phyla_long$Arm)

phyla_plot <- ggplot(mouse_phyla_long, aes(x=Phyla, y=relabd, color=Type)) +
    geom_boxplot(position = position_dodge(width = 0.8), width = 0.7, outlier.shape = NA) +
    geom_jitter(
        aes(color = Type),
        position = position_jitterdodge(jitter.width = 0.15, dodge.width = 0.8),
        size = 1, alpha = 0.7
    ) +
    theme_bw() +
    scale_y_continuous(breaks=seq(0, 100, 20), limits=c(0, 100))+
    theme(axis.text.x = element_text(angle = 45, hjust = 1))+
    labs(x = "Phyla", y = "Relative Abundance (%)") +
    scale_color_manual(values=c("#45A2D8", "#B3AE00", "#333333", "#ef3F34"))





# colored bar for relative abundance of species

## get prevalences
prevalences <- colMeans(mouse_species > 0)
prev_taxa <- names(which(prevalences >= 0.9)) # prevalent taxa with occurrence in at least 90% of all the samples
mean_relabd <- colMeans(mouse_species[, prev_taxa])
mean_relabd <- sort(mean_relabd, decreasing = T)
sorted_taxa <- names(mean_relabd) # sort by mean relative abundance
## merge other taxa
mouse_species_others <- mouse_species %>% select(-any_of(prev_taxa)) |> rowSums()

mouse_species_concise <- cbind(mouse_species[, sorted_taxa],
                               mouse_species_others)
colnames(mouse_species_concise)[ncol(mouse_species_concise)] <- "Others"
mouse_species_concise$Ovary <- mouse_meta$ovx_or_intact
mouse_species_concise$Arm <- mouse_meta$Rx

## take mean of relative abundance grouped by treatment arm and ovary
mouse_species_concise_mean <- mouse_species_concise %>% group_by(Arm, Ovary) %>%
    summarise_all(mean)
write.csv(mouse_species_concise_mean, "data/mouse_species_mean.csv",
          row.names=FALSE, quote=FALSE)
mouse_species_concise_long <- reshape2::melt(mouse_species_concise_mean, id.vars=c("Arm", "Ovary"),
                                   variable.name="Taxa")


mouse_species_concise_long$Group <- sprintf("%s_%s", mouse_species_concise_long$Ovary,
                                            mouse_species_concise_long$Arm)



mouse_species_concise_long$value <- mouse_species_concise_long$value*100
species_relabd_plot <- ggplot(mouse_species_concise_long, aes(x = Group, y = value, fill = Taxa)) +
    geom_bar(stat = "identity") +
    ylab("Relative Abundance (%)") +
    scale_y_continuous(breaks=seq(0, 100, 20), limits=c(0, 100.1))+
    scale_fill_manual(values=c("#E69F00","#56B4E9",
                               "#009E73","#F0E442",
                               "#0072B2", "#D55E00",
                               "#CC79A7", "#999999",
                               "#1B9E77", "#D95F02",
                               "#7570B3", "#E7298A",
                               "#66A61E", "#E6AB02",
                               "#A6761D", "#666666",
                               "#FB9A99", "#B2DF8A",
                               "#CAB2D6", "#000000"))+
    theme_classic()+
    theme(axis.text.x = element_text(angle = 45, hjust = 1))+
    guides(color = guide_legend(ncol = 2),
           fill = guide_legend(ncol = 2))



# Pcoa plot based on bray curtis distance

library(vegan)


## species
pcoa_analysis <- function(mydata){

    bcdist <- vegdist(mydata, method = "bray")
    pcoa <- cmdscale(bcdist, k=2, eig=TRUE)
    coords <- pcoa$points |> as.data.frame()

    prop_var_1 <- pcoa$eig[1] / sum(pcoa$eig) * 100
    prop_var_2 <- pcoa$eig[2] / sum(pcoa$eig) * 100
    title_x <- sprintf("PCOA1 (%.2f%%)", prop_var_1)
    title_y <- sprintf("PCOA2 (%.2f%%)", prop_var_2)

    return(list(coords=coords, title_x=title_x, title_y=title_y))

}

mouse_species <- mouse_metagenomics[, 7:118]

pcoa_results <- pcoa_analysis(mouse_species)
pcoa_df <- pcoa_results$coords
pcoa_df$Type <- sprintf("%s_%s",
                        mouse_meta$ovx_or_intact,
                        mouse_meta$Rx)

pcoa_plot_treatment_intact <- ggplot(pcoa_df, aes(x=V1, y=V2, color=Type)) +
    geom_point(alpha=0.8, size=1.5) +
    xlab(pcoa_results$title_x) +
    ylab(pcoa_results$title_y) +
    scale_color_manual(values=c("#45A2D8", "#B3AE00", "#333333", "#ef3F34")) +
    theme_bw(base_size=10)+
    ggtitle("MDS Plot based on Bray-curtis distance")


permanova_result <- adonis2(mouse_species ~ mouse_meta$Rx+mouse_meta$ovx_or_intact+mouse_meta$Adipo,
                                      method = "bray", permutations = 999)

