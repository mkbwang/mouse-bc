
library(dplyr)
library(tidyr)
library(reshape2)
library(patchwork)

rm(list=ls())
mouse_microbiome <- read.csv("microbiome_data.csv", row.names=1)
mouse_meta <- mouse_microbiome[, seq(1,3)]

mouse_metagenomics <- mouse_microbiome[, -seq(1,3)]
mouse_phyla <- mouse_metagenomics[, 1:6]
mouse_species <- mouse_metagenomics[, 7:118]


# repeat analysis

# sample types
table(mouse_meta$ovx_or_intact)
table(mouse_meta$Adipo)
table(mouse_meta$Rx)

sample_type_count <- mouse_meta %>% group_by(ovx_or_intact, Rx) %>%
    summarise(count=n())


# Visualize relative abundances of phyla
mouse_phyla$Ovary <- mouse_meta$ovx_or_intact
mouse_phyla$Arm <- mouse_meta$Rx
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



##TODO: differential abundance analysis



# colored bar for relative abundance of species

## get prevalences
prevalences <- colMeans(mouse_species > 0)
prev_taxa <- names(which(prevalences >= 0.9)) # prevalent taxa with occurrence in at least 90% of all the samples
mean_relabd <- colMeans(mouse_species[, prev_taxa]) |> sort(decreasing=T)
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



## compare duavee treatment vs control for ovary-intact mice
mouse_species_intact <- mouse_species[mouse_meta$ovx_or_intact == "intact", ]
mouse_meta_intact <- mouse_meta %>% filter(ovx_or_intact == "intact")
species_intact_pcoa <- pcoa_analysis(mouse_species_intact)
intact_pcoa_df <- as.data.frame(species_intact_pcoa$coords)
intact_pcoa_df$Treatment <- mouse_meta_intact$Rx
permanova_treatment_intact <- adonis2(mouse_species_intact ~ mouse_meta_intact$Rx,
                                      method = "bray", permutations = 999)
pval <- permanova_treatment_intact$`Pr(>F)`[1]
title_name <- sprintf("Duavee vs Control for Ovary-Intact Mice (p value %.2f)", pval)
pcoa_plot_treatment_intact <- ggplot(intact_pcoa_df, aes(x=V1, y=V2, color=Treatment)) +
    geom_point(alpha=0.8, size=1.5) +
    xlab(species_intact_pcoa$title_x) +
    ylab(species_intact_pcoa$title_y) +
    scale_color_manual(values=c("#45A2D8", "#B3AE00")) +
    theme_bw(base_size=8)+
    ggtitle(title_name)


## compare duavee treatment vs control for ovx mice
mouse_species_ovx <- mouse_species[mouse_meta$ovx_or_intact == "ovx", ]
mouse_meta_ovx <- mouse_meta %>% filter(ovx_or_intact == "ovx")
species_ovx_pcoa <- pcoa_analysis(mouse_species_ovx)
ovx_pcoa_df <- as.data.frame(species_ovx_pcoa$coords)
ovx_pcoa_df$Treatment <- mouse_meta_ovx$Rx
permanova_treatment_ovx <- adonis2(mouse_species_ovx ~ mouse_meta_ovx$Rx,
                                      method = "bray", permutations = 999)
pval <- permanova_treatment_ovx$`Pr(>F)`[1]
title_name <- sprintf("Duavee vs Control for OVX Mice (p value %.2f)", pval)
pcoa_plot_treatment_ovx <- ggplot(ovx_pcoa_df, aes(x=V1, y=V2, color=Treatment)) +
    geom_point(alpha=0.8, size=1.5) +
    xlab(species_ovx_pcoa$title_x) +
    ylab(species_ovx_pcoa$title_y) +
    scale_color_manual(values=c("#333333", "#ef3F34")) +
    theme_bw(base_size=8)+
    ggtitle(title_name)


## Compare ovx over intact for the control group
mouse_species_control <- mouse_species[mouse_meta$Rx == "control", ]
mouse_meta_control <- mouse_meta %>% filter(Rx == "control")
species_control_pcoa <- pcoa_analysis(mouse_species_control)
control_pcoa_df <- as.data.frame(species_control_pcoa$coords)
control_pcoa_df$Ovary <- mouse_meta_control$ovx_or_intact
permanova_ovary_control <- adonis2(mouse_species_control ~ mouse_meta_control$ovx_or_intact,
                                   method = "bray", permutations = 999)
pval <- permanova_ovary_control$`Pr(>F)`[1]
title_name <- sprintf("OVX vs Intact for Control Mice (p value %.2f)", pval)
pcoa_plot_ovary_control <- ggplot(control_pcoa_df, aes(x=V1, y=V2, color=Ovary)) +
    geom_point(alpha=0.8, size=1.5) +
    xlab(species_control_pcoa$title_x) +
    ylab(species_control_pcoa$title_y) +
    scale_color_manual(values=c("#45A2D8", "#333333")) +
    theme_bw(base_size=8)+
    ggtitle(title_name)


## compare ovx over intact for the duavee group
mouse_species_duavee <- mouse_species[mouse_meta$Rx == "duavee", ]
mouse_meta_duavee <- mouse_meta %>% filter(Rx == "duavee")
species_duavee_pcoa <- pcoa_analysis(mouse_species_duavee)
duavee_pcoa_df <- as.data.frame(species_duavee_pcoa$coords)
duavee_pcoa_df$Ovary <- mouse_meta_duavee$ovx_or_intact
permanova_ovary_duavee <- adonis2(mouse_species_duavee ~ mouse_meta_duavee$ovx_or_intact,
                                   method = "bray", permutations = 999)
pval <- permanova_ovary_duavee$`Pr(>F)`[1]
title_name <- sprintf("OVX vs Intact for Duavee Mice (p value %.2f)", pval)
pcoa_plot_ovary_duavee <- ggplot(duavee_pcoa_df, aes(x=V1, y=V2, color=Ovary)) +
    geom_point(alpha=0.8, size=1.5) +
    xlab(species_duavee_pcoa$title_x) +
    ylab(species_duavee_pcoa$title_y) +
    scale_color_manual(values=c("#B3AE00", "#ef3F34")) +
    theme_bw(base_size=8)+
    ggtitle(title_name)


pcoa_plot_species <- wrap_plots(pcoa_plot_treatment_intact,
                                pcoa_plot_treatment_ovx,
                                pcoa_plot_ovary_control,
                                pcoa_plot_ovary_duavee,
                                ncol=2)

