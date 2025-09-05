# this file processes Erin's mouse microbiome data  (GL, 8/25/2025)
setwd('C:\\Users\\ligen\\Dropbox (Personal)\\Research File\\Collaborative Projects\\WithErinGiles\\Microbiome project')
library(phyloseq)
library(microbiome)
library(ggplot2)
library(plyr)
library(vegan)
# if (!require("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# 
# BiocManager::install("microbiome")



####################################### PART I: Load and Prepare ###########################################


##### load data  ###############
alldata=read.table(file='microbiome_data.csv',header = TRUE,sep=',')
dim(alldata) # 67*122, each row is a sample (67 samples, 118 taxa)
unique(alldata$Rat_ID)
metadata=alldata[,1:4]
rownames(metadata)=metadata[,1]
OTUdata=alldata[,-c(1:4)]
rownames(OTUdata)=metadata[,1]



# ######### Load taxonomy data ####################
# taxonomy=read.table(file='Prospective Relative Abundance 022122 taxonomy.csv',header = TRUE,sep=',')
# dim(taxonomy) # 13341*9
# colnames(taxonomy)
# sum(is.na(taxonomy)) # no missing
# taxonomy[taxonomy=='0']=NA
# rownames(taxonomy)=taxonomy$Feature.ID
# taxonomy=taxonomy[,2:8]
# for(level in 1:7){
#   print(paste('#missing/#unique in',colnames(taxonomy)[level],':',sum(is.na(taxonomy[,level])),'/',length(unique(taxonomy[,level]))))
# } # note: many low-level (e.g., species) taxa are grouped together into NA or unknown or unidentified
# dim(taxonomy) # from high to low levels
# 
# # note: taxonomy is a matrix (of tree) with NA; if needed, it can be converted to a same-depth tree with redundant nodes (no NA)






# convert to phyloseq data type
OTU = otu_table(OTUdata, taxa_are_rows = FALSE) # n=67 p=118 
## NOTE: it seems the first 6 taxa almost always add up to 1, and the remaining 112 taxa add up to 1 (waiting for Erin's clarification)
OTU1 = otu_table(OTUdata[,1:6], taxa_are_rows = FALSE)
OTU2 = otu_table(OTUdata[,-c(1:6)], taxa_are_rows = FALSE)
TAX = NA #tax_table(as.matrix(phase1.tax))
samples = sample_data(metadata)

phylodata <- phyloseq(OTU, TAX, samples)
phylodata1 <- phyloseq(OTU1, TAX, samples)
phylodata2 <- phyloseq(OTU2, TAX, samples)

save(OTU,samples,phylodata,phylodata1,phylodata2, file="processed_data.Rdata")








####################################### PART II: PCoA and PERMANOVA #######################################
# plot composition
plot_composition(phylodata1, otu.sort = "abundance", group_by='ovx_or_intact', x.label='sample', plot.type = 'barplot')
plot_composition(phylodata2, otu.sort = "abundance", group_by='ovx_or_intact', x.label='sample', plot.type = 'barplot')


# ordination plot (ovx vs intact, lean vs obese, duavee vs control)
phylo1.ord <- ordinate(phylodata1, "PCoA", "bray")
colors=samples$Adipo
#colors=samples$Rx


plot_ordination(phylodata1, phylo1.ord, type="samples") +
  geom_point(size=2,color='white',shape=16)+ # mask up default points from ordination plot
  geom_point(aes(color=colors, shape=ovx_or_intact),size=2) + 
  scale_shape_manual(values = c(17,19)) + 
  ggtitle(paste("Ordination Plot (Bray-Curtis)"))
#ggsave(paste0('PCoA.pdf'),width=8, height=8) 




### permanova between groups
# PERMANOVA
adonis2(phyloseq::distance(phylodata1, method="bray") ~ sample_data(phylodata1)$Adipo)   # pooled data; p=0.102
adonis2(phyloseq::distance(phylodata1, method="bray") ~ sample_data(phylodata1)$Rx)   # pooled data; p=0.44



## compare alpha diversity