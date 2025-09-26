
rm(list=ls())
library(oligo)
library(clariomsrattranscriptcluster.db)
library(AnnotationDbi)
library(openxlsx)
library(dplyr)
library(pd.clariom.s.rat)
library(readr)

files <- list.files("data/microarray_cel", pattern="*mg.CEL")
files_fullname <- sprintf("data/microarray_cel/%s", files)

# load metadata
microarray_metadata <- read.xlsx("data/giles_duavee_microarray_metadata.xlsx",
                                 sheet=1)
microarray_metadata <- microarray_metadata[microarray_metadata$Tissue == "mg", ]
rownames(microarray_metadata) <- microarray_metadata$Sample
microarray_metadata <- microarray_metadata[, seq(1, 8)]


mouse_ovx <- read.celfiles(filenames=files_fullname)
sns <- sampleNames(mouse_ovx)
sns <- substr(sns, 1, 3)
sampleNames(mouse_ovx) <- sns


# normalize with rma
eset <- rma(mouse_ovx)
abundance <- exprs(eset)
featureIDs <- featureNames(eset)


# find annotations of sequences
full_annotations <- read_tsv("data/GPL23040-70516.txt",
                             comment="#")
reference_df <- sapply(full_annotations$mrna_assignment, function(longinfo){
    info_vec <- strsplit(longinfo, split=" // ")[[1]]
    info_vec[c(1,2)]
}) |> unname() |> t()
full_annotations$DB <- reference_df[, 2]
full_annotations$number <- reference_df[, 1]

# retain features which have annotations in known databases
subset_annotations <- full_annotations %>%
    filter(DB %in% c("ENSEMBL", "RefSeq", "Rat Genome Database", "Mm2Rn5")) %>%
    select(ID, DB, number, seqname, strand, start, stop, total_probes) |> as.data.frame()
rownames(subset_annotations) <- subset_annotations$ID


abundance_subset <- abundance[rownames(subset_annotations), ]
metadata_subset <- microarray_metadata[colnames(abundance), ]

library(SummarizedExperiment)

microarray_ovx <- SummarizedExperiment(
    assays=list(intensity=abundance_subset),
    colData=metadata_subset,
    rowData=subset_annotations
)

saveRDS(microarray_ovx, "data/processed_microarray.rds")








