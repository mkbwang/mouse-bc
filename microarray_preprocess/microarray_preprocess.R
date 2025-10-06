
rm(list=ls())
library(oligo)
library(clariomsrattranscriptcluster.db)
library(AnnotationDbi)
library(openxlsx)
library(dplyr)
library(pd.clariom.s.rat)
library(readr)

files <- list.files("data/microarray/cel", pattern="*mg.CEL")
files_fullname <- sprintf("data/microarray/cel/%s", files)

# load metadata
microarray_metadata <- read.xlsx("data/microarray/giles_duavee_microarray_metadata.xlsx",
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
full_annotations <- read_tsv("data/microarray/GPL23040-70516.txt",
                             comment="#")
reference_df <- sapply(full_annotations$mrna_assignment, function(longinfo){
    info_vec <- strsplit(longinfo, split=" // ")[[1]]
    info_vec[c(1,2)]
}) |> unname() |> t()
full_annotations$DB <- reference_df[, 2]
full_annotations$accession <- reference_df[, 1]

# retain features which have annotations in known databases
subset_annotations <- full_annotations %>%
    filter(DB %in% c("ENSEMBL", "RefSeq", "Mm2Rn5")) %>%
    select(ID, DB, accession, seqname, strand, start, stop, total_probes) |> as.data.frame()
rownames(subset_annotations) <- subset_annotations$ID


abundance_subset <- abundance[rownames(subset_annotations), ]
metadata_subset <- microarray_metadata[colnames(abundance), ]

library(SummarizedExperiment)

microarray_ovx <- SummarizedExperiment(
    assays=list(intensity=abundance_subset),
    colData=metadata_subset,
    rowData=subset_annotations
)

saveRDS(microarray_ovx, "data/microarray/processed_microarray.rds")


# after annotating the accessions, run the script below
microarray_ovx <- readRDS("data/microarray/processed_microarray.rds")
feature_metadata <- rowData(microarray_ovx) |> as.data.frame()
clean_accession <- function(accession){
  
  split_names <- strsplit(accession, split="_")[[1]]
  if(length(split_names) > 2){ # remove redundant IDs
    newname <- paste(split_names[c(1,2)], collapse="_")
  } else{
    newname <- accession
  }
  return(newname)
}

feature_metadata$accession <- sapply(feature_metadata$accession,
                                     clean_accession)

# load annotations
probe_ensembl_annotation <- read.csv("data/microarray/probe_ensembl_annotation.txt")
probe_ensembl_annotation$Gene <- probe_ensembl_annotation$Genename
probe_ensembl_annotation$Gene[probe_ensembl_annotation$Genename %in% c("Error Fetching", "No gene symbol")] <- "No abbreviation"


probe_refseq_annotation <- read.csv("data/microarray/probe_refseq_annotation.txt")
library(stringr)
extract_genename <- function(description){
  matches <- str_match_all(description, "\\((.*?)\\)")
  if (length(matches) >0 && ncol(matches[[1]]) > 1){
    info <- matches[[1]][, 2]
    if (length(info) != 0){
      return(info[length(info)])
    } else{
      return("No abbreviation")
    }
  } else{
    return("No abbreviation")
  }
}

refseq_genenames <- sapply(probe_refseq_annotation$Genename, extract_genename) |> unlist()
probe_refseq_annotation$Gene <- refseq_genenames
probe_annotation_combined <- rbind(probe_ensembl_annotation,
                                   probe_refseq_annotation) %>% select(Probe, Gene) %>%
  rename(accession=Probe) |> unique()

feature_metadata <- feature_metadata %>% 
  left_join(probe_annotation_combined, by="accession")

rownames(feature_metadata) <- feature_metadata$ID
rowData(microarray_ovx) <- feature_metadata


saveRDS(microarray_ovx, "data/microarray/processed_microarray.rds")


