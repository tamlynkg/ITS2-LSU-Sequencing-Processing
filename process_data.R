#!/usr/bin/env Rscript
.libPaths(Sys.getenv("R_LIBS_USER"))

if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("Biostrings")
BiocManager::install("dplyr")
BiocManager::install("argparse")
BiocManager::install("tidyr")
library(Biostrings)
library(dplyr)
library(tidyr)

# Read FASTA file
read_fasta_file <- function(fasta_file) {
  sequences <- readDNAStringSet(fasta_file)
#  new_headers <- sub("_repseq_0", "", names(sequences))
#  names(sequences) <- new_headers
  return(sequences)
}

fasta_file <- snakemake@input[['fasta_file']]

sequences <- read_fasta_file(fasta_file) 

#Getting ITS Length 

#ITSx_positions <- snakemake@input[['ITSx_positions]]

df <- read.delim("Atnarova.positions.txt", header = FALSE)

colnames(df) <- c("scata_ref", "length", "SSU", "ITS1", "5.8S", "ITS2", "LSU", "sequence")

#so ITSx positions only includes fungal
df$scata_ref <- gsub("_repseq_0", "", df$scata_ref)
df <- df[df$scata_ref %in% names(sequences), ]

#Remove reads where a LSU is not found
df <- df %>%
  filter(LSU != "LSU: Not found")

read_data_and_labels <- function(data_file, labels_file) {
  data <- read.delim(data_file, header = TRUE, sep = ";")
  labels <- read.delim(labels_file)
  selected_labels <- labels$Tag
  data <- data[data$Tag %in% selected_labels, ] # Select rows with specified labels
  return(data)
}


labels_file <- snakemake@input[['labels_file']]

data_file <- "data/all_tag_by_cluster_counts.txt"
 
data <- read_data_and_labels(data_file, labels_file)

dataunfiltered <- data

labels <- read.delim(labels_file)

filtered.counts.data <- data[,-1]
filtered.counts.data <- filtered.counts.data[, colSums(filtered.counts.data) != 0]
filtered.counts.data$Sample <- labels$Sample
filtered.counts.data$Sample <- trimws(filtered.counts.data$Sample)


# Remove columns not present in FASTA file as they are non-fungal
remove_columns_not_in_fasta <- function(data, sequences) {
  columns_to_remove <- setdiff(names(data), names(sequences))
  data <- data[, !names(data) %in% columns_to_remove]
  data <- data[, names(data) %in% df$scata_ref]
  return(data)
}

filtered.counts.data <- remove_columns_not_in_fasta(filtered.counts.data, sequences)

matching_sequences <- sequences[names(sequences) %in% names(filtered.counts.data)]

# Write matching sequences to a new FASTA file
output_fasta_file <- "source_11.fasta"
writeXStringSet(matching_sequences, file=output_fasta_file)

filtered.counts.data$Sample <- labels$Sample
filtered.counts.data$Sample <- substr(filtered.counts.data$Sample, 1, 3)
filtered.counts.data <- aggregate(. ~ Sample, filtered.counts.data, FUN = sum)

write.csv(filtered.counts.data, "countsdata_filtered.csv")




