#!/usr/bin/env Rscript
.libPaths(Sys.getenv("R_LIBS_USER"))

if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("Biostrings")
BiocManager::install("dplyr")
BiocManager::install("argparse")
library(Biostrings)
library(dplyr)

read_data_and_labels <- function(data_file, labels_file) {
  data <- read.delim(data_file, header = TRUE, sep = ";")
  labels <- read.delim(labels_file)
  selected_labels <- labels$Tag
  data <- data[data$Tag %in% selected_labels, ] # Select rows with specified labels
  return(data)
}


fasta_file <- snakemake@input[['fasta_file']]
labels_file <- snakemake@input[['labels_file']]


data_file <- "data/all_tag_by_cluster_counts.txt"
 

data <- read_data_and_labels(data_file, labels_file)

dataunfiltered <- data

# Read FASTA file
read_fasta_file <- function(fasta_file) {
  sequences <- readDNAStringSet(fasta_file)
  new_headers <- sub("_repseq_0", "", names(sequences))
  names(sequences) <- new_headers
  return(sequences)
}

sequences <- read_fasta_file(fasta_file) 

# Remove columns not present in FASTA file
remove_columns_not_in_fasta <- function(data, sequences) {
  columns_to_remove <- setdiff(names(data), names(sequences))
  data <- data[, !names(data) %in% columns_to_remove]
  return(data)
}

data <- remove_columns_not_in_fasta(data, sequences)

#proportion of data that was non-fungal
#reads_lost <- (sum(dataunfiltered[,-1])-sum(data))/sum(dataunfiltered[,-1])*100

#Calculate relative abundances
calculate_rel_abundance <- function(data) {
  sample_totals <- rowSums(data)
  data <- data / sample_totals
  return(data)
}

data <- calculate_rel_abundance(data)

# Function to process a row 
process_row <- function(row) {
  sorted_row <- sort(row, decreasing = TRUE)  # Sort the row in descending order
  cumsum_row <- cumsum(sorted_row)  # Calculate the cumulative sum
  index <- which.max(cumsum_row >= 0.80)  # Find the index where the cumulative sum first reaches 90% or more
  selected_columns <- names(sorted_row)[1:index]   # Select columns up to the identified index
  return(selected_columns)
}

# Process rows and write matching sequences to a new FASTA file
process_rows_and_write_fasta <- function(data, sequences, fasta_file) {
  result <- apply(data[, -1], 1, process_row) #Apply the process_row function to each row in the data frame
  unique_values <- unique(unlist(result)) # Combine and filter unique values
  matching_sequences <- sequences[names(sequences) %in% unique_values]
  writeXStringSet(matching_sequences, file = fasta_file) # Write matching sequences to a new FASTA file
}

write.csv(data, file = "processed_data.csv", row.names = FALSE)
finalfasta <- "80relabundITS.fasta"
process_rows_and_write_fasta(data, sequences, finalfasta)
