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

blastdf <- read.csv("bestmatches.csv", header = TRUE)

read_fasta_file <- function(fasta_file) {
  sequences <- readDNAStringSet(fasta_file)
 # new_headers <- sub("_repseq_0", "", names(sequences))
 # names(sequences) <- new_headers
  return(sequences)
}

fasta_data <- readLines("sequences-LSU-atnarova.fasta")

# Initialize variables to store the modified FASTA sequences
modified_fasta <-  character(0)

# Loop through the FASTA data and update the headers
i <- 1

while (i <= length(fasta_data)) {
  if (startsWith(fasta_data[i], ">")) {
    # Extract the header from the current line
    fasta_header <- substring(fasta_data[i], 2)  # Remove the ">" character

    # Search for a matching qseqid in the blast_df
    matching_row <- blastdf[blastdf$qseqid == fasta_header, ]

    if (nrow(matching_row) > 0) {
      # Extract qseqid and sseqid from the matching row in blast_df
      qseqid <- matching_row$qseqid
      species <- matching_row$species
      # Modify the header to include qseqid, sseqid, and Cluster.Size
      fasta_data[i] <- paste0(">", qseqid, "_", species)

    }
  }
  # Add the line to the modified FASTA data
  modified_fasta <- c(modified_fasta, fasta_data[i])
  i <- i + 1
}

# Write the modified FASTA data to a new file
fasta_output <- "annotated_LSU_full.fasta"
writeLines(modified_fasta, fasta_output)
