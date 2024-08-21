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
  new_headers <- sub("_repseq_0", "", names(sequences))
  names(sequences) <- new_headers
  return(sequences)
}

sequences <- read_fasta_file("Atnarova.LSU.fasta")

columns_to_remove <- setdiff(names(sequences), blastdf$qseqid)
exact_sequences <- intersect(names(sequences), blastdf$qseqid)

matching_sequences <- sequences[!names(sequences) %in% columns_to_remove]

output_fasta_file <- "sequences-LSU-atnarova.fasta"
writeXStringSet(matching_sequences, file=output_fasta_file)
