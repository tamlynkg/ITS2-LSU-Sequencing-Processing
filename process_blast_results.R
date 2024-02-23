#!/usr/bin/env Rscript
.libPaths(Sys.getenv("R_LIBS_USER"))
library(Biostrings)
library(argparse)
library(dplyr)

##### Process BLAST results
update_blast_results <- function(blast_input) {
  blast_df <- read.table(blast_input, header = FALSE, sep = "\t")
  # Rename the columns for clarity
  colnames(blast_df) <- c("qseqid", "sseqid", "slen", "pident", "length", "mismatch", "gapopen", "qlen", "qstart", "qend", "sstart", "send", "evalue", "bitscore")
  return(blast_df)
}

blast_input <- "blast_results.txt"
blast_file <- update_blast_results(blast_input)

#duplicates <- blast_file %>% 
#  group_by(sseqid) %>% 
#  filter(n() > 1) 

data <- read.csv("processed_data.csv", sep = ",", header = TRUE)

#Remove non-fungal 
if (any(grepl("NF", blast_file$sseqid))) {
    # Subset the blast_file data frame to select rows where "sseqid" contains "NF"
    subset_df <- subset(blast_file, grepl("NF", sseqid))
    # Extract qseqid values from the subsetted data frame
    qseqid_values <- subset_df$qseqid
    # Remove columns from the data data frame where column names match qseqid_values
    data <- data[, !(names(data) %in% qseqid_values)]
}


process_blast_results <- function(blast_file) {
  blast_file$qlen_98_percent <- 0.98 * blast_file$qlen   # Calculate the 98% of slen
  blast_file$sseqid[blast_file$length < blast_file$qlen_98_percent | blast_file$pident < 98] <- "UNKNOWN"
  scatadf <- blast_file[, 1:2]
  return(scatadf)
}

scatadf <- process_blast_results(blast_file)

if ("Lichenized" %in% scatadf$sseqid) {
  scatadf <- scatadf[!grepl("Lichenized", scatadf$sseqid, ignore.case = TRUE), ]
}

read_fasta_file <- function(fasta_file) {
  sequences <- readDNAStringSet(fasta_file)
  new_headers <- sub("_repseq_0", "", names(sequences))
  names(sequences) <- new_headers
  return(sequences)
}

sequences <- read_fasta_file(snakemake@input[['fasta_file']])

unknowns <- scatadf[grepl("UNKNOWN", scatadf$sseqid, ignore.case = FALSE), ]
noblast <- names(sequences)[!names(sequences) %in% scatadf$qseqid]

sequences_unknown <- union(sequences[names(sequences) %in% unknowns], sequences[names(sequences) %in% noblast])

writeXStringSet(sequencesunknown, file="Unknownsequences.fasta")

write.csv(scatadf, file = "processed_SIblast_results.csv", row.names = FALSE)

