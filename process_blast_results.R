#!/usr/bin/env Rscript
.libPaths(Sys.getenv("R_LIBS_USER"))
library(Biostrings)
library(argparse)
library(dplyr)
library(data.table)
library(stringr)

##### Process BLAST results
update_blast_results <- function(blast_input) {
  blast_df <- fread(blast_input, header = FALSE, sep = "\t")
  # Rename the columns for clarity
  colnames(blast_df) <- c("qseqid", "sseqid", "slen", "pident", "length", "mismatch", "gapopen", "qlen", "qstart", "qend", "sstart", "send", "evalue", "bitscore")
  return(blast_df)
}

blast_input <- "blast_source11_newdb10.txt"
blast_file <- update_blast_results(blast_input)

blast_file98 <- blast_file %>%
  filter(pident >= 98)

filtered_df <- blast_file98[grepl("refs", blast_file98$sseqid), ]

setorder(filtered_df, evalue, -bitscore, -pident)
best_matches_newdb_refs <- filtered_df[, .SD[1], by = qseqid]

speciesvalue_refs <- str_extract(best_matches_newdb_refs$sseqid, "s__([^;]+)")
speciesvalueSH_refs <- gsub("^s__", "", speciesvalue_refs)
best_matches_newdb_refs$species <- speciesvalueSH_refs

phylumvalue_refs <- str_extract(best_matches_newdb_refs$sseqid, "p__([^;]+)")
phylumvalueSH_refs <- gsub("^p__", "", phylumvalue_refs)
best_matches_newdb_refs$phylum <- phylumvalueSH_refs

ordervalue_refs <- str_extract(best_matches_newdb_refs$sseqid, "o__([^;]+)")
ordervalueSH_refs <- gsub("^o__", "", ordervalue_refs)
best_matches_newdb_refs$order <- ordervalueSH_refs


filtered_df1 <- blast_file[!blast_file$qseqid %in% best_matches_newdb_refs$qseqid, ]

filtered_df1 <- filtered_df1 %>%
  filter(pident >= 95)

setorder(filtered_df1, evalue, -bitscore, -pident)
best_matches_newdb_reps <- filtered_df1[, .SD[1], by = qseqid]

speciesvalue <- str_extract(best_matches_newdb_reps$sseqid, "s__([^;]+)")
speciesvalueSH <- gsub("^s__", "", speciesvalue)
best_matches_newdb_reps$species <- speciesvalueSH
best_matches_newdb_reps$species <- sub("_.*", "", best_matches_newdb_reps$species)
best_matches_newdb_reps$species <- paste0(best_matches_newdb_reps$species, "_sp_REPS")

phylumvalue <- str_extract(best_matches_newdb_reps$sseqid, "p__([^;]+)")
phylumvalueSH <- gsub("^p__", "", phylumvalue)
best_matches_newdb_reps$phylum <- phylumvalueSH

ordervalue <- str_extract(best_matches_newdb_reps$sseqid, "o__([^;]+)")
ordervalueSH <- gsub("^o__", "", ordervalue)
best_matches_newdb_reps$order <- ordervalueSH


best_matches_newdb <- rbind(best_matches_newdb_refs, best_matches_newdb_reps)

write.csv(best_matches_newdb, "bestmatches.csv")

read_fasta_file <- function(fasta_file) {
  sequences <- readDNAStringSet(fasta_file)
  return(sequences)
}

sequences <- read_fasta_file("source_11.fasta")

noblast <- names(sequences)[!names(sequences) %in% best_matches_newdb$qseqid]

fasta_data <- readLines("source_11.fasta")

# Initialize variables to store the modified FASTA sequences
modified_fasta <-  character(0)

# Loop through the FASTA data and update the headers
i <- 1

while (i <= length(fasta_data)) {
  if (startsWith(fasta_data[i], ">")) {
    # Extract the header from the current line
    fasta_header <- substring(fasta_data[i], 2)  # Remove the ">" character
    
    # Search for a matching qseqid in the blast_df
    matching_row <- best_matches_newdb[best_matches_newdb$qseqid == fasta_header, ]
    
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
fasta_output <- "annotated_atnarova_full.fasta"

writeLines(modified_fasta, fasta_output)
