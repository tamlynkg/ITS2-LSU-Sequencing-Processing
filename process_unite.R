#Blast output from unite
unite_blast <- read.table("unknowns.txt", header = FALSE, sep = "\t")
colnames(unite_blast) <- c("qseqid", "sseqid", "slen", "pident", "length", "mismatch", "gapopen", "qlen", "qstart", "qend", "sstart", "send", "evalue", "bitscore")
unite_blast$qlen_98_percent <- 0.98 * unite_blast$qlen    # Calculate the 98% of slen

#Removing the rows that didn't make the criteria
unite_blast <- unite_blast[unite_blast$length >= unite_blast$qlen_98_percent & unite_blast$pident >= 98, , drop = FALSE]
unite_blast <- unite_blast[!duplicated(unite_blast$qseqid), ]

lichen <- c(
  "Cladonia arbuscula subsp. beringiana",
  "Cladonia borealis",
  "Cladonia borealis",
  "Cladonia cenotea",
  "Cladonia deformis",
  "Cladonia gracilis (coll.)",
  "Cladonia gracilis (coll.)",
  "Cladonia phyllophora",
  "Cladonia pleurota",
  "Cladonia portentosa",
  "Cladonia rangiferina",
  "Cladonia rangiformis",
  "Bryoria kockiana",
  "Parmelia sulcata",
  "Hypogymnia physodes",
  "Lecanorales",
  "Leptogium",
  "Placynthiella",
  "Placynthiella dasaea/icmalea",
  "Placynthiella oligotropha",
  "Trapeliopsis granulosa",
  "Trapeliopsis",
  "Trapeliales",
  "Trapeliales",
  "Trapeliales",
  "Lichenomphalia ericetorum"
)

if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("stringr")
library(stringr)

#extract species values from UNITE results
speciesvalue <- str_extract(unite_blast$sseqid, "s__([^;]+)")
speciesvalueSH <- gsub("^s__", "", speciesvalue)
speciesvalue <- gsub("\\|.*", "", speciesvalueSH)
#species <- sub("\\.g__.*", "", speciesvalue)
unite_blast$sseqid <- speciesvalue

#Get species names into acceptable fromat
unite_blast$sseqid <- gsub("_sp", " ", unite_blast$sseqid)
unite_blast$sseqid <- gsub("_", " ", unite_blast$sseqid)

#Remove lichens
unite_blast <- unite_blast[!unite_blast$sseqid %in% lichen, ]

# Update the 'sseqid' column in SIUNITEtaxa with values from unite_dftest

SIUNITEtaxa <- read.csv("processed_SIblast_results.csv", header = TRUE)

merged_data <- merge(SIUNITEtaxa, unite_massblast, by = "qseqid", all.x = TRUE)

merged_data$Update <- ifelse(!is.na(merged_data$sseqid.y), 
                                merged_data$sseqid.y, 
                                merged_data$sseqid.x)

qseqid <- merged_data$qseqid

fasta_data <- readLines(snakemake@input[['scatafasta_file']])

# Modify sequence identifiers
fasta_data <- gsub("_repseq_0", "", fasta_data)

# Write the modified data back to a new file
writeLines(fasta_data, snakemake@input[['scatafasta_file']])

fullsequences <- readDNAStringSet(snakemake@input[['scatafasta_file']])

matching_fullsequences <- fullsequences[names(fullsequences) %in% qseqid]

# Read the BLAST output file into a df

# Read modified FASTA file

output_fasta_file <- "80relabunfullseq.fasta"
writeXStringSet(matching_fullsequences, file=output_fasta_file)
fasta_data <- readLines(output_fasta_file)

# Initialize variables to store the modified FASTA sequences
modified_fasta <-  character(0)

# Loop through the FASTA data and update the headers
i <- 1

while (i <= length(fasta_data)) {
  if (startsWith(fasta_data[i], ">")) {
    # Extract the header from the current line
    fasta_header <- substring(fasta_data[i], 2)  # Remove the ">" character

    # Search for a matching qseqid in the blast_df
    matching_row <- scatadf[scatadf$qseqid == fasta_header, ]

    if (nrow(matching_row) > 0) {
      # Extract qseqid and sseqid from the matching row in blast_df
      qseqid <- matching_row$qseqid
      sseqid <- matching_row$sseqid
      # Modify the header to include qseqid, sseqid, and Cluster.Size
      fasta_data[i] <- paste0(">", qseqid, "_", sseqid)

    }
  }
  # Add the line to the modified FASTA data
    modified_fasta <- c(modified_fasta, fasta_data[i])
  i <- i + 1
}

# Write the modified FASTA data to a new file
fasta_output <- "annotated80relabunfullseq.fasta"
writeLines(modified_fasta, fasta_output)



