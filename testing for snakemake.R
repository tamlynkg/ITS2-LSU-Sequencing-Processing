#!/usr/bin/env Rscript
.libPaths(Sys.getenv("R_LIBS_USER"))

library(dplyr)
library(ggplot2)
library(ggsignif)
library(tibble)
library(ape)

# Load the Newick file
tree <- read.tree("updated_treespecies.nwk")

#previous.assignments <- read.csv("scataitslongsequencesmatched-atnarovanew.csv", header = TRUE, sep = ",")

#read in fun guild assignments
funguild <- read.delim("otutable_funguild.guilds.txt", header = TRUE)

#extract the percentage identity from the file
pident <- sub("^([0-9]+\\.?[0-9]*)\\|.*", "\\1", funguild$taxonomy)
funguild$pident <- as.numeric(pident)

#extract clusters with less than 98% match
funguildsp <- funguild %>%
  filter(pident < 98)

#extract the species value
speciesvalue <- str_extract(funguildsp$taxonomy, "s__([^;]+)")
speciesvalueSH <- gsub("^s__", "", speciesvalue)

#update species value since we only want it at genus level
name <- sub("_.+", "_sp", speciesvalueSH)

#funguild$species <- paste0(funguild$OTU_ID, "_", name)

funguildsp$species <- paste0(funguildsp$OTU_ID, "_", name) 

#TO UPDATE THE NAMES ON THE PHYLOGENETIC TREE
#extracting the scata number
extract_numeric <- function(name) {
  strsplit(name, "_")[[1]][2]
}

tree_numeric <- sapply(tree$tip.label, extract_numeric)

# Apply the function to the dataframe
df_numeric <- sapply(funguildsp$species, extract_numeric)

#Updating the branch name if it is different
for (i in seq_along(tree_numeric)) {
  match_index <- which(df_numeric == tree_numeric[i])
  if (length(match_index) > 0) {
    # Update the branch name with the new name from the dataframe
    tree$tip.label[i] <- funguild$species[match_index]
  }
}

#update the tree to a new file
write.tree(tree, file="updated_treespecies.nwk")

#check tip labels
tree$tip.label

#select for only ectomycorrhizal matches from the funguild file
filtered_ecto_guild_only <- funguild[funguild$Guild == '|Ectomycorrhizal|', ]
#select for matches containing ECM in the name
filtered_ecto_guild <- funguild[grepl('Ectomycorrhizal', funguild$Guild), ]

#check the difference
checkecmassing <- filtered_ecto_guild[!filtered_ecto_guild$OTU_ID %in% filtered_ecto_guild_only$OTU_ID, ]

#select for ericoid
filtered_ericoid_guild <- funguild[grepl('Ericoid', funguild$Guild), ]
#filter for agarico
filtered_agarico_counts <- funguild[grepl('Agaricomycetes', funguild$taxonomy), ]

#list of all potential ecm species to grab
#headers <- c(
#  "scata6384_913", "scata6384_383", "scata6384_456", "scata6384_2164", "scata6384_1611", 
#  "scata6384_2420", "scata6384_2661", "scata6384_2674", "scata6384_2691", "scata6384_1820", 
#  "scata6384_285", "scata6384_75", "scata6384_671", "scata6384_327", "scata6384_31", 
#  "scata6384_212", "scata6384_420", "scata6384_226", "scata6384_439", "scata6384_468", 
#  "scata6384_457", "scata6384_412", "scata6384_880", "scata6384_2806", "scata6384_1092", 
#  "scata6384_2893", "scata6384_2493", "scata6384_2505", "scata6384_474", "scata6384_2093", 
#  "scata6384_1875", "scata6384_2272", "scata6384_2392", "scata6384_621", "scata6384_935", 
#  "scata6384_314", "scata6384_2480", "scata6384_727", "scata6384_2606", "scata6384_2159", 
#  "scata6384_2647", "scata6384_2060"
#)

#remove this scata cluster - it is a lichen
countsfiltered <- countsfiltered[!rownames(countsfiltered) %in% "scata6384_2196", ]

#identified as ecm from unite
ecm_species <- c(
  "scata6384_897", "scata6384_1357", "scata6384_2007", 
  "scata6384_3065", "scata6384_2691", "scata6384_2674", "scata6384_2647", "scata6384_75", 
  "scata6384_2480", "scata6384_2661", "scata6384_2480", "scata6384_2272", "scata6384_2420", 
  "scata6384_314", "scata6384_383", "scata6384_621", "scata6384_935", "scata6384_2093", 
  "scata6384_2159", "scata6384_1820",  "scata6384_2060", "scata6384_2164", "scata6384_783", "scata6384_2278", 
  "scata6384_2544", "scata6384_1409"
)

#unique ecm species
combined_speciestograb <- unique(c(filtered_ecto_guild_only$OTU_ID, ecm_species))

#potential saprotrophs
sapstocheck <- c(
  "scata6384_262", "scata6384_400", "scata6384_549", "scata6384_428", "scata6384_496", "scata6384_729",  
  "scata6384_972", "scata6384_777", "scata6384_849", "scata6384_1184", "scata6384_1084", 
  "scata6384_1180", "scata6384_1307", "scata6384_2260", "scata6384_2393", "scata6384_1508", 
  "scata6384_2526", "scata6384_2796", "scata6384_2958", "scata6384_788"
)

checking_species_saps <- c(
  "scata6384_690", "scata6384_711", "scata6384_971", "scata6384_753", "scata6384_880", 
  "scata6384_1126", "scata6384_2806", "scata6384_1203", "scata6384_2893", "scata6384_1930", 
  "scata6384_2026", "scata6384_2747", "scata6384_2849", "scata6384_2874", "scata6384_32", 
  "scata6384_285", "scata6384_268", "scata6384_592", "scata6384_830", "scata6384_420", 
  "scata6384_439", "scata6384_489"
)



#selecting for certain sequences given a list

read_fasta_file <- function(fasta_file) {
  sequences <- readDNAStringSet(fasta_file)
  return(sequences)
}

sequences <- read_fasta_file("source_11_andreas.fasta")

matching_sequences <- sequences[names(sequences) %in% checkecmassing$OTU_ID]

output_fasta_file <- "checkingandreas_saps.fasta"
writeXStringSet(matching_sequences, file=output_fasta_file)

#noblast <- names(sequences)[!names(sequences) %in% filtered_ericoid_guild$OTU_ID]

#filtered list of best matches from blast against unite ref database on command line
bestmatches <- read.csv("bestmatches.csv", header = TRUE, sep = ";")
bestmatches$taxonomy <- bestmatches$qseqid

#filter best match list for Agarico and Basidio
filtered_agarico <- bestmatches[grepl('Agaricomycetes', bestmatches$sseqid), ]
filtered_basidio <- bestmatches[grepl('Basidiomycota', bestmatches$sseqid), ]

#setdiff(ecm_species, filtered_ecto_guild_only$OTU_ID)

#counts table 
countsfiltered <- read.delim("countsdata_filtered_andreas.csv", header = TRUE, sep = ",")
#metadata
metadata <- read.delim("metadata.txt", header = TRUE)
#add metadata data to counts dataframe
countsfiltered_metadata <- merge(countsfiltered, metadata, by = "Sample")

#calculate the total reads per sample
sample.totals <- rowSums(countsfiltered[,-c(1:2)])

#make a filtered counts table for ectos
common_columns <- intersect(combined_speciestograb, colnames(countsfiltered_metadata))
selected_df <- countsfiltered_metadata[ ,common_columns]
selected_df$Sample <- countsfiltered_metadata$Sample
selected_df$Treatment <- countsfiltered_metadata$Treatment

#filter agarico funguild table to remove ectomycorrhizal 
minusecto <- setdiff(filtered_agarico_counts$OTU_ID, names(selected_df))
filtered_rows <- filtered_agarico_counts[filtered_agarico_counts$OTU_ID %in% minusecto, ]

#setdiff(combined_speciestograb, names(countsfiltered_metadata))

#calculate ecm total per sample
ecm.sample.total <- rowSums(selected_df[, 1:81])
ecm.sample.total.df <- as.data.frame(ecm.sample.total)
ecm.sample.total.df$Sample <- countsfiltered_metadata$Sample
ecm.sample.total.df$Treatment <- countsfiltered_metadata$Treatment
ecm.sample.total.df$Density <- countsfiltered_metadata$Density

#calculate total of each ecm species
ecm.species.total <- colSums(selected_df[, 1:81])
#calculate total ecm across all samples
ecm.total <- sum(selected_df[, 1:81])

#proportion of ecm in each sample as a proportion of the total fungal community
ecm.proportion.sample <- (ecm.sample.total / sample.totals)*100  # Might need adjustment for zero sample_totals1

ECM.proportion.sample.df <- data.frame(Sample = countsfiltered_metadata$Sample,
                                       ECMproportion = ecm.proportion.sample,
                                       Treatment = countsfiltered_metadata$Treatment)

#make treatment a factor
ECM.proportion.sample.df$Treatment <- as.factor(ECM.proportion.sample.df$Treatment)

#relevel for stats to have control as the reference point
ECM.proportion.sample.df$Treatment <- relevel(ECM.proportion.sample.df$Treatment, ref = "Control")

#linear model
ECMprop <- lm(ECMproportion ~ Treatment, data = ECM.proportion.sample.df)
summary(ECMprop)

#reorder groups for plot

ECM.proportion.sample.df$Treatment <- factor(ECM.proportion.sample.df$Treatment, 
                                        levels = c("Clearcut", "30% Retained", "60% Retained", "Control"))


pdf("ecmproportion-sig-new.pdf", width = 5.43, height = 3.71) 

#barplot showing signficance and SD
ggplot(ECM.proportion.sample.df, aes(x = Treatment, y = ECMproportion)) +
  theme_minimal() +
  labs(x = "Treatment Group", y = "ECM Proportion (%)") +
  stat_summary(geom = "bar", fun = mean,  position = position_dodge(width = 0.2), width = 0.4, color = "black", size = 0,  aes(fill = "lightblue")) +
  stat_summary(geom = "errorbar", fun.data = mean_se, position = position_dodge(width = 0.1), width = 0.2) +
  geom_signif(comparisons = list(c("Clearcut", "Control")), 
              annotations = "***")   +
  scale_fill_manual(values = "lightgreen")  + 
  geom_signif(stat="identity",
              data=data.frame(x=c(2, 3), xend=c(4, 4),
                              y=c(10, 11), annotation=c("**", "***")),
              aes(x=x,xend=xend, y=y, yend=y, annotation=annotation)) +
  guides(fill = "none") +
  theme(
    axis.line = element_line(color = "black"),
    panel.grid = element_blank(),
    panel.border = element_blank(),
    axis.ticks = element_line(color = "black")
  ) +
  ylim(0, NA) 

dev.off()

#to get average across treatment groups
collapsed.df.ecm <- ECM.proportion.sample.df %>%
    group_by(Treatment) %>%
    summarise(across(-Sample, mean))
  
ggplot(collapsed.df.ecm, aes(x = Treatment, y = ECMproportion)) +
  geom_bar(stat = "identity") +
  labs(title = "ECM Proportion per Sample with SD", x = "Treatment", y = "ECM Proportion") +
  theme_bw() 

ggplot(ECM.proportion.sample.df, aes(x = Treatment, y = ECMproportion)) +
  geom_point(aes(alpha = 5)) +  # Adjust point transparency with 'alpha'
  geom_smooth(method = "lm", se = FALSE) +  # Adds the line of best fit
  labs(title = "Sample ECM Proportion vs Treatment",  # Add a title 
       x = "Treatment",
       y = "ECM Proportion")


############################ SAPROTROPHS #########################

#make counts table of just saprotrophs (might chance when I assign ecm species from tree)
selected_df_sap_basid <- countsfiltered_metadata[, minusecto]
selected_df_sap_basid$Sample <- countsfiltered_metadata$Sample
selected_df_sap_basid$Treatment <- countsfiltered_metadata$Treatment

#calculate saprotroph total per sample
sap.sample.total <- rowSums(selected_df_sap_basid[, 1:223])
sap.sample.total.df <- as.data.frame(sap.sample.total)
sap.sample.total.df$Sample <- countsfiltered_metadata$Sample
sap.sample.total.df$Treatment <- countsfiltered_metadata$Treatment
sap.sample.total.df$Density <- countsfiltered_metadata$Density

#calculate total of each saprotroph species
sap.species.total <- colSums(selected_df_sap_basid[, 1:223])
sap.total <- sum(selected_df_sap_basid[, 1:20])

#saprotroph proportion per sample
sap.proportion.sample <- (sap.sample.total / sample.totals)*100  # Might need adjustment for zero sample_totals1


SAP.proportion.sample.df <- data.frame(Sample = countsfiltered_metadata$Sample,
                                       SAPproportion = sap.proportion.sample,
                                       Treatment = countsfiltered_metadata$Treatment)

#factor treatment level
SAP.proportion.sample.df$Treatment <- as.factor(SAP.proportion.sample.df$Treatment)

#reset reference for statistics
SAP.proportion.sample.df$Treatment <- relevel(SAP.proportion.sample.df$Treatment, ref = "30% Retained")

#linear model
SAPprop <- lm(SAPproportion ~ Treatment, data = SAP.proportion.sample.df)
summary(SAPprop)

#average across treatment groups
collapsed.df.sap <- SAP.proportion.sample.df %>%
  group_by(Treatment) %>%
  summarise(across(-Sample, mean))

#order for plots

SAP.proportion.sample.df$Treatment <- factor(SAP.proportion.sample.df$Treatment, levels = c("Clearcut", "30% Retained", "60% Retained", "Control"))

# Plot the bar plot with the ordered Treatment variable
pdf("sapproportion-sig.pdf",width = 5.43, height = 3.71) 

ggplot(collapsed.df.sap, aes(x = Treatment, y = SAPproportion)) +
  geom_bar(stat = "identity") +
  labs(title = "SAP Proportion per Sample with SD", x = "Treatment", y = "SAP Proportion") +
 theme_bw() 

pdf("sapproportion-sig.pdf",width = 5.43, height = 3.71) 

ggplot(SAP.proportion.sample.df, aes(x = Treatment, y = SAPproportion)) +
  theme_minimal() +
  labs(x = "Treatment Group", y = "Saprotroph Proportion (%)") +
  stat_summary(geom = "bar", fun = mean,  position = position_dodge(width = 0.2), width = 0.4, color = "black", size = 0,  aes(fill = "lightblue")) +
  stat_summary(geom = "errorbar", fun.data = mean_se, position = position_dodge(width = 0.1), width = 0.2) +
  scale_fill_manual(values = "pink")  + 
  guides(fill = "none") +
  theme(
    axis.line = element_line(color = "black"),
    panel.grid = element_blank(),
    panel.border = element_blank(),
    axis.ticks = element_line(color = "black")
  ) +
  ylim(0, NA) 


dev.off()

##############################################
ecm.species.intotal <- ecm.species.total/ecm.totalcountsfiltered_trans <- t(countsfiltered)
countsfiltered_trans <- as.data.frame(countsfiltered_trans)
colnames(countsfiltered_trans) <- countsfiltered$Sample
countsfiltered_trans <- countsfiltered_trans[-c(1,2),]

#round percentage identity
bestmatches$pident <- round(bestmatches$pident, 2)

#to get into the correct format for funguil
bestmatches <- bestmatches %>%
  mutate(taxonomy = paste(pident, sseqid, sep = "|"))

df1 <- rownames_to_column(countsfiltered_trans, var = "qseqid")

merged_df <- merge(df1, bestmatches, by = "qseqid")

merged_df <- merged_df[, -c(19:35)]

names(merged_df)[names(merged_df) == "qseqid"] <- "OTU_ID"

write.table(merged_df, "otutable_funguild_andreas.txt", sep = "\t", row.names = FALSE, quote = FALSE)

#### to get length of ITS and LSU ############################

df <- read.delim("Atnarova.positions.txt", header = FALSE)
colnames(df) <- c("scata_ref", "length", "SSU", "ITS1", "5.8S", "ITS2", "LSU", "sequence")
df$scata_ref <- gsub("_repseq_0", "", df$scata_ref)
df <- df[df$scata_ref %in% blast_file1$V1,]

numbers <- gsub("ITS2: (\\d+)-(\\d+)", "\\1 \\2", df$ITS2)
numbers_matrix <- matrix(unlist(strsplit(numbers, " ")), ncol = 2, byrow = TRUE)

# Convert the matrix to a dataframe and subtract the second number from the first
df$ITS2_difference <- as.numeric(numbers_matrix[,2]) - as.numeric(numbers_matrix[,1])
LSUnumbers <- gsub("LSU: (\\d+)-(\\d+)", "\\1 \\2", df$LSU)

# Split the extracted numbers into two columns
LSU_matrix <- matrix(unlist(strsplit(LSUnumbers, " ")), ncol = 2, byrow = TRUE)

#Calculate length of LSU
df$LSU_difference <- as.numeric(LSU_matrix[,2]) - as.numeric(LSU_matrix[,1])

#Full length of sequences
df$length <- gsub(" bp.", "", df$length)
df$length <- as.numeric(df$length)
df$ITScoverage <- df$ITS2_difference/df$length

##### Process BLAST results ###################
library(data.table)
library(stringr)
library(Biostrings)
library(argparse)
library(dplyr)

update_blast_results <- function(blast_input) {
  blast_df <- fread(blast_input, header = FALSE, sep = "\t")
  # Rename the columns for clarity
  colnames(blast_df) <- c("qseqid", "sseqid", "slen", "pident", "length", "mismatch", "gapopen", "qlen", "qstart", "qend", "sstart", "send", "evalue", "bitscore")
  return(blast_df)
}

blast_input <- "blast.txt"
blast_file <- update_blast_results(blast_input)

#to check coverage
blast_file$qlen_98_percent <- 0.98 * blast_file$qlen   # Calculate the 98% of slen
blast_file$coverage <- blast_file$length/blast_file$slen

#top hit blast file
blast_file1 <- fread("blast_source11_newdb.txt", header = FALSE)
colnames(blast_file1) <- c("qseqid", "sseqid", "slen", "pident", "length", "mismatch", "gapopen", "qlen", "qstart", "qend", "sstart", "send", "evalue", "bitscore")

#blast file with 10 hits
blast_file2 <- fread("blast_source11_newdb10.txt", header = FALSE)
colnames(blast_file2) <- c("qseqid", "sseqid", "slen", "pident", "length", "mismatch", "gapopen", "qlen", "qstart", "qend", "sstart", "send", "evalue", "bitscore")

blast_file2$qlen_98_percent <- 0.98 * blast_file2$qlen   # Calculate the 98% of slen

blast_file2 <- blast_file2 %>%
  filter(length >= qlen_98_percent)

#select above 98% for species identification of refs
blast_file98 <- blast_file2 %>%
  filter(pident >= 98)

filtered_df <- blast_file98[grepl("refs", blast_file98$sseqid), ]

setorder(filtered_df, evalue, -bitscore, -pident)
best_matches_newdb_refs <- filtered_df[, .SD[1], by = qseqid]

#extract species values
speciesvalue_refs <- str_extract(best_matches_newdb_refs$sseqid, "s__([^;]+)")
speciesvalueSH_refs <- gsub("^s__", "", speciesvalue_refs)
best_matches_newdb_refs$species <- speciesvalueSH_refs

#extract phylum values
phylumvalue_refs <- str_extract(best_matches_newdb_refs$sseqid, "p__([^;]+)")
phylumvalueSH_refs <- gsub("^p__", "", phylumvalue_refs)
best_matches_newdb_refs$phylum <- phylumvalueSH_refs

#extract order values
ordervalue_refs <- str_extract(best_matches_newdb_refs$sseqid, "o__([^;]+)")
ordervalueSH_refs <- gsub("^o__", "", ordervalue_refs)
best_matches_newdb_refs$order <- ordervalueSH_refs

andreasassign <- read.delim("taxa-table-nodup.txt")

#filter to get the remainder of the rows not present in best_matches
filtered_df1 <- blast_file2[!blast_file2$qseqid %in% best_matches_newdb_refs$qseqid, ]

#filtered_df1 <- anti_join(blast_file2, best_matches_newdb_refs, by = "sseqid")

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

names(previous.assignments)[names(previous.assignments) == "qseqid"] <- "assingment"
names(previous.assignments)[names(previous.assignments) == "sseqid"] <- "qseqid"

names(andreasassign)[names(andreasassign) == "sseqid"] <- "qseqid"

merged_df_nondup <- merge(best_matches_newdb, andreasassign, by = "qseqid")
merged_df_nondup <- merged_df_nondup[, c("qseqid", "sseqid", "Genus_Species")]

merged_df <- merge(best_matches_newdb, previous.assignments, by = "qseqid")

merged_df <- merged_df[, c("qseqid", "sseqid", "assingment")]

andreasassign$qseqid <- andreasassign$sseqid

counts.data2 <- ecm.its.long[, c("scata6384_1257",  "scata6384_897")]

anti_join(andreasassign, merged_df, by="qseqid")
#merged_df <- merge(best_matches_newdb_notrefs, best_matches_newdb, by = "qseqid", suffixes = c("_notrefs", "_refs"))
#diff_sseqid <- merged_df[merged_df$sseqid_notrefs != merged_df$sseqid_refs, ]
#diff_sseqid <- diff_sseqid[, c("qseqid", "sseqid.x", "Genus_Species")]

#scatadf <- best_matches[, 1:2]

read_fasta_file <- function(fasta_file) {
  sequences <- readDNAStringSet(fasta_file)
  return(sequences)
}

sequences <- read_fasta_file("source_11.fasta")

noblast <- names(sequences)[!names(sequences) %in% best_matches_newdb$qseqid]

#filtered_fasta <- sequences[!names(sequences) %in% noblast]

#writeXStringSet(filtered_fasta, filepath = "source_11.fasta")


#to annotate the ITS fasta file  with sequence names
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


warnings()
# Write the modified FASTA data to a new ile
fasta_output <- "annotated_atnarova_full.fasta"

writeLines(modified_fasta, fasta_output)

#to grab LSU sequences that were not annotated because they were unknown

read_fasta_file <- function(fasta_file) {
  sequences <- readDNAStringSet(fasta_file)
  new_headers <- sub("_repseq_0", "", names(sequences))
  names(sequences) <- new_headers
  return(sequences)
}

sequences <- read_fasta_file("Andreas.LSU.fasta")

columns_to_remove <- setdiff(colnames(countsfiltered[, 3:2261]), names(sequences))

#exact_sequences <- intersect(names(sequences), colnames(countsfiltered))
#matching_sequences <- sequences[!names(sequences) %in% exact_sequences]
#output_fasta_file <- "sequences-LSU-atnarova.fasta"
#writeXStringSet(matching_sequences, file=output_fasta_file)
#olumns_to_remove <- setdiff(names(matching_sequences), names(sequences))

matching_sequences <- sequences[names(sequences) %in% columns_to_remove]

output_fasta_file <- "sequences-LSU-unknown-Andreas.fasta"
writeXStringSet(matching_sequences, file=output_fasta_file)

#to annotate LSU file   
#read file
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
fasta_output <- "annotated_atnarova_LSU.fasta"
writeLines(modified_fasta, fasta_output)
