#Load Necessary Packages

#if (!require("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")
#BiocManager::install("fastmap")
#install.packages("fastmap")
#install.packages("phytools")
#install.packages("ape")
#devtools::install_github("pmartinezarbizu/pairwiseAdonis/pairwiseAdonis")
library(phyloseq)
library(ggplot2)
library(ape)
library(phytools)
library(fastmap)
library(pairwiseAdonis)

#Load tree data 
phylo_tree <- read.tree("L4BIdfqK-lA7ilRY-IMybA_newick.txt")

#Root tree
rooted_tree <- root(phylo_tree, "scata6384_2606", resolve.root = TRUE)
#common_ancestor <- getMRCA(phylo_tree, "scata6384_2544")
#physeq <- set_tree_root(physeq, "node_label")

#Plot tree
plot(rooted_tree)

# Write the rooted tree to a file
write.tree(rooted_tree, file="path_to_rooted_tree_file")

#Load tree as phyloseq object
treephy <- phy_tree(rooted_tree)
treephy$tip.label

#Create Sample Data Object From Metadata
#Needs same order across files
custom_order <- c("S1", "S10", "S11", "S12", "S13", "S14", "S15", "S16", "S17", "S18", "S19", "S2", "S20", "S3", "S4", "S5", "S6", "S7", "S8", "S9")

#Load metadata
metadata <- read.delim("metadata.txt", header = TRUE)

#Order metadata
metadata$Sample <- factor(metadata$Sample, levels = custom_order)
metadata <- metadata[order(metadata$Sample), ]
#metadata <- row.names(metadata$Sample)
row_names <- metadata[, 1]
# Remove the first row
metadata <- metadata[, -1]

# Set the row names
rownames(metadata) <- row_names
sample_data <- sample_data(metadata)

# Make OTU table by transposing scata table
df.phyloseq <- selected_df[, 1:78]

#transpose and turn back into dataframe object
df_transposed <- t(df.phyloseq)
df_transposed <- as.data.frame(df_transposed)
#add samples as column names
colnames(df_transposed) <- selected_df$Sample

#df_transposed <- df_transposed[-1, ]
#colnames(df_transposed) <- gsub("\\(", "", colnames(df_transposed))
#df_transposed <- as.data.frame(lapply(df_transposed, as.numeric))
#ecm.its.long2 <- ecm.its.long[,-1]
#colnames(df_transposed) <- rownames(ecm.its.long)
#row.names2 <- as.numeric(sub(".*_", "", rownames(df_transposed)))
#sorted_row_names <- rownames(df_transposed)[order(row.names2)]
#df_transposed <- df_transposed[sorted_row_names, ]

#data type as ,atrox
df_matrix <- as.matrix(df_transposed)

#Load otu table as a phyloseq object
otu_table <- otu_table(df_matrix, taxa_are_rows = T)


###############################
#Create OTU table for saprotrophs

df.phyloseq.saps <- selected_df_sap_basid[, 1:223]

#transpose and turn back into dataframe object
df_transposed.saps <- t(df.phyloseq.saps)
df_transposed.saps <- as.data.frame(df_transposed.saps)

#add samples as column names
colnames(df_transposed.saps) <- selected_df$Sample
df_matrix.saps <- as.matrix(df_transposed.saps)

#Load otu table as a phyloseq object
otu_table.saps <- otu_table(df_matrix.saps, taxa_are_rows = T)

#extract species values
speciesvalue <- str_extract(filtered_rows$taxonomy, "s__([^;]+)")
speciesvalueSH <- gsub("^s__", "", speciesvalue)

#best_matches_newdb_reps$species <- sub("_.*", "", best_matches_newdb_reps$species)
#best_matches_newdb_reps$species <- paste0(best_matches_newdb_reps$species, "_sp_REPS")

#extract phylum value
phylumvalue <- str_extract(filtered_rows$taxonomy, "p__([^;]+)")
phylumvalueSH <- gsub("^p__", "", phylumvalue)

#extract order value
ordervalue <- str_extract(filtered_rows$taxonomy, "o__([^;]+)")
ordervalueSH <- gsub("^o__", "", ordervalue)

#extract kingdom value
kingdomvaluevalue <- str_extract(filtered_rows$taxonomy, "k__([^;]+)")
kingdomvaluevalueSH <- gsub("^k__", "", kingdomvaluevalue)

#extract class value
classvaluevalue <- str_extract(filtered_rows$taxonomy, "c__([^;]+)")
classvaluevalueSH <- gsub("^c__", "", classvaluevalue)

#extract genus value
genusvaluevalue <- str_extract(filtered_rows$taxonomy, "g__([^;]+)")
genusvaluevalueSH <- gsub("^g__", "", genusvaluevalue)

#extract family value
familyvaluevalue <- str_extract(filtered_rows$taxonomy, "f__([^;]+)")
familyvaluevalueSH <- gsub("^f__", "", familyvaluevalue)

#make into dataframe
taxa_table <- data.frame(OTU_ID = OTU, Kingdom = kingdomvaluevalueSH, Phylum = phylumvalueSH, Class = classvaluevalueSH, Order = ordervalueSH, Family = familyvaluevalueSH, Genus = genusvaluevalueSH, Species = speciesvalueSH)
OTU <- filtered_rows$OTU_ID

rownames(taxa_table) <- taxa_table$OTU_ID
taxa_table <- taxa_table[, c(2:8)]

#setdiff(ecm_species, OTU)


taxa.table <- read.delim("taxa_table_new.txt.txt", header = TRUE, row.names = 1)

taxa.table <- as.matrix(taxa_table)

#Load Taxa Table as Phyloseq object
tax_table <- tax_table(taxa.table)

#checks and balances
rownames(tax_table)
colnames(tax_table)

#setdiff(treephy$tip.label, rownames(tax_table))

rownames(otu_table)
colnames(otu_table)

#MAKE THE PHYLOSEQ OBJECT - add tree here

physeq <- phyloseq(otu_table.saps, sample_data, tax_table)

plot_bar(physeq, fill = "Genus")

#caclulate relative abundance
physeq_rel_abundance <- transform_sample_counts(physeq, function(x) x / sum(x) *100)
otu_table(physeq_rel_abundance)

#plot bar graph
plot_bar(physeq_rel_abundance, fill = "Genus")

genus_rel_abundance <- tax_glom(physeq_rel_abundance, taxrank = "Genus")

# Transform to relative abundances if not already done
genus_rel_abundance <- transform_sample_counts(genus_rel_abundance, function(x) x / sum(x))


# Calculate the average relative abundance for each genus
# Extract the OTU table as a data frame for ggplot2
df <- psmelt(genus_rel_abundance)

# Calculate the mean relative abundance for each genus
mean_abundance <- aggregate(Abundance ~ Genus, data = df, mean)

# Plot the bar plot showing average relative abundance of each genus
#pdf("plotrelabund-top15-saps.pdf", width = 12, height = 7) 
#ggplot(mean_abundance, aes(x = reorder(Genus, -Abundance), y = Abundance, fill = Genus)) +
#  geom_bar(stat = "identity") +
#  labs(x = "Species", y = "Average Relative Abundance", title = "Average Relative Abundance of Genera") +
#  theme(axis.text.x = element_text(angle = 90, hjust = 1))

#top 10 species/genus
top_10_abundance <- mean_abundance %>%
  arrange(desc(Abundance)) 

top_10_abundance <- top_10_abundance[1:5,]

# Create the bar plot
pdf("plotrelabund-top15-saps.pdf", width = 12, height = 7) 

ggplot(top_10_abundance, aes(x = reorder(Genus, -Abundance), y = Abundance, fill = Genus)) +
  geom_bar(stat = "identity") +
  labs(x = "Genus", y = "Average Relative Abundance (%)") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1),
panel.background = element_blank(),
axis.line = element_line(color = "black"),
axis.ticks = element_line(color = "black"))

dev.off()

ggplot(mean_abundance, aes(x = "", y = Abundance, fill = Genus)) +
  geom_bar(stat = "identity", width = 1) +
  coord_polar("y", start = 0) +
  labs(title = "Average Relative Abundance of Genera") +
  theme_void() +
  theme(legend.position = "right")

carbom_fraction <- merge_samples(physeq, "Treatment")
physeq_rel_abundance_pergroup <- transform_sample_counts(carbom_fraction, function(x) x / sum(x))

#pdf("plotrelabund.pdf", width = 12, height = 7) 

plot_bar(GP, fill = "Genus") + 
 geom_bar(aes(color=Genus, fill=Genus), stat="identity")


#dev.off()
#perform hellinger transformation
physeq_rel_abundance_hellinger <- transform_sample_counts(physeq, function(x) sqrt(x))

# Check the transformed OTU table
otu_table(physeq_rel_abundance_hellinger)

### INFORMATION FROM PHYLOSEQ ####

rank_names(tax_table)
ntaxa(physeq_rel_abundance)
nsamples(physeq)
sample_names(physeq)
taxa_names(physeq)
sample_variables(physeq)
get_taxa_unique(physeq, "Genus")

#Alpha diversity - but not necessary

(p = plot_richness(physeq, x = "Treatment"))
p + geom_boxplot(data = p$data, aes(x = Treatment, y = value, color = NULL), alpha = 0.1)


#physeq_rel_abundance_hellinger <- transform_sample_counts(physeq_rel_abundance, function(x) sqrt(x))
#otu_table(physeq_rel_abundance_hellinger)

weighted_unifrac_dist <- UniFrac(physeq_rel_abundance, weighted = TRUE)

#make nmds
nmds <- ordinate(physeq_rel_abundance, method = "NMDS", distance = "weighted_unifrac_dist")
plot_ordination(physeq_rel_abundance, nmds, color = "Treatment", label = "Density") +
  geom_text(aes(label = Density), size = 3, vjust = -0.5)

#most abundant taxa across all
barplot(sort(taxa_sums(physeq_rel_abundance), TRUE)[1:50]/nsamples(physeq_rel_abundance), las=2)

#to show abundance of top genera across treatments
topsp <- names(sort(taxa_sums(physeq_rel_abundance), TRUE)[1:50])
GP    <- prune_taxa(topsp, physeq_rel_abundance)
top5ph <- sort(tapply(taxa_sums(GP), tax_table(GP)[, "Species"], sum), decreasing=TRUE)[1:5]
GP     <- subset_taxa(GP, Species %in% names(top5ph))
plot_bar(GP, x="Treatment", facet_grid= ~ Species)
plot_bar(GP, x="Treatment", facet_grid= ~ Genus_Species)

#DESeq

#install.packages("DESeq2")
library("DESeq2")

prunedf <- prune_samples(sample_sums(physeq) > 100, physeq)
sample_names(prunedf)
sample_data_df <- data.frame(sample_data(physeq_rel_abundance))

library(vegan)

# Perform PERMANOVA
permanova_result <- adonis2(weighted_unifrac_dist ~ Treatment, data = sample_data_df)
#adonis2(weighted_unifrac_dist ~ Density, data = sample_data_df)

#perform non-parametric tests on the axes
nmds_scores <- as.data.frame(scores(nmds))

kruskal.test(nmds_scores$NMDS1 ~ Treatment, data = sample_data_df)
kruskal.test(nmds_scores$NMDS2 ~ Treatment, data = sample_data_df)

#kruskal.test(weighted_unifrac_dist ~ Treatment, data = sample_data_df)

library(emmeans)

#Weighted average to calculate the influence of each axis

rel.abund.df <- data.frame(otu_table(physeq_rel_abundance))
rel.abund.df.t <- as.data.frame(t(rel.abund.df))

WAEaxis1 <- colSums(rel.abund.df.t*nmds_scores$NMDS1)/colSums(rel.abund.df.t)
WAEaxis2 <- colSums(rel.abund.df.t*nmds_scores$NMDS2)/colSums(rel.abund.df.t)

#make dataframe of species scores
species_scores <- data.frame(
  Species = row.names(rel.abund.df),
  NMDS1 = WAEaxis1,  # Species scores on axis 1
  NMDS2 = WAEaxis2  # Species scores on axis 2
)

#merge species data frame with taxa table to get species names
merged_df_speciescore <- merge(species_scores, taxa.table.2, by.x = "Species", by.y = "sseqid", all.x = TRUE)
merged_df_speciescore <- merged_df_speciescore[, -c(4:8,10)]

#add metadata to nmds score dataframe 
nmds_scores$Treatment <- metadata$Treatment
nmds_scores$TreeDensity <- metadata$Density

#plot_ordination(physeq_rel_abundance_hellinger, nmds, color = "", label = "Density") +
#  geom_text(aes(label = Density), size = 3, vjust = -0.5)

#Species scores at genus level

merged_df_speciescore_aggregated <- merged_df_speciescore %>%
  group_by(Genus) %>%
  summarise(
    NMDS1 = mean(NMDS1, na.rm = TRUE),
    NMDS2 = mean(NMDS2, na.rm = TRUE),
  )

#create nmdGenus#create nmds with weighted average species scores
#fontface = "bold"

pdf("plot.pdf", width = 12, height = 7) 

ggplot() +
  geom_point(data = nmds_scores, aes(x = NMDS1, y = NMDS2, size = TreeDensity, fill = Treatment), shape = 21, color = "black", stroke = 0.5) +
  #  stat_ellipse() +
  #  geom_text(data = nmds_scores, aes(x = NMDS1, y = NMDS2),, label = nmds_scores$Treatment, vjust = 1.5, size = 5) +
  #  geom_point(data = merged_df_speciescore_aggregated, aes(x = NMDS1, y = NMDS2), color = "black", size = 2, shape = 15) +
  geom_text(data = merged_df_speciescore_aggregated, aes(x = NMDS1, y = NMDS2, label = Genus), color = "black", size = 6, vjust = 1.5) +
  theme_minimal() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(size = 0.5, color = "black"),
    legend.text = element_text(size = 12),    # Increase legend text size
    legend.title = element_text(size = 14)
  )  +
  scale_fill_manual(
    values = c(
      "30% Retained" = "orange",
      "60% Retained" = "yellow",
      "Clearcut" = "red",
      "Control" = "green"
    )) +
  xlab("NMDS1") +
  ylab("NMDS2")  + scale_size(range = c(5, 20))  


dev.off()

averagetrees <- metadata %>%
  group_by(Treatment) %>%
  summarise(
    Density = mean(Density, na.rm = TRUE)
  )

plot_bar(ps, x = "Sample", fill = "Genus") +
  theme_minimal() +
  labs(x = "Sample", y = "Abundance") +
  theme(
    axis.text.x = element_text(size = 10, angle = 45, hjust = 1),  # Adjust x-axis text size and rotation
    axis.text.y = element_text(size = 10),  # Adjust y-axis text size
    axis.title = element_text(size = 12),  # Adjust axis title size
    plot.title = element_text(size = 14, face = "bold")  # Adjust plot title size and style
  ) +
  scale_x_discrete(expand = c(0, 0))  # Adjust x-axis limits or expand

# Perform pairwise PERMANOVA
sample_data_df$Treatment <- as.character(sample_data_df$Treatment)

pairwise_results <- pairwise.adonis(weighted_unifrac_dist,sample_data_df$Treatment)

sample_data_df$Sample <- row.names(sample_data_df)

#change multiple test correction to FDR
adjusted_p_values <- p.adjust(pairwise_results$p.value, method = "fdr")

# Perform post hoc test
mod <- betadisper(weighted_unifrac_dist, metadata$Treatment)
permutest(mod)
plot(mod)
boxplot(mod)
mod.HSD <- TukeyHSD(mod)
plot(mod.HSD)

nmds <- metaMDS(weighted_unifrac_dist)


