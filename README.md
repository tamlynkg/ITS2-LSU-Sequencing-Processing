# ITS2-LSU-Sequencing-Processing
Processing downstream ITS2 and LSU data from the output of the SCATA pipeline

1)	Perform ITSx: to extract LSU and ITS2 from all_clusters_scataXXXX.fasta
2)	Remove chimeras: these are duplicate ITS sequences – first identify ITS sequences that are the same and then select the one with the highest cluster size
3)	MUMU to remove daughter OTU
4)	Process scata cluster file – remove any extra misprimed scata clusters from the sequences and remove any 0 abundance sequences (YES) + chimeras + daughter OTUs identified above and non-fungal clusters
5)	Perform blast using UNITE database (MASSblast)
6)	Process blast results – filter to get 98% assignments and extract the species, phylum assignments etc…
7)	Select corresponding LSU sequences
8)	Annotate LSU sequences
9)	 Align
10)	Make a phylogenetic tree
11)	Perform phyloseq
