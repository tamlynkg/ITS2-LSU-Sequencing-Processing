ITS2-LSU Processing Pipeline
A  bioinformatics pipeline for processing downstream ITS2 and LSU sequencing data from SCATA pipeline output, designed for taxonomic assignment, fungal community analysis and phylogenetic analyses.

Pipeline Workflow
1. ITS Region Extraction
Tool: ITSx
Input: all_clusters_scataXXXX.fasta
Purpose: Extract LSU and ITS2 regions from clustered sequences

bash# ITSx extraction command
ITSx -i all_clusters_scataXXXX.fasta -o output_prefix --preserve T


2. Chimera Removal
Script: chimeras.py
Purpose: Identify and remove duplicate ITS sequences, retaining the sequence with the highest cluster size


The script identifies ITS sequences that are identical, compares cluster sizes for duplicates, selects the representative with maximum abundance

3. Daughter OTU Removal
Tool: MUMU
Purpose: Remove daughter OTUs (sequences that are subsequences of larger, more abundant sequences)

4. Cluster Processing
Script: process_data.R
Purpose: Clean and filter the SCATA cluster file by removing:

- Extra misprimed SCATA clusters
- Zero abundance sequences
- Identified chimeras
- Daughter OTUs
- Non-fungal clusters

5. Taxonomic Assignment
Tool: MASSblast with UNITE database
Purpose: Perform BLAST search against the UNITE fungal database for taxonomic identification


6. BLAST Results Processing
Script: process_blast_results.R

Purpose:
- Filter assignments to 98% similarity threshold
- Extract taxonomic information (species, phylum, etc.)
- Generate assignment summary statistics

7. LSU Sequence Selection
Script: selectLSUsequences.R
Purpose: Select corresponding LSU sequences based on processed ITS2 results

8. LSU Annotation
Script: annotatingLSUsequences.R
Purpose: Annotate selected LSU sequences with taxonomic information

9. Multiple Sequence Alignment
Tool: MAFFT
Script: mafft.sh
Purpose: Generate multiple sequence alignment of LSU sequences for phylogenetic analysis
bash# Example MAFFT command
mafft --auto input_sequences.fasta > aligned_sequences.fasta

10. Phylogenetic Tree Construction
Tool: RAxML
Script: raxml.sh
Purpose: Construct maximum likelihood phylogenetic tree from aligned LSU sequences
bash# Example RAxML command

11. Phyloseq Analysis
Script: phyloseq-analyses.R
Purpose: Generate phyloseq object for downstream ecological and statistical analyses

```
File Structure
pipeline/
├── README.md
├── config.yaml                    # Pipeline configuration
├── snakefile                     # Snakemake workflow file
├── scripts/
│   ├── chimeras.py              # Chimera detection and removal
│   ├── process_data.R           # Cluster file processing
│   ├── process_blast_results.R  # BLAST output processing
│   ├── selectLSUsequences.R     # LSU sequence selection
│   ├── annotatingLSUsequences.R # LSU annotation
│   ├── phyloseq-analyses.R      # Phyloseq object creation
│   ├── mafft.sh                 # Multiple sequence alignment
│   └── raxml.sh                 # Phylogenetic tree construction
├
└── results/
    ├── processed_sequences/
    ├── blast_results/
    ├── alignments/
    ├── trees/
    └── phyloseq_objects/
```


Requirements
Software Dependencies

```
ITSx (>= 1.1.3)
Python (>= 3.7) with BioPython
R (>= 4.0) with packages:

phyloseq
dplyr
tidyverse
Biostrings


MUMU
BLAST+ (>= 2.10)
MAFFT (>= 7.0)
RAxML (>= 8.2)

Databases

UNITE database (latest release recommended)
Configure database paths in config.yaml

```

Clone this repository:

```sh
git clone https://github.com/yourusername/its2-lsu-pipeline.git
cd its2-lsu-pipeline
```


Install dependencies (recommended using conda):

```sh
conda env create -f environment.yml
conda activate its2-lsu-pipeline
```

Configure database paths in config.yaml

Usage
Basic Usage

```sh
snakemake --cores 8
```

Output
The pipeline generates:

Filtered and annotated ITS2 sequences
Corresponding LSU sequences
Taxonomic assignments with confidence scores
Multiple sequence alignments
Phylogenetic trees
Phyloseq objects ready for ecological analysis


Citation
If you use this pipeline in your research, please cite:

