# ITS2-LSU Sequencing Analysis Pipeline

[![Snakemake](https://img.shields.io/badge/snakemake-≥6.0-brightgreen.svg?style=flat)](https://snakemake.readthedocs.io)
[![Python](https://img.shields.io/badge/python-3.7+-blue.svg?style=flat)](https://python.org)
[![R](https://img.shields.io/badge/R-≥4.0-blue.svg?style=flat)](https://r-project.org)
[![License](https://img.shields.io/badge/license-MIT-green.svg?style=flat)](LICENSE)

> **A comprehensive bioinformatics pipeline for high-throughput fungal community analysis from SCATA pipeline output**

Transform your raw fungal sequencing data into publication-ready phylogenetic insights with this robust, automated pipeline designed for ITS2 and LSU region analysis.

---

## Key Features

- **Automated Workflow**: Complete pipeline from raw sequences to phylogenetic trees
- **Quality Control**: Advanced chimera detection and daughter OTU removal
- **High Precision**: 98% similarity threshold for taxonomic assignments
- **Phylogenetic Analysis**: Maximum likelihood trees with bootstrap support
- **Ecological Ready**: Direct integration with phyloseq for downstream analysis
- **Scalable**: Optimized for HPC environments with parallel processing

---

## 🚀 Quick Start

```bash
# Clone the repository
git clone https://github.com/yourusername/ITS2-LSU-Sequencing-Processing.git
cd ITS2-LSU-Sequencing-Processing

# Configure your data paths in config.yaml
nano config.yaml

# Run the complete pipeline
snakemake --cores 8
```

---

## Pipeline Workflow

Our pipeline implements a sophisticated 11-step analytical workflow optimized for fungal community analysis:

### Step 1: ITS Region Extraction
- **Tool**: ITSx (v≥1.1.3)
- **Function**: Precisely extract ITS2 and LSU regions from clustered sequences
- **Output**: Region-specific FASTA files

```bash
ITSx -i all_clusters_scataXXXX.fasta -o output_prefix --preserve T -t Fungi
```

### Step 2: Chimera Detection & Removal
- **Script**: `scripts/chimeras.py`
- **Function**: Identify and eliminate duplicate ITS sequences
- **Algorithm**: Retains highest abundance representative for each unique sequence

### Step 3: Daughter OTU Removal
- **Tool**: MUMU
- **Function**: Remove subsequences of larger, more abundant sequences
- **Impact**: Reduces artificial diversity inflation

### Step 4: Advanced Cluster Processing
- **Script**: `scripts/process_data.R`
- **Function**: Comprehensive quality filtering removing:
  - Extra misprimed SCATA clusters
  - Zero abundance sequences  
  - Identified chimeras
  - Daughter OTUs
  - Non-fungal sequences

### Step 5: Taxonomic Assignment
- **Tool**: BLAST+ with UNITE database
- **Database**: Latest UNITE fungal taxonomy
- **Parameters**: Optimized for fungal ITS sequences

### Step 6: Results Processing
- **Script**: `scripts/process_blast_results.R`
- **Threshold**: 98% similarity cutoff
- **Output**: High-confidence taxonomic assignments with statistics

### Step 7-8: LSU Sequence Processing
- **Selection**: `scripts/selectingLSUsequences`
- **Annotation**: `scripts/annotatingLSUsequences.R`
- **Function**: Match and annotate LSU sequences to filtered ITS2 results

### Step 9: Multiple Sequence Alignment
- **Tool**: MAFFT (v≥7.0)
- **Method**: Auto-detection of optimal alignment strategy

```bash
mafft --auto input_sequences.fasta > aligned_sequences.fasta
```

### Step 10: Phylogenetic Reconstruction
- **Tool**: RAxML-NG
- **Model**: GTR with 1000 bootstrap replicates
- **Support**: Transfer Bootstrap Expectation (TBE) values

```bash
raxml-ng --msa aligned.fasta --model GTR --threads 16 --bs-trees 1000 --all
```

### Step 11: Ecological Analysis Preparation
- **Script**: `scripts/phyloseq-analyses.R`
- **Output**: Publication-ready phyloseq objects

---

## Repository Structure

```
ITS2-LSU-Sequencing-Processing/
├── README.md                        # This comprehensive guide
├── CLAUDE.md                        # AI assistant instructions  
├── config.yaml                      # Pipeline configuration
├── snakefile                        # Snakemake workflow definition
├── mergingreplicatesanalysis.Rmd    # Analysis notebook
├── scripts/                         # Core pipeline scripts
│   ├── chimeras.py                  # Chimera detection algorithm
│   ├── process_data.R               # Data quality control
│   ├── process_blast_results.R      # Taxonomic assignment processing
│   ├── selectingLSUsequences        # LSU sequence selection
│   ├── annotatingLSUsequences.R     # Sequence annotation
│   ├── phyloseq-analyses.R          # Ecological analysis prep
│   ├── process_unite.R              # UNITE database processing
│   ├── process_sequences.py         # Additional sequence utilities
│   └── testing for snakemake.R     # Pipeline testing utilities
└── data/                           # Input data directory
```

---

## System Requirements

### Core Dependencies
| Tool | Version | Purpose |
|------|---------|---------|
| **ITSx** | ≥1.1.3 | ITS region extraction |
| **Python** | ≥3.7 | Script execution + BioPython |
| **R** | ≥4.0 | Statistical analysis |
| **BLAST+** | ≥2.10 | Taxonomic assignment |
| **MAFFT** | ≥7.0 | Multiple sequence alignment |
| **RAxML-NG** | ≥1.1.0 | Phylogenetic reconstruction |

### R Package Dependencies
```r
# Essential packages
phyloseq, dplyr, tidyverse, Biostrings
```

### Database Requirements
- **UNITE Database**: Latest fungal taxonomy release
- Configure paths in `config.yaml`

---

## Installation & Setup

### 1. Environment Setup
```bash
# Create conda environment (recommended)
conda env create -f environment.yml
conda activate its2-lsu-pipeline

# Or install dependencies manually
conda install -c bioconda itsx blast mafft raxml-ng
conda install -c conda-forge r-base python biopython
```

### 2. Database Configuration
```bash
# Download UNITE database
wget https://files.plutof.ut.ee/public/orig/...
# Configure path in config.yaml
```

### 3. Data Preparation
Place your SCATA output files in the `data/` directory and update paths in `config.yaml`

---

## Usage Examples

### Basic Execution
```bash
# Run complete pipeline
snakemake --cores 8

# Dry run (see execution plan)
snakemake -n

# High-performance execution
snakemake --cores 32 --cluster "sbatch -N 1 -n {threads}"
```

### Module Loading (HPC)
```bash
module load bioinfo-tools ITSx blast MAFFT RAxML-NG R python3
```

---

## Output Products

Your analysis generates a comprehensive suite of results:

- **Quality-filtered ITS2 sequences** with taxonomic assignments
- **Corresponding LSU sequences** with full annotations  
- **High-confidence taxonomic assignments** (≥98% similarity)
- **Multiple sequence alignments** ready for phylogenetic analysis
- **Maximum likelihood phylogenetic trees** with bootstrap support
- **Phyloseq objects** optimized for ecological analysis

---

## Contributing

We welcome contributions! Please see our [contribution guidelines](CONTRIBUTING.md) for:
- Bug reports
- Feature requests  
- Code improvements
- Documentation enhancements

---

## Citation

If this pipeline contributes to your research, please cite:

```bibtex
@software{its2_lsu_pipeline,
  title={ITS2-LSU Sequencing Analysis Pipeline},
  author={[Tamlyn Gangiah]},
  year={2024},
  url={https://github.com/yourusername/ITS2-LSU-Sequencing-Processing}
}
```
