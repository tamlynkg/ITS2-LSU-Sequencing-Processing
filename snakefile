
# Load configuration
configfile: "config.yaml"
itsx_output_prefix = "Atnarova" 

rule all:
    input:
        "AtnarovaLSU.raxml.support"

rule run_raxml:
    input:
        aligned_fasta="aligned_LSU.fasta"
    output:
        "AtnarovaLSU.raxml.bestTreeCollapsed",
        "AtnarovaLSU.raxml.bestTree",
        "AtnarovaLSU.raxml.mlTrees",
        "AtnarovaLSU.raxml.support",
        "AtnarovaLSU.raxml.bestModel",
        "AtnarovaLSU.raxml.bootstraps",
        "AtnarovaLSU.raxml.log"
    threads: 16
    shell:
        """
        module load bioinfo-tools RAxML-NG/1.1.0
        raxml-ng --msa {input.aligned_fasta} --model GTR --threads {threads} --bs-metric tbe --seed 2 --bs-trees 1000 --all --prefix AtnarovaLSU
        """

# Define the rule for running MAFFT
rule run_mafft:
    input:
        LSU_file="annotated_atnarova_LSU.fasta"
    output:
        aligned_fasta="aligned_LSU.fasta"
    threads: 1
    shell:
        """
        module load bioinfo-tools MAFFT
        mafft --auto {input.LSU_file} > {output.aligned_fasta}
        """

#Read in LSU and select LSU sequences that pass
rule annotate_LSU:
    input:
        data_file= "bestmatches.csv",
        filteredLSU_file= "sequences-LSU-atnarova.fasta"
    output:
        annotatedLSU="annotated_LSU_full.fasta"
    threads: 1
    script:
        "scripts/annotatingLSUsequences.R"

#Read in LSU and select LSU sequences that pass
rule select_LSU:
    input:
        data_file= "bestmatches.csv",
        LSU_file= "Atnarova.LSU.fasta"
    output:
        filteredLSU="sequences-LSU-atnarova.fasta"
    threads: 1
    script:
        "scripts/selectingLSUsequences.R"

# Define the rule for processing BLAST results
rule process_blast_result:
    input:
        blast_input="blast_source11_newdb10.txt",
        fasta_file="source_11.fasta"
    output:
        annotated_file="annotated_atnarova_full.fasta",
	matches="bestmatches.csv"
    threads: 1
    script:
        "scripts/process_blast_results.R"

# Define the rule for running massblast
rule run_massblast:
    input:
        filteredITSx="source_11.fasta"
    output:
        blast_output="blast.txt"
    threads: 4
    shell:
        """
	module load bioinfo-tools blast
	blastn -db scripts/plutof_fungi_its -query {input.filteredITSx} -out {output.blast_output} -outfmt "6 qseqid sseqid slen pident length mismatch gapopen qlen qstart qend sstart send evalue bitscore" -max_target_seqs 10 -reward 1 -gapextend 2 -penalty -2 -word_size 28 -gapopen 0
        """

# Define the rule for R processing
rule process_data:
    input:
        data_file= "data/all_tag_by_cluster_counts.txt",
        ITSx_positions= "Atnarova.positions.txt",
        labels_file=config["labels_file"],
        fasta_file="nonchimera_Atnarova_ITS2.fasta"
    output:
        filteredcounts_data="countsdata_filtered.csv",
        filteredITSx="source_11.fasta"
    threads: 1
    script:
        "scripts/process_data.R" 

# Define the rule for running duplicate check for chimera removal
rule run_duplicate_check:
    input:
        fungalfasta_file=config["fungalITSx_file"],
        scata_file=config["scata_file"]
    output:
        nonduplicate=config["nonduplicate_ITSx"],
        selectedduplicate=config["selectedduplicate"],
        duplicate=config["duplicate"],
	final_output_fasta="nonchimera_Atnarova_ITS2.fasta"
    threads: 1
    shell:
        """
        module load bioinfo-tools R/4.3.1
        module load bioinfo-tools python3
        module load bioinfo-tools biopython
        python3 scripts/chimeras.py {input.fungalfasta_file} {input.scata_file} {output.nonduplicate} {output.selectedduplicate} {output.duplicate} {output.final_output_fasta}
        """
# Define the rule for running ITSx
rule run_itsx:
    input:
        scatafasta_file=config["scatafastaclusters_file"]      
    output:
        "{itsx_output_prefix}.graph",
        "{itsx_output_prefix}.full.fasta",
        "{itsx_output_prefix}.ITS1.fasta",
        "{itsx_output_prefix}_no_detections.txt",
        "{itsx_output_prefix}.problematic.txt",
        "{itsx_output_prefix}.positions.txt",
        "{itsx_output_prefix}.ITS2.fasta",
        "{itsx_output_prefix}.LSU.fasta",
        "{itsx_output_prefix}_no_detections.fasta",
        "{itsx_output_prefix}.summary.txt"
    threads: 1
    shell:
        """
        module load bioinfo-tools ITSx
        ITSx -i {input.scatafasta_file} -o {itsx_output_prefix} --preserve T -t Fungi --save_regions all
        """






