
# Load configuration
configfile: "config.yaml"
itsx_output_prefix = "Mineral" 

rule all:
    input:
        "80relabunfullseq.raxml.support"

# Define the rule for running RAxML
rule run_raxml:
    input:
        aligned_fasta=config["aligned_fasta_fullseq"]
    output:
        "80relabunfullseq.raxml.bestTreeCollapsed",
        "80relabunfullseq.raxml.bestTree",
        "80relabunfullseq.raxml.mlTrees",
        "80relabunfullseq.raxml.support",
        "80relabunfullseq.raxml.bestModel",
        "80relabunfullseq.raxml.bootstraps",
        "80relabunfullseq.raxml.log"
    shell:
        """
        module load bioinfo-tools RAxML-NG/1.1.0
        raxml-ng --msa {input.aligned_fasta} --model GTR --threads 16 --bs-metric tbe --seed 2 --bs-trees 1000 --all --prefix 80relabunfullseq
        """

# Define the rule for running MAFFT
rule run_mafft:
    input:
        annotated_file=config["updateannotated_fasta_output"]
    output:
        aligned_fasta=config["aligned_fasta_fullseq"]
    shell:
        """
        module load bioinfo-tools MAFFT
        mafft --auto {input.annotated_file} > {output.aligned_fasta}
        """

rule remove_characters:
    input:
        annotated_fasta=config["annotated_fasta_fullseq"]
    output:	
        remcharoutput=config["updateannotated_fasta_output"]
    shell:
        """
        module load bioinfo-tools python3
        python3 removecharcs.py {input.annotated_fasta} {output.remcharoutput}
        """

# Define the rule for processing BLAST results
rule process_unite_result:
    input:
        blast_input=config["unite_results_file"],
        scatafasta_file=config["scatafastaclusters_file"],
	SIblastresult=config["processed_blast_output"]
    output:
        fasta_output=config["annotated_fasta_fullseq"]
    script:
        "scripts/process_unite.R"

# Define the rule for running BLAST
rule run_uniteblast:
    input:
        fasta_file=config["unknown_sequences"],
        database=config["UNITE_database"]
    output:
        blast_output=config["unite_results_file"]
    shell:
        """
        module load bioinfo-tools blast
        makeblastdb -in {input.database} -dbtype nucl -out UNITEdb
 blastn -db UNITEdb -query {input.fasta_file} -out {output.blast_output} -outfmt "6 qseqid sseqid slen pident length mismatch gapopen qlen qstart qend sstart send evalue bitscore" -max_target_seqs 10
        """

# Define the rule for processing BLAST results
rule process_blast_result:
    input:
        blast_input=config["blast_results_file"],
        processdata="processed_data.csv",
	fasta_file=config["eightyfungal_fasta"]
    output:
        processed_blast=config["processed_blast_output"],
	unknown_fasta=config["unknown_sequences"]
    script:
        "scripts/process_blast_results.R"

# Define the rule for running BLAST
rule run_blast:
    input:
        fasta_file=config["eightyfungal_fasta"],
        database=config["SI_database"]
    output:
        blast_output=config["blast_results_file"]
    shell:
        """
        module load bioinfo-tools blast
        makeblastdb -in {input.database} -dbtype nucl -out SIfungaldatabase
 blastn -db SIfungaldatabase -query {input.fasta_file} -out {output.blast_output} -outfmt "6 qseqid sseqid slen pident length mismatch gapopen qlen qstart qend sstart send evalue bitscore" -max_target_seqs 1
        """

rule process_data:
    input:
        data_file= "data/all_tag_by_cluster_counts.txt",
        labels_file=config['labels_file'],
        fasta_file=config["nonduplicate_ITSx"]
    output:
        processed_data="processed_data.csv",
        ITSfasta_file=config["eightyfungal_fasta"]
    params:
        fasta_file=config["nonduplicate_ITSx"],
        labels_file=config['labels_file']
    script:
        "scripts/process_data.R" 

# Define the rule for running duplicate check for chimera removal
rule run_duplicate_check:
    input:
        fungalfasta_file=config["fungalITSx_file"]
    output:
        nonduplicate=config["nonduplicate_ITSx"]
    shell:
        """
        module load bioinfo-tools R/4.3.1
        module load bioinfo-tools python3
        python3 removeduplicatesunspecific.py {input.fungalfasta_file} {output.nonduplicate}
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
        "{itsx_output_prefix}_no_detections.fasta",
        "{itsx_output_prefix}.summary.txt"
    shell:
        """
        module load bioinfo-tools ITSx
        ITSx -i {input.scatafasta_file} -o {itsx_output_prefix} --preserve T -t Fungi
        """






