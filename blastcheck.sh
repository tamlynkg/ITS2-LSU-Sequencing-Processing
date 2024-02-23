#!/bin/bash -l

#SBATCH -A naiss2023-5-350
#SBATCH -p core
#SBATCH -n 1
#SBATCH -t 00:50:00
#SBATCH -J blast

module load bioinfo-tools blast/2.14.1+

#make blast database
makeblastdb -in Andreas.ITS2.fasta -dbtype nucl -out Andreas.ITS2

#ITS Eukaryote Database
blastn -db Andreas.ITS2 -query Andreas.ITS2.fasta -out dupseqcheck.txt -outfmt "6 qseqid sseqid slen pident length mismatch gapopen qlen qstart qend sstart send evalue bitscore" -max_target_seqs 1
