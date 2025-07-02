import sys
import pandas as pd
from Bio import SeqIO, SeqRecord, Seq

def identify_sequences(input_fasta):
    sequences = {}
    duplicates = {}
    # Parse the FASTA file
    for record in SeqIO.parse(input_fasta, "fasta"):
        seq = str(record.seq)
        identifier = record.id.replace("_repseq_0", "") 
        if seq in sequences:
            if seq not in duplicates:
                duplicates[seq] = [sequences[seq]]
            duplicates[seq].append(identifier)
        else:
            sequences[seq] = identifier
    # Separate non-duplicate sequences
    non_duplicates = {seq: identifier for seq, identifier in sequences.items() if seq not in duplicates}
    return duplicates, non_duplicates

def find_and_select(input_csv, duplicates):
    # Load the dataframe
    df = pd.read_csv(input_csv, delimiter=';', skiprows=42, on_bad_lines='skip')
    selected_identifiers = []
    for seq, ids in duplicates.items():
        filtered_df = df[df['Cluster ID'].isin(ids)]
        if not filtered_df.empty:
            largest_cluster_row = filtered_df.loc[filtered_df['Cluster Size'].idxmax()]
            selected_identifiers.append(largest_cluster_row['Cluster ID'])
    return selected_identifiers

def save_fasta(sequences, output_fasta):
    records = []
    for seq, identifier in sequences.items():
        record = SeqRecord.SeqRecord(Seq.Seq(seq), id=identifier, description="")
        records.append(record)
    SeqIO.write(records, output_fasta, "fasta")

def save_selected_duplicates(duplicates, selected_identifiers, output_fasta):
    records = []
    for seq, ids in duplicates.items():
        for identifier in ids:
            if identifier in selected_identifiers:
                record = SeqRecord.SeqRecord(Seq.Seq(seq), id=identifier, description="")
                records.append(record)
    SeqIO.write(records, output_fasta, "fasta")

def save_all_duplicates(duplicates, output_fasta):
    records = []
    for seq, ids in duplicates.items():
        for identifier in ids:
            record = SeqRecord.SeqRecord(Seq.Seq(seq), id=identifier, description="")
            records.append(record)
    SeqIO.write(records, output_fasta, "fasta")

def merge_fasta_files(file1, file2, output_file):
    records = list(SeqIO.parse(file1, "fasta")) + list(SeqIO.parse(file2, "fasta"))
    SeqIO.write(records, output_file, "fasta")

if __name__ == "__main__":
    if len(sys.argv) != 7:
        print("Usage: python script.py input_fasta input_csv output_non_duplicates_fasta output_selected_duplicates_fasta output_all_duplicates_fasta final_output_fasta")
        sys.exit(1)
    
    input_fasta = sys.argv[1]
    input_csv = sys.argv[2]
    output_non_duplicates_fasta = sys.argv[3]
    output_selected_duplicates_fasta = sys.argv[4]
    output_all_duplicates_fasta = sys.argv[5]
    final_output_fasta = sys.argv[6]
    
    # Identify duplicates and non-duplicates in the FASTA file
    duplicates, non_duplicates = identify_sequences(input_fasta)
    
    # Find and select identifiers in the CSV file for duplicates
    selected_identifiers = find_and_select(input_csv, duplicates)
    
    # Save non-duplicate sequences to a new FASTA file
    save_fasta(non_duplicates, output_non_duplicates_fasta)
    
    # Save selected duplicate sequences to a new FASTA file
    save_selected_duplicates(duplicates, selected_identifiers, output_selected_duplicates_fasta)
    
    # Save all duplicate sequences to a new FASTA file
    save_all_duplicates(duplicates, output_all_duplicates_fasta)

    # Merge selected duplicates and non-duplicates into a final FASTA file
    merge_fasta_files(output_non_duplicates_fasta, output_selected_duplicates_fasta, final_output_fasta)

    print(f"Merged FASTA file saved as {final_output_fasta}")
    

