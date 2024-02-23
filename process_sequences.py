from collections import defaultdict

# Read data from file
data = {}
with open('Andreas.ITS2.fasta', 'r') as file:
    lines = file.readlines()

# Process the data
current_id = ''
for line in lines:
    line = line.strip()
    if line.startswith('>'):
        current_id = line[1:]
        data[current_id] = ''
    else:
        data[current_id] += line

# Group sequences by content
groups = defaultdict(list)
for identifier, sequence in data.items():
    groups[sequence].append(identifier)

# Output grouped sequences and their identifiers
for sequence, identifiers in groups.items():
    print(f"Sequence: {sequence[:50]}...")  # Print a part of the sequence
    print("Identifiers: ", ", ".join(identifiers))
    print()

