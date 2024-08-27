import csv

# Load enzymes from TSV file
enzymes = []
with open('processed_enzymes.tsv', 'r') as tsvfile:
    reader = csv.DictReader(tsvfile, delimiter='\t')
    for row in reader:
        enzyme = {
            'name': row['Enzyme'],
            'recognition_sequence': row['Recognition Sequence'],
            'cut_position': int(row['Cut Position'])
        }
        enzymes.append(enzyme)

# Load sequences from CSV file
sequences = []
with open('output.csv', 'r') as csvfile:
    reader = csv.DictReader(csvfile)
    for row in reader:
        sequences.append(row)

# Prepare output list
output_data = []

# Find sequences where Sequence 1 is cut but Sequence 2 is not, and vice versa
for enzyme in enzymes:
    for seq in sequences:
        sequence1 = seq['Sequence 1']
        sequence2 = seq['Sequence 2']
        
        # Check if Sequence 1 is cut but not Sequence 2
        if enzyme['recognition_sequence'] in sequence1 and enzyme['recognition_sequence'] not in sequence2:
            found_sequence = enzyme['recognition_sequence']
            start_position_in_sequence1 = sequence1.find(found_sequence)
            cut_position_in_sequence1 = start_position_in_sequence1 + enzyme['cut_position']
            output_data.append({
                'Chromosome': seq['Chromosome'],
                'Position': seq['Position'],
                'Enzyme Name': enzyme['name'],
                'Cut Sequence Number': 1,
                'Cut Position in Sequence': cut_position_in_sequence1,
                'Found Recognition Sequence': found_sequence,
                'All Recognition Sequences': ','.join([e['recognition_sequence'] for e in enzymes if e['name'] == enzyme['name']]),
                'Sequence 1': sequence1,
                'Sequence 2': sequence2
            })
        
        # Check if Sequence 2 is cut but not Sequence 1
        if enzyme['recognition_sequence'] in sequence2 and enzyme['recognition_sequence'] not in sequence1:
            found_sequence = enzyme['recognition_sequence']
            start_position_in_sequence2 = sequence2.find(found_sequence)
            cut_position_in_sequence2 = start_position_in_sequence2 + enzyme['cut_position']
            output_data.append({
                'Chromosome': seq['Chromosome'],
                'Position': seq['Position'],
                'Enzyme Name': enzyme['name'],
                'Cut Sequence Number': 2,
                'Cut Position in Sequence': cut_position_in_sequence2,
                'Found Recognition Sequence': found_sequence,
                'All Recognition Sequences': ','.join([e['recognition_sequence'] for e in enzymes if e['name'] == enzyme['name']]),
                'Sequence 1': sequence1,
                'Sequence 2': sequence2
            })

# Write output to CSV file
with open('cut_sequences_output.csv', 'w', newline='') as csvfile:
    fieldnames = [
        'Chromosome', 
        'Position', 
        'Enzyme Name', 
        'Cut Sequence Number', 
        'Cut Position in Sequence', 
        'Found Recognition Sequence', 
        'All Recognition Sequences',
        'Sequence 1', 
        'Sequence 2'
    ]
    writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
    writer.writeheader()
    writer.writerows(output_data)

print("CSV file 'cut_sequences_output.csv' has been created with the cut sequence information.")
