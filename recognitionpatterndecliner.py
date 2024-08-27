import pandas as pd
from Bio.Seq import Seq
from Bio.Data import IUPACData
import itertools

def reverse_complement(seq):
    complement = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}
    return ''.join(complement[base] for base in reversed(seq))

def generate_iupac_possibilities(sequence):
    """Generate all possible sequences from an IUPAC encoded sequence."""
    iupac_dict = IUPACData.ambiguous_dna_values
    possibilities = [''.join(p) for p in itertools.product(*(iupac_dict[base] for base in sequence))]
    return possibilities

def process_enzyme_sites(input_file, output_file):
    results = []
    processed_sequences = set()

    # Read the input TSV file
    df = pd.read_csv(input_file, sep='\t')

    for index, row in df.iterrows():
        cut_sequence = row['Cut Sequence']
        enzyme = row['Enzyme']
        recognition_sequence = row['Recognition Sequence']

        # Calculate the cut position
        if '(' in cut_sequence:
            # Skip the cases with '(' in Cut Sequence
            continue
        elif '/' in cut_sequence:
            cut_pos = cut_sequence.index('/')
        else:
            cut_pos = None

        # Generate all possible sequences from the IUPAC code
        possibilities = generate_iupac_possibilities(recognition_sequence)

        for seq in possibilities:
            seq_obj = Seq(seq)
            rc_seq = str(seq_obj.reverse_complement())
            seq_size = len(seq)

            # Skip if the sequence or its reverse complement has already been processed
            if seq in processed_sequences or rc_seq in processed_sequences:
                continue

            # Determine cut position for original sequence
            if cut_pos is not None:
                cut_pos_value = int(cut_pos)
                reverse_cut_pos_value = int(seq_size - cut_pos)
            else:
                cut_pos_value = None
                reverse_cut_pos_value = None

            # Add the original sequence to results
            results.append({
                'Enzyme': enzyme,
                'Recognition Sequence': seq,
                'Recognition Sequence Size': seq_size,
                'Cut Position': cut_pos_value if cut_pos_value is not None else ''
            })
            processed_sequences.add(seq)

            # Only add the reverse complement if it's different from the original sequence
            if seq != rc_seq:
                results.append({
                    'Enzyme': enzyme,
                    'Recognition Sequence': rc_seq,
                    'Recognition Sequence Size': seq_size,
                    'Cut Position': reverse_cut_pos_value if reverse_cut_pos_value is not None else ''
                })
                processed_sequences.add(rc_seq)

    # Convert results to a DataFrame
    output_df = pd.DataFrame(results)
    
    # Save the output as a new TSV file
    output_df.to_csv(output_file, sep='\t', index=False)
    print(f"Processed TSV file saved as {output_file}")

# Example usage
input_file = 'sites_RE.txt'  # Replace with your input file path
output_file = 'processed_enzymes.tsv'  # Replace with your desired output file path

process_enzyme_sites(input_file, output_file)
