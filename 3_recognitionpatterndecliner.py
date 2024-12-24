import pandas as pd
from Bio.Seq import Seq
from Bio.Data import IUPACData
import itertools
import argparse

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
    processed_sequences = set()  # Track (enzyme, sequence) to allow duplicates with different enzymes

    # Read the input TSV file
    df = pd.read_csv(input_file, sep='\t')

    for index, row in df.iterrows():
        cut_sequence = row['Cut Sequence']
        enzyme = row['Enzyme']
        recognition_sequence = row['Recognition Sequence']

        # Ensure cut sequence contains '/'
        if '/' not in cut_sequence:
            print(f"Warning: Enzyme '{enzyme}' skipped because its Cut Sequence does not contain '/'.")
            continue

        # Calculate the cut position
        cut_pos = cut_sequence.index('/')

        # Generate all possible sequences from the IUPAC code
        possibilities = generate_iupac_possibilities(recognition_sequence)

        for seq in possibilities:
            seq_obj = Seq(seq)
            rc_seq = str(seq_obj.reverse_complement())
            seq_size = len(seq)

            # Skip if the enzyme-sequence pair or its reverse complement has already been processed
            if (enzyme, seq) in processed_sequences or (enzyme, rc_seq) in processed_sequences:
                continue

            # Add the original sequence to results
            results.append({
                'Enzyme': enzyme,
                'Recognition Sequence': seq,
                'Recognition Sequence Size': seq_size,
                'Cut Position': cut_pos  # For sequences starting with '/', cut_pos is 0
            })
            processed_sequences.add((enzyme, seq))

            # Add the reverse complement if different
            if seq != rc_seq:
                reverse_cut_pos_value = seq_size - cut_pos  # Cut position on the reverse complement
                results.append({
                    'Enzyme': enzyme,
                    'Recognition Sequence': rc_seq,
                    'Recognition Sequence Size': seq_size,
                    'Cut Position': reverse_cut_pos_value
                })
                processed_sequences.add((enzyme, rc_seq))

    # Convert results to a DataFrame
    output_df = pd.DataFrame(results)

    # Save the output as a new TSV file
    output_df.to_csv(output_file, sep='\t', index=False)
    print(f"Processed TSV file saved as {output_file}")



if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Process enzyme recognition sites from an input TSV file.")
    parser.add_argument("-input_file", help="Path to the input TSV file.")
    parser.add_argument("-output_file", help="Path to the output TSV file.")
    
    args = parser.parse_args()
    
    process_enzyme_sites(args.input_file, args.output_file)
