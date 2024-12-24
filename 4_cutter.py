import csv
import argparse
import sys
import primer3
from Bio.Seq import Seq

# Function to read enzyme data from a file and store in a dictionary
def read_enzymes(enzyme_file_path):
    enzyme_dict = {}
    try:
        with open(enzyme_file_path, 'r') as file:
            # Determine the delimiter (tab or comma)
            first_line = file.readline()
            file.seek(0)  # Reset file pointer to the beginning
            if '\t' in first_line:
                delimiter = '\t'
            else:
                delimiter = ','
            reader = csv.DictReader(file, delimiter=delimiter)
            for row in reader:
                enzyme = row['Enzyme'].strip()
                sequences = row['Recognition Sequence']
                seq_list = [seq.strip().upper() for seq in sequences.split(',')]
                enzyme_dict[enzyme] = seq_list
        return enzyme_dict
    except FileNotFoundError:
        print(f"Error: Enzyme file '{enzyme_file_path}' not found.")
        sys.exit(1)
    except KeyError as e:
        print(f"Error: Missing expected column in enzyme file: {e}")
        sys.exit(1)

# Function to find cuts and the recognition sequences that occur in the sequence
def find_cuts_with_sequences(sequence, recognition_sequences):
    cuts = []
    sequence = sequence.upper()
    for recog_seq in recognition_sequences:
        recog_seq = recog_seq.upper()
        start = 0
        while True:
            index = sequence.find(recog_seq, start)
            if index == -1:
                break
            # Positions are recorded as 1-based indices
            cuts.append({
                'Position': index + 1,
                'RecognitionSequence': recog_seq
            })
            start = index + 1  # Move past this match
    return cuts

# Function to compute fragment sizes given amplicon length and cut positions
def compute_fragment_sizes(amplicon_length, cut_positions):
    positions = [0] + sorted(cut_positions) + [amplicon_length]
    fragment_sizes = []
    for i in range(len(positions) - 1):
        fragment_size = positions[i+1] - positions[i]
        fragment_sizes.append(fragment_size)
    return fragment_sizes

# Function to run Primer3 and design primers for a given sequence
def design_primers(sequence, primer_params):
    # Default Primer3 settings
    primer3_settings = {
        'PRIMER_OPT_SIZE': 20,
        'PRIMER_MIN_SIZE': 18,
        'PRIMER_MAX_SIZE': 25,
        'PRIMER_OPT_TM': 60.0,
        'PRIMER_MIN_TM': 57.0,
        'PRIMER_MAX_TM': 63.0,
        'PRIMER_MAX_POLY_X': 4,
        'PRIMER_SALT_MONOVALENT': 50.0,
        'PRIMER_DNA_CONC': 50.0,
        'PRIMER_NUM_RETURN': 1,
    }

    # Update settings with user-specified parameters
    primer3_settings.update(primer_params)

    seq_args = {
        'SEQUENCE_TEMPLATE': sequence,
    }

    # Use design_primers instead of the deprecated designPrimers
    primer_results = primer3.bindings.design_primers(seq_args, primer3_settings)

    # Optional: Check if there's a general error message
    general_error = primer_results.get('PRIMER_ERROR', '')

    if primer_results['PRIMER_PAIR_NUM_RETURNED'] > 0:
        left_primer = primer_results['PRIMER_LEFT_0_SEQUENCE']
        right_primer = primer_results['PRIMER_RIGHT_0_SEQUENCE']
        return left_primer, right_primer, general_error
    else:
        # If no primers returned, gather any available explanation
        left_explain = primer_results.get('PRIMER_LEFT_EXPLAIN', '')
        right_explain = primer_results.get('PRIMER_RIGHT_EXPLAIN', '')
        pair_explain = primer_results.get('PRIMER_PAIR_EXPLAIN', '')

        # Combine the explanations into one string
        # This might help you see if there's a specific reason (e.g., too small, no valid primer, etc.)
        explanation = f"PRIMER_LEFT_EXPLAIN: {left_explain}; PRIMER_RIGHT_EXPLAIN: {right_explain}; PRIMER_PAIR_EXPLAIN: {pair_explain}; PRIMER_ERROR: {general_error}"

        # Return None for both primers but include an explanation
        return None, None, explanation

# Function to extract amplicons from sequences using given primers
def extract_amplicon(sequence, left_primer_seq, right_primer_seq):
    sequence = sequence.upper()
    left_primer_seq = left_primer_seq.upper()
    right_primer_seq_rc = str(Seq(right_primer_seq).reverse_complement())
    # Find left primer position
    left_pos = sequence.find(left_primer_seq)
    if left_pos == -1:
        return None  # Left primer not found
    # Find right primer position
    right_pos = sequence.find(right_primer_seq_rc)
    if right_pos == -1:
        return None  # Right primer not found
    right_pos_end = right_pos + len(right_primer_seq_rc)
    # Extract amplicon
    amplicon_sequence = sequence[left_pos:right_pos_end]
    return amplicon_sequence

# Function to compare enzyme cuts between amplicons and record positions and sequences
def compare_enzyme_cuts(snp_file, enzyme_dict, primer_params):
    try:
        with open(snp_file, mode='r') as csvfile:
            reader = csv.DictReader(csvfile)
            results = []
            for row in reader:
                seq1 = row['Sequence 1'].upper()
                seq2 = row['Sequence 2'].upper()

                # Design primers based on Sequence 1 only
                left_primer_seq, right_primer_seq, primer_explanation = design_primers(seq1, primer_params)

                if left_primer_seq is None or right_primer_seq is None:
                    # Print a more detailed warning message
                    print(
                        f"Warning: Could not design primers for SNP at Chromosome {row['Chromosome']}, "
                        f"Position {row['Position']}. Explanation: {primer_explanation}"
                    )
                    continue

                # Extract amplicons from both sequences using the same primers
                amplicon_seq1 = extract_amplicon(seq1, left_primer_seq, right_primer_seq)
                amplicon_seq2 = extract_amplicon(seq2, left_primer_seq, right_primer_seq)

                if amplicon_seq1 is None or amplicon_seq2 is None:
                    print(f"Warning: Primers do not anneal properly in one of the sequences at Chromosome {row['Chromosome']}, Position {row['Position']}")
                    continue

                # Check if amplicons are the same size
                if len(amplicon_seq1) != len(amplicon_seq2):
                    print(f"Warning: Amplicon sizes differ for SNP at Chromosome {row['Chromosome']}, Position {row['Position']}")
                    continue

                amplicon_size = len(amplicon_seq1)

                # Perform RFLP analysis on the amplicons
                for enzyme, recognition_seqs in enzyme_dict.items():
                    cuts_seq1 = find_cuts_with_sequences(amplicon_seq1, recognition_seqs)
                    cuts_seq2 = find_cuts_with_sequences(amplicon_seq2, recognition_seqs)

                    # Get positions as sets for set operations
                    positions_seq1_set = set(cut['Position'] for cut in cuts_seq1)
                    positions_seq2_set = set(cut['Position'] for cut in cuts_seq2)
                    positions_seq1 = sorted(positions_seq1_set)
                    positions_seq2 = sorted(positions_seq2_set)

                    # Compute fragment sizes
                    fragment_sizes_seq1 = compute_fragment_sizes(amplicon_size, positions_seq1)
                    fragment_sizes_seq2 = compute_fragment_sizes(amplicon_size, positions_seq2)

                    # Convert cut positions to strings
                    cuts_positions_seq1_str = ';'.join(map(str, positions_seq1))
                    cuts_positions_seq2_str = ';'.join(map(str, positions_seq2))

                    # Check if the difference in the number of cuts is exactly one
                    num_cuts_seq1 = len(positions_seq1)
                    num_cuts_seq2 = len(positions_seq2)
                    cut_difference = abs(num_cuts_seq1 - num_cuts_seq2)

                    if cut_difference == 1:
                        # Find unique cuts
                        unique_cuts_positions = positions_seq1_set.symmetric_difference(positions_seq2_set)

                        unique_cuts_details = []

                        for pos in unique_cuts_positions:
                            if pos in positions_seq1_set:
                                # Cut is unique to Amplicon 1
                                cut_amplicon = 'Amplicon 1'
                                cut_seq = next(cut['RecognitionSequence'] for cut in cuts_seq1 if cut['Position'] == pos)
                                corresponding_seq = amplicon_seq2[pos - 1: pos - 1 + len(cut_seq)]
                            else:
                                # Cut is unique to Amplicon 2
                                cut_amplicon = 'Amplicon 2'
                                cut_seq = next(cut['RecognitionSequence'] for cut in cuts_seq2 if cut['Position'] == pos)
                                corresponding_seq = amplicon_seq1[pos - 1: pos - 1 + len(cut_seq)]

                            # Ensure sequences are uppercase
                            cut_seq = cut_seq.upper()
                            corresponding_seq = corresponding_seq.upper()

                            unique_cuts_details.append(f"Position: {pos}, {cut_amplicon} Cut: {cut_seq}, Other Amplicon Sequence: {corresponding_seq}")

                        # Combine unique cuts details into a single string
                        unique_cuts_str = ' | '.join(unique_cuts_details)

                        results.append({
                            'Chromosome': row['Chromosome'],
                            'Position': row['Position'],
                            'Enzyme': enzyme,
                            'Left Primer': left_primer_seq,
                            'Right Primer': right_primer_seq,
                            'Amplicon Size': amplicon_size,
                            'Amplicon Sequence 1': amplicon_seq1,
                            'Amplicon Sequence 2': amplicon_seq2,
                            'Cuts in Amplicon 1': num_cuts_seq1,
                            'Cuts in Amplicon 2': num_cuts_seq2,
                            'Positions in Amplicon 1': cuts_positions_seq1_str,
                            'Positions in Amplicon 2': cuts_positions_seq2_str,
                            'Fragment Sizes Amplicon 1': ';'.join(map(str, fragment_sizes_seq1)),
                            'Fragment Sizes Amplicon 2': ';'.join(map(str, fragment_sizes_seq2)),
                            'Unique Cuts': unique_cuts_str,
                        })
                    # Else, skip enzymes that do not meet the criteria

            return results
    except FileNotFoundError:
        print(f"Error: SNP file '{snp_file}' not found.")
        sys.exit(1)
    except KeyError as e:
        print(f"Error: Missing expected column in SNP file: {e}")
        sys.exit(1)

# Main execution
def main():
    # Set up argument parsing
    parser = argparse.ArgumentParser(description='Perform RFLP analysis on hypothetical amplicons generated by Primer3.')
    parser.add_argument('snp_file', help='Path to the SNP CSV file.')
    parser.add_argument('enzyme_file', help='Path to the enzyme data file.')
    parser.add_argument('-o', '--output', default='enzyme_cut_amplicons.csv', help='Output CSV file name.')
    # Primer3 parameters
    parser.add_argument('--PRIMER_OPT_SIZE', type=int, default=20, help='Optimal primer size.')
    parser.add_argument('--PRIMER_MIN_SIZE', type=int, default=18, help='Minimum primer size.')
    parser.add_argument('--PRIMER_MAX_SIZE', type=int, default=25, help='Maximum primer size.')
    parser.add_argument('--PRIMER_OPT_TM', type=float, default=60.0, help='Optimal primer melting temperature (Tm).')
    parser.add_argument('--PRIMER_MIN_TM', type=float, default=57.0, help='Minimum primer Tm.')
    parser.add_argument('--PRIMER_MAX_TM', type=float, default=63.0, help='Maximum primer Tm.')
    parser.add_argument('--PRIMER_PRODUCT_SIZE_RANGE', type=str, default='100-300', help='Desired PCR product size range (e.g., "100-300").')
    args = parser.parse_args()

    # Collect Primer3 parameters
    primer_params = {
        'PRIMER_OPT_SIZE': args.PRIMER_OPT_SIZE,
        'PRIMER_MIN_SIZE': args.PRIMER_MIN_SIZE,
        'PRIMER_MAX_SIZE': args.PRIMER_MAX_SIZE,
        'PRIMER_OPT_TM': args.PRIMER_OPT_TM,
        'PRIMER_MIN_TM': args.PRIMER_MIN_TM,
        'PRIMER_MAX_TM': args.PRIMER_MAX_TM,
        'PRIMER_PRODUCT_SIZE_RANGE': [[int(x) for x in args.PRIMER_PRODUCT_SIZE_RANGE.split('-')]]
    }

    # Read enzyme data
    enzyme_dict = read_enzymes(args.enzyme_file)

    # Process SNP data
    results = compare_enzyme_cuts(args.snp_file, enzyme_dict, primer_params)

    # Output the results
    try:
        with open(args.output, mode='w', newline='') as csvfile:
            fieldnames = [
                'Chromosome', 'Position', 'Enzyme',
                'Left Primer', 'Right Primer',
                'Amplicon Size',
                'Amplicon Sequence 1', 'Amplicon Sequence 2',
                'Cuts in Amplicon 1', 'Cuts in Amplicon 2',
                'Positions in Amplicon 1', 'Positions in Amplicon 2',
                'Fragment Sizes Amplicon 1', 'Fragment Sizes Amplicon 2',
                'Unique Cuts',
            ]
            writer = csv.DictWriter(csvfile, fieldnames=fieldnames)

            writer.writeheader()
            for result in results:
                writer.writerow(result)

        print(f"Analysis complete. Results saved to '{args.output}'.")
    except IOError as e:
        print(f"Error writing to output file '{args.output}': {e}")
        sys.exit(1)

if __name__ == "__main__":
    main()
