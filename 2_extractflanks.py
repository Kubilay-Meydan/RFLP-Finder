import os
import csv
import multiprocessing as mp
from functools import partial
from tqdm import tqdm
import argparse

def read_snp_file(snp_file_path):
    snp_data = []
    with open(snp_file_path, 'r') as snp_file:
        for line in snp_file:
            parts = line.strip().split()
            chromosome = parts[0]
            position1 = int(parts[1])
            nucleotide1 = parts[3]
            nucleotide2 = parts[4]
            snp_data.append((chromosome, position1, nucleotide1, nucleotide2))
    print(f"Read {len(snp_data)} SNP entries from {snp_file_path}")
    return snp_data

def get_chromosome_from_filename(fasta_file_name):
    base_name = os.path.basename(fasta_file_name)
    chromosome = base_name.split('.fa')[0].split('.')[-1]
    return chromosome

def extract_sequence_from_chromosome(genome_fasta_folder, snp_info, sequence_size):
    chromosome, position1, nucleotide1, nucleotide2 = snp_info
    fasta_files = [f for f in os.listdir(genome_fasta_folder) if f.endswith('.fa') and get_chromosome_from_filename(f) == chromosome]
    
    if not fasta_files:
        print(f"No FASTA file found for chromosome {chromosome}. Skipping.")
        return None
    
    fasta_file_path = os.path.join(genome_fasta_folder, fasta_files[0])
    
    try:
        with open(fasta_file_path, 'r') as fasta_file:
            sequence = []
            for line in fasta_file:
                if not line.startswith(">"):
                    sequence.append(line.strip())
                    
            sequence = ''.join(sequence)
            start = max(0, position1 - sequence_size - 1)
            end = position1 + sequence_size
            
            if end > len(sequence):
                print(f"Position {position1} is out of bounds for chromosome {chromosome}. Skipping.")
                return None
            
            extracted_sequence = sequence[start:end]
            
            sequence2 = list(extracted_sequence)
            if position1 - start - 1 >= 0 and position1 - start - 1 < len(sequence2):
                sequence2[position1 - start - 1] = nucleotide2
                sequence2 = ''.join(sequence2)
            else:
                sequence2 = extracted_sequence
            
            return (chromosome, position1, nucleotide1, nucleotide2, extracted_sequence, sequence2)
    except Exception as e:
        print(f"Error processing chromosome {chromosome} at position {position1}: {e}")
        return None

def process_snp_data(genome_fasta_folder, snp_data, sequence_size):
    with mp.Pool(mp.cpu_count()) as pool:
        func = partial(extract_sequence_from_chromosome, genome_fasta_folder, sequence_size=sequence_size)
        results = list(tqdm(pool.imap(func, snp_data), total=len(snp_data)))
    return [result for result in results if result is not None]

def write_to_csv(output_file_path, extracted_data):
    with open(output_file_path, 'w', newline='') as csvfile:
        csvwriter = csv.writer(csvfile)
        csvwriter.writerow(['Chromosome', 'Position', 'Nucleotide1', 'Nucleotide2', 'Sequence 1', 'Sequence 2'])
        
        for data in extracted_data:
            csvwriter.writerow(data)
    print(f"Data written to {output_file_path}")

def main():
    parser = argparse.ArgumentParser(description="Extract sequences around SNP positions.")
    parser.add_argument('-gf', '--genome_fasta_folder', required=True, help="Path to the folder containing the genome FASTA files.")
    parser.add_argument('-s', '--snp_file', required=True, help="Path to the SNP file.")
    parser.add_argument('-o', '--output_file', required=True, help="Path to the output CSV file.")
    parser.add_argument('-z', '--sequence_size', type=int, default=501, help="Number of nucleotides before and after the SNP position (default: 501).")
    
    args = parser.parse_args()

    genome_fasta_folder = args.genome_fasta_folder
    snp_file_path = args.snp_file
    output_file_path = args.output_file
    sequence_size = args.sequence_size

    snp_data = read_snp_file(snp_file_path)
    extracted_data = process_snp_data(genome_fasta_folder, snp_data, sequence_size)
    write_to_csv(output_file_path, extracted_data)
    print(f"Extraction complete. Data saved to {output_file_path}.")

if __name__ == "__main__":
    main()
