import os
import ftplib
import gzip
import shutil
from multiprocessing import Pool, cpu_count
import argparse


def download_file(ftp, output_dir, fasta_file):
    """Download a single file from FTP."""
    output_file_gz = os.path.join(output_dir, fasta_file)

    with open(output_file_gz, "wb") as f:
        ftp.retrbinary(f"RETR {fasta_file}", f.write)
    print(f"Successfully downloaded {fasta_file}.")

    return output_file_gz


def download_genome_files(species, assembly):
    """Download genome files for the specified species and assembly."""
    # Setup FTP connection to Ensembl
    ftp_host = "ftp.ensembl.org"
    ftp = ftplib.FTP(ftp_host)
    ftp.login()

    # Format species name to lowercase and underscore format
    species = species.lower().replace(" ", "_")
    species_capitalized = species.capitalize()

    # Directory for storing genomes
    genomes_dir = "GENOMES"
    os.makedirs(genomes_dir, exist_ok=True)

    # Check if genome is already downloaded
    output_dir = os.path.join(
        genomes_dir, f"{species}_{assembly}_release-112_dna_rm_filtered"
    )
    if os.path.exists(output_dir):
        print(
            f"Genome for {species} assembly {assembly} is already "
            f"downloaded in {output_dir}."
        )

        return output_dir

    # Accumulate all available assemblies across releases
    all_assemblies = {}
    found_assembly = False

    # Loop through each release from 19 to 112
    for release in range(112, 18, -1):
        base_path = f"/pub/release-{release}/fasta/{species}/dna/"
        try:
            ftp.cwd(base_path)
            print(f"Checking release {release}...")
        except ftplib.error_perm:
            print(f"Release {release} not found for species '{species}'.")
            continue

        # List all files in the directory
        files = ftp.nlst()

        # List all available assemblies by extracting the text before 'dna'
        assemblies = set()
        for file in files:
            if "dna_rm" in file:
                assembly_version = file.split(".dna_rm")[0].split(
                    species_capitalized + "."
                )[1]
                assemblies.add(assembly_version)
                if assembly_version not in all_assemblies:
                    all_assemblies[assembly_version] = release

        # Check if the specified assembly exists in this release
        if assembly in assemblies:
            print(f"Assembly '{assembly}' found in release {release}.")
            found_assembly = True

            # Define explicit valid file patterns for chromosomes and scaffolds
            valid_patterns = [
                f"{species_capitalized}.{assembly}"
                f".dna_rm.chromosome.{chrom}.fa.gz"
                for chrom in list(map(str, range(1, 23))) + ["X", "Y", "MT"]
            ]
            valid_patterns.append(
                f"{species_capitalized}.{assembly}.dna_rm.scaffold.fa.gz"
            )

            # Filter files that match the valid patterns
            fasta_files = [f for f in files if f in valid_patterns]

            # If no chromosome files found, check for primary assembly files
            if not fasta_files:
                print(
                    f"No 'chromosome' files found for assembly '{assembly}'."
                    f" Checking for 'primary_assembly' files..."
                )
                valid_patterns = [
                    f"{species_capitalized}.{assembly}"
                    ".dna_rm.primary_assembly.{chrom}.fa.gz"
                    for chrom in list(map(str, range(1, 23))) + ["X", "Y", "MT"]
                ]
                fasta_files = [f for f in files if f in valid_patterns]

                if not fasta_files:
                    print(
                        f"No 'primary_assembly' files found for assembly '{assembly}' for species '{species}' in release {release}."
                    )
                    continue

            # Create directory for the species and assembly
            output_dir = os.path.join(
                genomes_dir, f"{species}_{assembly}_release-{release}_dna_rm_filtered"
            )
            os.makedirs(output_dir, exist_ok=True)

            # Download all files sequentially
            for fasta_file in fasta_files:
                print(f"Downloading {fasta_file} from release {release}...")
                download_file(ftp, output_dir, fasta_file)

            ftp.quit()
            print(f"All files from release {release} downloaded successfully.")

            # Return the directory where the files were saved
            return output_dir

        else:
            print(f"Assembly '{assembly}' not found in release {release}.")

    # If the specified assembly was not found in any release
    if not found_assembly:
        print(
            f"Assembly '{assembly}' was not found for species '{species}' in any release from 19 to 112."
        )
        print("Available assemblies for this species across all releases are:")
        for assembly_version in sorted(all_assemblies, key=lambda x: all_assemblies[x]):
            if all_assemblies[assembly_version] == 112:
                print(f" - {assembly_version} (current)")
            else:
                print(f" - {assembly_version}")

    ftp.quit()
    return None


def decompress_genome_files(output_dir):
    """Decompress all .gz files in the specified directory."""
    # List all .gz files in the output directory
    gz_files = [
        os.path.join(output_dir, f) for f in os.listdir(output_dir) if f.endswith(".gz")
    ]

    # Unzip files in parallel using multiprocessing
    print(f"Unzipping files using {cpu_count()} cores...")
    with Pool(cpu_count()) as pool:
        pool.map(unzip_file, gz_files)

    print("All files unzipped successfully.")


def unzip_file(file_path):
    """Unzip a single .gz file and remove the original .gz file."""
    output_file = file_path[:-3]  # Remove the .gz extension

    with gzip.open(file_path, "rb") as f_in:
        with open(output_file, "wb") as f_out:
            shutil.copyfileobj(f_in, f_out)
    os.remove(file_path)
    print(f"Successfully unzipped {file_path}.")

    return output_file


def main():
    parser = argparse.ArgumentParser(
        description="Download and process genome files from Ensembl FTP server."
    )
    parser.add_argument(
        "-s",
        "--species",
        type=str,
        required=True,
        help="Species name (e.g., 'Homo sapiens').",
    )
    parser.add_argument(
        "-a",
        "--assembly",
        type=str,
        required=True,
        help="Assembly name (e.g., 'ARS-UCD1.3').",
    )
    args = parser.parse_args()

    # Step 1: Download the files
    output_dir = download_genome_files(args.species, args.assembly)

    # Step 2: Decompress the files
    if output_dir:
        decompress_genome_files(output_dir)
        print(
            f"FASTA folder: {output_dir}"
        )  # Explicitly print for the main script to capture


if __name__ == "__main__":
    main()
