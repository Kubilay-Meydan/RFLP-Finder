import os
import subprocess
import json

def get_input(prompt_text, default_value=None):
    """Prompt the user for input, with an optional default value."""
    if default_value:
        prompt_text += f" [{default_value}]: "
    else:
        prompt_text += ": "
    user_input = input(prompt_text)
    return user_input if user_input else default_value

# Ensure RESULTS directory exists
os.makedirs("RESULTS", exist_ok=True)

# Ask if the user wants to use autopilot mode
autopilot_mode = get_input("Do you want to use autopilot mode? (yes/no)", "no").strip().lower() == "yes"

if autopilot_mode:
    config_file_path = get_input("Enter the path to the config file", "config.json")
    try:
        with open(config_file_path, "r") as config_file:
            config = json.load(config_file)
        print("Configuration loaded successfully.")
    except Exception as e:
        print(f"Error loading config file: {e}")
        exit(1)

    use_provided_genome = config.get("use_provided_genome", False)
    if use_provided_genome:
        fasta_folder = config.get("fasta_folder")
        if not fasta_folder:
            print("Error: 'fasta_folder' must be specified in the config file.")
            exit(1)
    else:
        species = config.get("species")
        assembly = config.get("assembly")
        if not species or not assembly:
            print("Error: 'species' and 'assembly' must be specified in the config file.")
            exit(1)

        print("=== Step 1: Download Genome Files ===")
        result = subprocess.run(
            ["python", "1_get_genome.py", "-s", species, "-a", assembly],
            stdout=subprocess.PIPE,  # Capture stdout
            stderr=subprocess.STDOUT,  # Redirect stderr to stdout
            text=True
        )

        # Print all captured output
        print(result.stdout)

        # Capture the FASTA folder from the output of Step 1
        fasta_folder = None
        for line in result.stdout.splitlines():
            if line.startswith("FASTA folder:"):
                fasta_folder = line.split(":")[1].strip()
                break

        if not fasta_folder:
            print("Error: Could not determine the FASTA folder from Step 1. Exiting.")
            exit(1)

        print(f"FASTA folder determined as: {fasta_folder}")

    # Collect additional inputs from config
    snp_file = config.get("snp_file")
    sequence_size = config.get("sequence_size", "501")
    input_file_enzyme = config.get("input_file_enzyme")
    output_file_rflp_name = config.get("output_file_rflp_name", "output.csv")
    output_file_rflp = os.path.join("RESULTS", output_file_rflp_name)
    primer_opt_size = config.get("primer_opt_size", "20")
    primer_min_size = config.get("primer_min_size", "18")
    primer_max_size = config.get("primer_max_size", "25")
    primer_opt_tm = config.get("primer_opt_tm", "60.0")
    primer_min_tm = config.get("primer_min_tm", "57.0")
    primer_max_tm = config.get("primer_max_tm", "63.0")
    primer_product_size_range = config.get("primer_product_size_range", "100-300")
else:
    # Ask the user if they are providing the genome
    use_provided_genome = get_input("Do you provide the genome FASTA folder? (yes/no)", "no").strip().lower() == "yes"

    if use_provided_genome:
        # If the user provides the genome, ask for the FASTA folder path
        fasta_folder = get_input("Enter the path to the genome FASTA folder (repeatmasked)")
    else:
        # Collect inputs for genome download
        print("=== Collecting Inputs for Genome Download ===")
        species = get_input("Enter species name (e.g., Homo sapiens)")
        assembly = get_input("Enter assembly name (e.g., ARS-UCD1.3)")

        print("=== Step 1: Download Genome Files ===")
        result = subprocess.run(
            ["python", "1_get_genome.py", "-s", species, "-a", assembly],
            stdout=subprocess.PIPE,  # Capture stdout
            stderr=subprocess.STDOUT,  # Redirect stderr to stdout
            text=True
        )

        # Print all captured output
        print(result.stdout)

        # Capture the FASTA folder from the output of Step 1
        fasta_folder = None
        for line in result.stdout.splitlines():
            if line.startswith("FASTA folder:"):
                fasta_folder = line.split(":")[1].strip()
                break

        if not fasta_folder:
            print("Error: Could not determine the FASTA folder from Step 1. Exiting.")
            exit(1)

        print(f"FASTA folder determined as: {fasta_folder}")

    # Collect additional inputs
    print("=== Collecting Inputs for Processing ===")
    snp_file = get_input("Enter path to the SNP file")
    sequence_size = get_input("Enter sequence size (default 501)", "501")
    input_file_enzyme = get_input("Enter path to enzyme recognition pattern input file")
    output_file_rflp_name = get_input("Name your final analysis file name", "(.csv)")
    output_file_rflp = os.path.join("RESULTS", output_file_rflp_name)
    primer_opt_size = get_input("Enter optimal primer size (default 20)", "20")
    primer_min_size = get_input("Enter minimum primer size (default 18)", "18")
    primer_max_size = get_input("Enter maximum primer size (default 25)", "25")
    primer_opt_tm = get_input("Enter optimal primer melting temperature (Tm) (default 60.0)", "60.0")
    primer_min_tm = get_input("Enter minimum primer Tm (default 57.0)", "57.0")
    primer_max_tm = get_input("Enter maximum primer Tm (default 63.0)", "63.0")
    primer_product_size_range = get_input("Enter desired PCR product size range (e.g., 100-300)", "100-300")

# Execute the steps
print("=== Step 2: Extract Flanking Sequences ===")
subprocess.run([
    "python", "2_extractflanks.py", "-gf", fasta_folder, "--snp_file", snp_file,
    "--output_file", "WORK/FLANKS", "--sequence_size", sequence_size
])

print("=== Step 3: Process Enzyme Recognition Sites ===")
subprocess.run([
    "python", "3_recognitionpatterndecliner.py", "-input_file", input_file_enzyme, "-output_file", "WORK/ENZYME_PATTERNS"
])

print("=== Step 4: Perform RFLP Analysis ===")
subprocess.run([
    "python", "4_cutter.py", "WORK/FLANKS", "WORK/ENZYME_PATTERNS", "--output", output_file_rflp,
    "--PRIMER_OPT_SIZE", primer_opt_size, "--PRIMER_MIN_SIZE", primer_min_size, "--PRIMER_MAX_SIZE", primer_max_size,
    "--PRIMER_OPT_TM", primer_opt_tm, "--PRIMER_MIN_TM", primer_min_tm, "--PRIMER_MAX_TM", primer_max_tm,
    "--PRIMER_PRODUCT_SIZE_RANGE", primer_product_size_range
])

print("All steps completed! Final analysis file saved in:", output_file_rflp)

# Example config file (config.json):
# {
#     "use_provided_genome": true,
#     "fasta_folder": "/path/to/repeatmasked/folder",
#     "snp_file": "/path/to/snp_file",
#     "sequence_size": "501",
#     "input_file_enzyme": "/path/to/enzyme/file",
#     "output_file_rflp_name": "final_output.csv",
#     "primer_opt_size": "20",
#     "primer_min_size": "18",
#     "primer_max_size": "25",
#     "primer_opt_tm": "60.0",
#     "primer_min_tm": "57.0",
#     "primer_max_tm": "63.0",
#     "primer_product_size_range": "100-300"
# }
