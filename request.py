def main():
    decisions = {}

    # Helper functions
    def get_input(prompt):
        return input(prompt)
    
    def store_value(key, value):
        decisions[key] = value

    # Start of the decision tree
    enzyme_list = get_input("Do you want to input enzyme list? (yes/no): ").strip().lower()
    store_value("Input enzyme list?", enzyme_list)
    
    if enzyme_list == 'yes':
        enzyme_list_file_path = get_input("Enter enzyme list file path: ").strip()
        store_value("Enzyme list file path", enzyme_list_file_path)
        
        max_restriction_size = get_input("Enter max size of restriction site: ").strip()
        store_value("Max size of Restriction site?", max_restriction_size)
    else:
        max_restriction_size = get_input("Enter max size of restriction site: ").strip()
        store_value("Max size of Restriction site?", max_restriction_size)
        
        enzyme_list_from_memory = "Use enzyme list from memory"
        store_value("Enzyme list file path", enzyme_list_from_memory)

    snp_input = get_input("Do you want to input SNP? (yes/no): ").strip().lower()
    store_value("Input SNP?", snp_input)
    
    if snp_input == 'yes':
        snp_list_file_path = get_input("Enter SNP list file path: ").strip()
        store_value("SNP list file path", snp_list_file_path)
    
    genome_input = get_input("Do you want to input Genome? (yes/no): ").strip().lower()
    store_value("Input Genome?", genome_input)
    
    if genome_input == 'yes':
        genome_file_path = get_input("Enter Genome file path: ").strip()
        store_value("Genome file path", genome_file_path)
        
        if snp_input == 'no':
            species_name_assembly = get_input("Enter Species Name and assembly (ENSEMBL): ").strip()
            store_value("Species Name and assembly? (ENSEMBL)", species_name_assembly)
            store_value("SNP list", f"Download SNPs From ENSEMBL for {species_name_assembly}")
    else:
        species_name_assembly = get_input("Enter Species Name and assembly (ENSEMBL): ").strip()
        store_value("Species Name and assembly? (ENSEMBL)", species_name_assembly)
        
        if snp_input == 'yes':
            store_value("SNP list file path", snp_list_file_path)  # Keeps the file path input for SNP
            store_value("Download Genome Fasta & SNPs From ENSEMBL", f"Download Genome Fasta & SNPs From ENSEMBL for {species_name_assembly}")
        else:
            store_value("Download Genome Fasta From ENSEMBL", f"Download Genome Fasta From ENSEMBL for {species_name_assembly}")
            store_value("SNP list", f"Download SNPs From ENSEMBL for {species_name_assembly}")
    
    max_amplicon_size = get_input("Enter max Amplicon size: ").strip()
    store_value("Max Amplicon size?", max_amplicon_size)
    
    species_for_masking_repeats = get_input("Enter Species to use for Masking Repeats: ").strip()
    store_value("Species to use for Masking Repeats?", species_for_masking_repeats)
    
    primer3_settings = get_input("Enter Primer3 Settings: ").strip()
    store_value("Primer3 Settings", primer3_settings)
    
    request_formulated = "Request formulated"
    store_value("Request formulated", request_formulated)
    
    return decisions

if __name__ == "__main__":
    decisions_made = main()
    for key, value in decisions_made.items():
        print(f"{key}: {value}")
