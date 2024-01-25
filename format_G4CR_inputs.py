# AUTHOR: Lucas Nelson
import os
import sys

ANCESTRAL_GENOMES = ['hg38', '_HP', '_HPG', '_HPGP', '_HPGPN', '_HPGPNRMPC', '_HPGPNRMPCCS', '_HPGPNRMPCCSO', '_HPGPNRMPCCSOT', '_HPGPNRMPCCSOTSJMCMMRHCCOOO', '_HPGPNRMPCCSOTSJMCMMRHCCOOOSVCTOPBOCECFCMAOLPPEMMESC', '_HPGPNRMPCCSOTSJMCMMRHCCOOOSVCTOPBOCECFCMAOLPPEMMESCLETCEOD']
ANCESTRAL_GENOMES.reverse()
NTS = ('A','C','G','T')
MAX_LOOP = 7
MAX_BULGE = 3
MIN_TEMP = 50

def get_subdirectories(directory):
    subdirectories = [d for d in os.listdir(directory) if os.path.isdir(os.path.join(directory, d))]
    return subdirectories

def format_input_files(chrom_filepath, output_filepath):
    '''
    read_input_files(chrom_filepath: str) -> 
        gene_names: List[str],
        num_G4CRs_per_gene: List[int],
        ancestral_G4CR_regions: List[str], 
        start_and_end_indices:  List[tuple(str, str)], 
        real_match_rates: List[float], 
        ancestral_seq_names: List[str]
        
    len(gene_names) = len(num_G4CRs_per_gene) # gene-indexed
    len(ancestral_G4CR_regions) = len(start_and_end_indices) = len(real_match_rates) = len(ancestral_seq_names) # G4CR-indexed
    
    Given a chromosome number and filepath to the input directory for that chromosome,
    Reads the contents of the directory and returns the data for each gene/G4CR (depending on output) in list format
    ''' 
    
    chroms = get_subdirectories(chrom_filepath)
    for chrom in chroms:
        input_filepath = f"{chrom_filepath}/{chrom}"
        for input_file in os.listdir(input_filepath):
            # file_name = input_filepath.split("/")[-1]
            gene_name = input_file.split(".")[0]

            with open(f"{input_filepath}/{input_file}", "r") as input_file:
                
                # Extract data for all the G4CRs and ancestral sequences + match rates in the given gene's input file
                
                modern_G4CRs = []
                ancestral_aligned_sequences = {}
                ancestral_match_rates = {}
                ancestral_polymorphism_vals = {}
                G4CR_length_in_NTs = {}
                g_percentages = {}
                ntot = {}
                ntand = {}
                tm_max = {}
                tm_median = {}
                tm_min = {}
                
                input_lines = input_file.readlines()
                
                for line in input_lines:
                    if line == "":
                        continue
                    line = line.split(",")
                    
                    if line[0] == "hg38" and line[3] != '' and line[4] != '': # Not interested in human promoters (files) containing no G4CRs
                        modern_G4CRs.append(line)
                        
                    if len(line) == 14:
                        seq_name = line[0]
                        poly_vals = line[2]
                        this_G4CR_length = line[5]
                        g_percentage = line[6]
                        this_ntot = line[7]
                        this_ntand = line[8]
                        this_tm_max = line[9]
                        this_tm_median = line[10]
                        this_tm_min = line[11]
                        match_rate = line[12]
                        aligned_sequence = line[13]
                        
                        ancestral_polymorphism_vals[seq_name] = poly_vals
                        G4CR_length_in_NTs[seq_name] = this_G4CR_length
                        g_percentages[seq_name] = g_percentage
                        ntot[seq_name] = this_ntot
                        ntand[seq_name] = this_ntand
                        tm_max[seq_name] = this_tm_max
                        tm_median[seq_name] = this_tm_median
                        tm_min[seq_name] = this_tm_min
                        ancestral_aligned_sequences[seq_name] = aligned_sequence
                        ancestral_match_rates[seq_name] = match_rate
                
                # Iterate through each line of information for the given gene
                for G4CR_info in modern_G4CRs:
                    start_index = int(G4CR_info[3])
                    end_index = int(G4CR_info[4])

                    # Expand the window being simulated to -20 and +20 beyond the start and end indices
                    if start_index >= 20 and end_index < 3979: 
                        # 4000 - 21 = 3979, i.e., 20 positions before the end of the ancestral sequence
                        start_index -= 20
                        end_index += 21
                    elif end_index < 3979:
                        start_index = 0
                        end_index += 21
                    elif start_index >= 20:
                        start_index -= 20
                        end_index = 3999
                    else:
                        start_index = 0
                        end_index = 3999

                    ancestral_promoter_name = ""
                    ancestral_promoter = ""
                    index = 0
                    
                    included_ancestral_genomes = {}
                    # Check that there are enough actual nucleotides in the REGION, not the WHOLE PROMOTER
                    while index < len(ANCESTRAL_GENOMES):
                        G4CR_seq = []
                        ancestral_promoter_name = ANCESTRAL_GENOMES[index]
                        try:
                            ancestral_promoter = ancestral_aligned_sequences[ancestral_promoter_name]
                        except KeyError: # Genome in question doesn't exist for this file
                            index += 1
                            continue
                        for i in range(start_index, end_index + 1):
                            if ancestral_promoter[i] in NTS:
                                G4CR_seq.append(ancestral_promoter[i])
                        if len(G4CR_seq) >= 15:
                            included_ancestral_genomes[ancestral_promoter_name] = "".join(G4CR_seq)
                        index += 1
                    
                    if ancestral_promoter == "": # RARE CASE: The modern G4CR was inserted after the human-chimp ancestor (_HP)
                        continue
                    
                    this_output_filepath = f"{output_filepath}/{chrom}/{gene_name}_{start_index}_{end_index}.csv"
                    with open(this_output_filepath, "w") as output_file:
                        for ancestral_promoter_name in sorted(included_ancestral_genomes.keys(), reverse=True, key=lambda x: '' if x=='hg38' else x):
                            real_match_rate = ancestral_match_rates[ancestral_promoter_name]
                            poly_vals = ancestral_polymorphism_vals[ancestral_promoter_name]
                            this_G4CR_length = G4CR_length_in_NTs[ancestral_promoter_name]
                            g_percentage = g_percentages[ancestral_promoter_name] 
                            this_ntot = ntot[ancestral_promoter_name] 
                            this_ntand = ntand[ancestral_promoter_name]
                            this_tm_max = tm_max[ancestral_promoter_name]
                            this_tm_median = tm_median[ancestral_promoter_name]
                            this_tm_min = tm_min[ancestral_promoter_name]
                            G4CR_seq = included_ancestral_genomes[ancestral_promoter_name]
                            output_file.write(f"{ancestral_promoter_name},{real_match_rate},{poly_vals},{this_G4CR_length},{g_percentage},{this_ntot},{this_ntand},{this_tm_max},{this_tm_median},{this_tm_min},{G4CR_seq}\n")
    

def main():
    input_filepath = "/home/mcb/users/lnelso12/evoGReg/outputs"
    output_filepath = "/home/mcb/users/lnelso12/evoGReg/outputs_formatted"
    format_input_files(input_filepath, output_filepath)

main()
