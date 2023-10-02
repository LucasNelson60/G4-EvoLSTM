import sys
import os
import random
import numpy as np
from simulate import simulate
from greg import main as greg_main
from get_GQs import get_GQs

ANCESTRAL_GENOMES = ['_HP', '_HPG', '_HPGP', '_HPGPN', '_HPGPNRMPC', '_HPGPNRMPCCS', '_HPGPNRMPCCSO', '_HPGPNRMPCCSOT', '_HPGPNRMPCCSOTSJMCMMRHCCOOO', '_HPGPNRMPCCSOTSJMCMMRHCCOOOSVCTOPBOCECFCMAOLPPEMMESC', '_HPGPNRMPCCSOTSJMCMMRHCCOOOSVCTOPBOCECFCMAOLPPEMMESCLETCEOD']
ANCESTRAL_GENOMES.reverse()
NTS = ('A','C','G','T')
MAX_LOOP = 7
MAX_BULGE = 3
MIN_TEMP = 50

def DNA(val): # function that takes in a number from 1 to 3 as input, outputs a corresponding insertion of random nts
    nts = ""
    for _ in range(val):
        val = random.random()
        if val < 0.25:
            nts += "A"
        elif val < 0.5:
            nts += "C"
        elif val < 0.75:
            nts += "G"
        else:
            nts += "T"
    return nts

def read_input_files(chrom, chrom_filepath):
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
    # One-to-one correspondence between the two following lists (same lengths)
    start_and_end_indices = []
    ancestral_G4CR_regions = []
    real_match_rates = []
    ancestral_seq_names = []
    
    # One-to-one correspondence between the two following lists (same lengths)
    num_G4CRs_per_gene = []
    gene_names = []    

    for input_file in os.listdir(chrom_filepath):
        
        input_filepath = f"{chrom_filepath}/{input_file}"
        file_name = input_filepath.split("/")[-1]
        gene_name = file_name.split(".")[0]
        G4CR_counter = 0

        with open(input_filepath, "r") as input_file:
            
            # Extract data for all the G4CRs and ancestral sequences + match rates in the given gene's input file
            
            modern_G4CRs = []
            ancestral_aligned_sequences = {}
            ancestral_match_rates = {}
            
            input_lines = input_file.readlines()
            
            for line in input_lines:
                if line == "":
                    continue
                line = line.split(",")
                
                if line[0] == "hg38" and line[3] != '' and line[4] != '': # Not interested in human promoters (files) containing no G4CRs
                    modern_G4CRs.append(line)
                    G4CR_counter += 1
                    
                elif len(line) == 14:
                    seq_name = line[0]
                    match_rate = line[12]
                    aligned_sequence = line[13]
                    
                    ancestral_aligned_sequences[seq_name] = aligned_sequence
                    ancestral_match_rates[seq_name] = match_rate
            
            # GENE-INDEXED DATA
            num_G4CRs_per_gene.append(G4CR_counter) # Will be important later when iterating through every G4CR of every gene
            gene_names.append(gene_name)
            
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
                
                # Check that there are enough actual nucleotides in the REGION, not the WHOLE PROMOTER
                while index < len(ANCESTRAL_GENOMES):
                    ancestral_promoter_name = ANCESTRAL_GENOMES[index]
                    try:
                        ancestral_promoter = ancestral_aligned_sequences[ancestral_promoter_name]
                    except KeyError: # Genome in question doesn't exist for this file
                        index += 1
                        continue
                    counter = 0
                    for i in range(start_index, end_index + 1):
                        if ancestral_promoter[i] in NTS:
                            counter += 1
                    if counter < 15:
                        index += 1
                    else:
                        break
                
                if ancestral_promoter == "": # RARE CASE: The modern G4CR was inserted after the human-chimp ancestor (_HP)
                    with open(f"/home/mcb/users/lnelso12/G4_evoLSTM/discarded_G4CRs/chr{chrom}.csv", "a") as discard_file:
                        discard_file.write(f"{gene_name},{start_index},{end_index},{G4CR_info[1]},recent G4CR insertion\n")
                    continue # These will be set aside and skipped to be studied individually in closer detail
                
                # G4CR-INDEXED DATA
                ancestral_G4CR_region = "".join([ancestral_promoter[i] for i in range(start_index, end_index) if ancestral_promoter[i] in NTS])
                
                ancestral_G4CR_regions.append(ancestral_G4CR_region) 
                start_and_end_indices.append((start_index,end_index))
                real_match_rates.append(match_rate)
                ancestral_seq_names.append(ancestral_promoter_name)
    
    G4CR_indexed_gene_names = [] # gene_names will now be G4CR-INDEXED to standardize the lengths of every list
    for i in range(len(num_G4CRs_per_gene)):
        gene_name = gene_names[i]
        num_G4CRs_for_this_gene = num_G4CRs_per_gene[i]
        G4CR_indexed_gene_names.extend([gene_name for _ in range(num_G4CRs_for_this_gene)])

    return G4CR_indexed_gene_names, ancestral_G4CR_regions, start_and_end_indices, real_match_rates, ancestral_seq_names

def align_simulated_outputs_with_ancestors(aligned_ancestral_G4CRs, sim_output_seqs):
    '''
    align_simulated_outputs(sim_input_seqs: str, sim_output_seqs: str) -> aligned_input_seqs: str, aligned_output_seqs: str
    
    Given the input/output simulations (in string format) of simulate.py, return the aligned + expanded input/output seqs
    '''
    
    # BELOW: Maintain delimiters, replace 1,2,3 insertion chars, re-align input/output seqs
    new_ancestor = []
    new_output = []
    index = 0

    while index < len(sim_output_seqs)-1: # Clip by one because EvoLSTM adds a redundant nt at the end of each output

        output_char = sim_output_seqs[index]
        ancestor_char = aligned_ancestral_G4CRs[index]
        new_ancestor.append(ancestor_char)
        
        # CRUCIAL CHECK: delimiters must remain aligned between parallel simulation outputs
        if ancestor_char == "3": # Need to add a delimiter, regardless of what was present in simulation output
            new_output.append("3")
        elif output_char == "1":
            new_ancestor.append("-")
            new_output.append(ancestor_char)
            new_output.append(DNA(1))
        elif output_char == "2":
            new_ancestor.append("--")
            new_output.append(ancestor_char)
            new_output.append(DNA(2))
        elif output_char == "3" and ancestor_char != "3": # Not present in ancestor, hence not a delimiter
            z = np.random.geometric(p=0.1)
            gap_len = z + 2
            new_ancestor.append("-"*gap_len)
            new_output.append(ancestor_char)
            new_output.append(DNA(gap_len))
        else:
            new_output.append(output_char)

        index += 1

    aligned_ancestral_G4CRs = "".join(new_ancestor)
    aligned_output_seqs = "".join(new_output)
    
    return aligned_ancestral_G4CRs, aligned_output_seqs
    
def get_new_match_rates(sim_input_seqs, sim_output_seqs):
    '''
    get_new_match_rates(sim_input_seqs: List[str], sim_output_seqs: List[str]) -> List[float]
    IMPORTANT: len(sim_input_seqs) = len(sim_output_seqs)
    
    Given the simulation inputs and outputs, in list format (split by "3333333333"),
    Returns the match rate of each output sim compared to its corresponding input
    '''
    new_match_rates = []
    
    for i in range(sim_output_seqs):
        
        sim_input_seq = sim_input_seqs[i]
        sim_output_seq = sim_output_seqs[i]
        
        # Get the match rate between the input and output seqs
        ungapped_length = 0
        k = 0
        matches = 0
        
        while k < len(sim_input_seq):
            if sim_output_seq[k] not in NTS:
                k += 1
                continue
            elif sim_output_seq[k] == sim_input_seq[k]:
                matches += 1
            ungapped_length += 1
            k += 1
            
        # Skip simulations that have zero ungapped_length, i.e., that consist of all "-" characters, these are useless
        if ungapped_length == 0:
            continue
        
        new_match_rate = matches/ungapped_length
        new_match_rates.append(new_match_rate)
    
    return new_match_rates

def write_G4CR_to_output_file(chrom, gene_name, ancestral_seq_name, start_index, end_index, sim_input_seqs_for_this_G4CR, aligned_ancestors_for_this_G4CR, GQs, scores, real_match_rate, num_recursions):
    
    output_filepath = f"/home/mcb/users/lnelso12/G4_evoLSTM/outputs/chr{chrom}/{gene_name}_{start_index}_{end_index}.csv"
    with open(output_filepath, "w") as output_file:
        
        # Headers outlining the REGIONAL, ISOMER, and UNIVERSAL fields
        output_file.write("SIMULATION FIELDS,num_g4crs,sim_match_rate,aligned_input_seq,aligned_output_seq,max_ntot,sum_ntot,max_ntand,sum_ntands,max_tm_max,max_tm_median,max_tm_min,sum_g4cr_length,max_g4cr_length\n")
        output_file.write("G4CR FIELDS,g4cr_seq,nt_polymorphism,start_index,end_index,g4cr_length,g_percentage,ntot,ntand,tm_max,tm_median,tm_min\n")
        output_file.write("UNIVERSAL FIELDS,ancestral_seq_name,real_match_rate,num_recursions,avg_num_g4crs,avg_sim_match_rate,avg_max_ntot,avg_sum_ntots,avg_max_ntand,avg_sum_ntands,avg_max_tm_max,avg_max_tm_median,avg_max_tm_min,avg_sum_g4cr_length,avg_max_g4cr_length,avg_g_percentage\n")
        output_file.write("\n")
        
        sum_num_G4CRs = 0
        sum_sim_match_rate = 0.0
        sum_max_ntot = 0
        sum_sum_ntots = 0
        sum_max_ntand = 0
        sum_sum_ntands = 0
        sum_max_Tm_max = 0
        sum_max_Tm_median = 0
        sum_max_Tm_min = 0
        sum_sum_G4CR_length = 0
        sum_max_G4CR_length = 0
        sum_avg_g_percentage = 0.0
        sum_max_g_percentage = 0
        
        for z in range(100):
            
            sim_input_seq = sim_input_seqs_for_this_G4CR[z] # Use input, not output, because its branch length is closer to the real one
            aligned_ancestor = aligned_ancestors_for_this_G4CR[z]
            
            ungapped_length = 0
            k = 0
            matches = 0
            while k < len(aligned_ancestor):
                if sim_input_seq[k] not in NTS:
                    k += 1
                    continue
                elif sim_input_seq[k] == aligned_ancestor[k]:
                    matches += 1
                ungapped_length += 1
                k += 1
            
            # Skip simulations that have zero ungapped_length, i.e., that consist of all "-" characters, to avoid ZeroDivisionError below
            if ungapped_length == 0:
                continue

            sim_match_rate = matches/ungapped_length
            
            # Run GReg on the simulation output
            gapless_output = "".join([nt for nt in sim_input_seq if nt != "-"])
            G4CR_fields, num_G4CRs, max_ntot, sum_ntots, max_ntand, sum_ntands, max_Tm_max, max_Tm_median, max_Tm_min, sum_G4CR_length, max_G4CR_length, avg_g_percentage, max_g_percentage = greg_main(gapless_output, MAX_LOOP, GQs, scores)
            
            output_file.write(f"SIMULATION,{num_G4CRs},{sim_match_rate},{aligned_ancestor},{sim_input_seq},{max_ntot},{sum_ntots},{max_ntand},{sum_ntands},{max_Tm_max},{max_Tm_median},{max_Tm_min},{sum_G4CR_length},{max_G4CR_length}\n")
            output_file.write(G4CR_fields)
            
            sum_num_G4CRs += int(num_G4CRs)
            sum_sim_match_rate += float(sim_match_rate)
            sum_max_ntot += int(max_ntot)
            sum_sum_ntots += int(sum_ntots)
            sum_max_ntand += int(max_ntand)
            sum_sum_ntands += int(sum_ntands)
            sum_max_Tm_max += int(max_Tm_max)
            sum_max_Tm_median += int(max_Tm_median)
            sum_max_Tm_min += int(max_Tm_min)
            sum_sum_G4CR_length += int(sum_G4CR_length)
            sum_max_G4CR_length += int(max_G4CR_length)
            sum_avg_g_percentage += float(avg_g_percentage)
            sum_max_g_percentage += int(max_g_percentage)
        
        output_file.write(f"UNIVERSAL,{ancestral_seq_name},{real_match_rate},{num_recursions},{sum_num_G4CRs/100},{sum_sim_match_rate/100},{sum_max_ntot/100},{sum_sum_ntots/100},{sum_max_ntand/100},{sum_sum_ntands/100},{sum_max_Tm_max/100},{sum_max_Tm_median/100},{sum_max_Tm_min/100},{sum_sum_G4CR_length/100},{sum_max_G4CR_length/100},{sum_avg_g_percentage/100},{sum_max_g_percentage/100}\n")

def desimulate_output_seq(aligned_ancestors_for_this_G4CR, sim_output_seqs_for_this_G4CR, desimulating_distance):
    
    for k in range(100): # iterate through each of the simulations
        
        # Split the lists to be able to replace mutations as needed
        aligned_ancestor = aligned_ancestors_for_this_G4CR[k]
        sim_output = sim_output_seqs_for_this_G4CR[k]
        
        new_aligned_ancestor = []
        new_sim_output = []
        
        reversing_indel = False
        
        for j in range(len(aligned_ancestor)):
            
            ancestor_char = aligned_ancestor[j]
            sim_char = sim_output[j]
            
            if ancestor_char == "-" or sim_char == "-": # insertion/deletion (indel)
                if reversing_indel:
                    continue
                elif random.random() <= desimulating_distance:
                        reversing_indel = True
                else:
                    new_aligned_ancestor.append(ancestor_char)
                    new_sim_output.append(sim_char)
                    
            elif ancestor_char != sim_char: # substitution
                if random.random() <= desimulating_distance:
                    new_sim_output.append(ancestor_char)
                else:
                    new_sim_output.append(sim_output)
                new_aligned_ancestor.append(ancestor_char)
                
            elif ancestor_char == sim_char:
                new_aligned_ancestor.append(ancestor_char)
                new_sim_output.append(sim_char)
                
                if reversing_indel:
                    reversing_indel = False
        
        # Replace simulations for this G4CR after de-simulating them by desired branch length
        aligned_ancestors_for_this_G4CR[k] = new_aligned_ancestor
        sim_output_seqs_for_this_G4CR[k] = new_sim_output

def main():
    '''
    BASIC IDEA -> start by assembling the input string of 100 * each G4CR input (from files)
    1. Run a single simulation
    2. Split output string, calculate whether each output G4CR is closer or further from real_match_rate than input
    3. Remove those that are FURTHER, run GReg on these, and output to files
    4. Build new string from those that are CLOSER, and repeat process from step 1 until string is empty
    '''
    
    chrom = sys.argv[1]
    gpu = sys.argv[2]
    
    GQs, scores = get_GQs(MAX_LOOP, MAX_BULGE, MIN_TEMP)
    
    chrom_filepath = f"/home/mcb/users/lnelso12/evoGReg/outputs/chr{chrom}"
    
    # NOTE THAT EVERY OUTPUT OF THIS FUNCTION IS G4CR-INDEXED AND OF THE SAME LENGTH
    gene_names, aligned_ancestral_G4CRs, start_and_end_indices, real_match_rates, ancestral_seq_names = read_input_files(chrom, chrom_filepath)
    print("Successfully read all input files and stored their information to begin simulations")
    
    sim_input_seqs = [] # Will be passed in as input to simulation across every recursion
    for i in range(len(aligned_ancestral_G4CRs)): # 100 simulation inputs for each G4CR
        sim_input_seqs.extend([aligned_ancestral_G4CRs[i] for _ in range(100)])
        
    aligned_ancestral_G4CRs = sim_input_seqs.copy() # Will be progressively aligned to simulation outputs at each recursion 
    sim_input_seqs = "3333333333".join(sim_input_seqs) # length = num_G4CRs * 100
    
    num_recursions = 0
    while sim_input_seqs != "":
        
        sim_output_seqs = simulate(gpu, sim_input_seqs) # TEST SIMULATION CALL
        aligned_ancestral_G4CRs, sim_output_seqs = align_simulated_outputs_with_ancestors("3333333333".join(aligned_ancestral_G4CRs), sim_output_seqs) # Always relative to prev sim
        
        aligned_ancestral_G4CRs = "3333333333".split(aligned_ancestral_G4CRs)
        sim_output_seqs = "3333333333".split(sim_output_seqs)
        new_match_rates = get_new_match_rates(aligned_ancestral_G4CRs, sim_output_seqs) # Always relative to real ancestors

        new_sim_output_seqs = []
        new_aligned_ancestral_G4CRs = []
        
        new_gene_names = []
        new_start_and_end_indices = []
        new_real_match_rates = []
        new_ancestral_seq_names = []
        
        for i in range(100, len(new_match_rates), 100):

            new_match_rates_for_this_G4CR = [new_match_rates[j] for j in range(i-100, i)]

            avg_new_match_rate = sum(new_match_rates_for_this_G4CR)/100
            real_match_rate = real_match_rates[int(i/100)]

            sim_output_seqs_for_this_G4CR = [sim_output_seqs[l] for l in range(i-100,i)]
            aligned_ancestors_for_this_G4CR = [aligned_ancestral_G4CRs[b] for b in range(i-100,i)]
            
            gene_name = gene_names[int(i/100)]
            ancestral_seq_name = ancestral_seq_names[int(i/100)]
            start_index = start_and_end_indices[int(i/100)][0]
            end_index = start_and_end_indices[int(i/100)][1]

            if real_match_rate == avg_new_match_rate: 
                write_G4CR_to_output_file(chrom, gene_name, ancestral_seq_name, start_index, end_index, sim_output_seqs_for_this_G4CR, aligned_ancestors_for_this_G4CR, GQs, scores, real_match_rate, num_recursions)
            
            elif real_match_rate > avg_new_match_rate:
                # BELOW: Preserve all the data for the G4CRs we will be analyzing in the next recursive call
                new_sim_output_seqs.extend(sim_output_seqs_for_this_G4CR)
                new_aligned_ancestral_G4CRs.extend(aligned_ancestors_for_this_G4CR)
                
                new_gene_names.append(gene_name)
                new_ancestral_seq_names.append(ancestral_seq_name)
                new_start_and_end_indices.append((start_index, end_index))
                new_real_match_rates.append(real_match_rate)

            elif avg_new_match_rate > real_match_rate: # DE-SIMULATION
                
                # For every single mutation (insertion/deletion/substitution), reject it with a probability of x
                # Where "x" is the difference between the real match rate and the average of the match rates of the simulation outputs
                real_sim_match_rate_diff = abs(avg_new_match_rate - real_match_rate) # use this to determine probability of reverting mutations
                
                # No output, modifies the aligned ancestors and sim output sequences in-place (mutable objects)
                desimulate_output_seq(aligned_ancestors_for_this_G4CR, sim_output_seqs_for_this_G4CR, real_sim_match_rate_diff)
                write_G4CR_to_output_file(chrom, gene_name, ancestral_seq_name, start_index, end_index, sim_output_seqs_for_this_G4CR, aligned_ancestors_for_this_G4CR, GQs, scores, real_match_rate, num_recursions)
              
                
        num_recursions += 1
        
        sim_input_seqs = "3333333333".join(new_sim_output_seqs)
        aligned_ancestral_G4CRs = new_aligned_ancestral_G4CRs
        real_match_rates = new_real_match_rates
        
        gene_names = new_gene_names
        ancestral_seq_names = new_ancestral_seq_names
        start_and_end_indices = new_start_and_end_indices
        
        
main()
