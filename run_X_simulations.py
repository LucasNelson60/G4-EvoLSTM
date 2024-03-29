# AUTHOR: Lucas Nelson
import sys
import os
import random
import numpy as np
from simulate import simulate
from greg import main as greg_main
from get_GQs import get_GQs
import copy

ANCESTRAL_GENOMES = ['hg38', '_HP', '_HPG', '_HPGP', '_HPGPN', '_HPGPNRMPC', '_HPGPNRMPCCS', '_HPGPNRMPCCSO', '_HPGPNRMPCCSOT', '_HPGPNRMPCCSOTSJMCMMRHCCOOO', '_HPGPNRMPCCSOTSJMCMMRHCCOOOSVCTOPBOCECFCMAOLPPEMMESC', ] # '_HPGPNRMPCCSOTSJMCMMRHCCOOOSVCTOPBOCECFCMAOLPPEMMESCLETCEOD'
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

def read_input_files(chrom_filepath):
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
    gene_names = [] 
    start_and_end_indices = []
    real_match_rates = []
    oldest_ancestral_G4CRs = []
    oldest_real_match_rates = []
    poly_vals = []
    G4CR_lengths = []
    g_percentages = []
    ntots = []
    ntands = []
    tm_maxes = []
    tm_medians = []
    tm_mins = []
    ancestral_G4CR_regions = []

    for input_file in os.listdir(chrom_filepath):
        
        input_filepath = f"{chrom_filepath}/{input_file}"
        file_name = input_filepath.split("/")[-1]
        clipped_file_name = file_name.split(".")[0].split("_")
        
        if len(clipped_file_name) == 4:
            gene_name = "_".join(clipped_file_name[:2])
            start_index = int(clipped_file_name[2])
            end_index = int(clipped_file_name[3])
        elif len(clipped_file_name) == 5:
            gene_name = "_".join(clipped_file_name[:3])
            start_index = int(clipped_file_name[3])
            end_index = int(clipped_file_name[4])
        else:
            continue
        

        with open(input_filepath, "r") as input_file:
            
            match_rates_for_this_G4CR = {}
            poly_vals_for_this_G4CR = {}
            this_G4CR_lengths = {}
            g_percentages_for_this_G4CR = {}
            ntots_for_this_G4CR = {}
            ntands_for_this_G4CR = {}
            tm_maxes_for_this_G4CR = {}
            tm_medians_for_this_G4CR = {}
            tm_mins_for_this_G4CR = {}
            ancestral_regions_for_this_G4CR = {}
            
            # Input file format: ancestral_promoter_name, real_match_rate, G4CR_seq
            lines = input_file.readlines()
            
            for i in range(len(lines)):
                line = lines[i].strip("\n").split(",")
                seq_name = line[0]

                if i == 0:
                    oldest_ancestor = line[10]
                    oldest_real_match_rate = float(line[11])
                    continue
                
                match_rates_for_this_G4CR[seq_name] = line[1]
                poly_vals_for_this_G4CR[seq_name] = line[2]
                this_G4CR_lengths[seq_name] = line[3]
                g_percentages_for_this_G4CR[seq_name] = line[4]
                ntots_for_this_G4CR[seq_name] = line[5]
                ntands_for_this_G4CR[seq_name] = line[6]
                tm_maxes_for_this_G4CR[seq_name] = line[7]
                tm_medians_for_this_G4CR[seq_name] = line[8]
                tm_mins_for_this_G4CR[seq_name] = line[9]
                ancestral_regions_for_this_G4CR[seq_name] = line[10]
                
            
            oldest_ancestral_G4CRs.append(oldest_ancestor)
            oldest_real_match_rates.append(oldest_real_match_rate)
            gene_names.append(gene_name)
            start_and_end_indices.append((start_index,end_index))
            real_match_rates.append(match_rates_for_this_G4CR)
            poly_vals.append(poly_vals_for_this_G4CR)
            G4CR_lengths.append(this_G4CR_lengths)
            g_percentages.append(g_percentages_for_this_G4CR)
            ntots.append(ntots_for_this_G4CR)
            ntands.append(ntands_for_this_G4CR)
            tm_maxes.append(tm_maxes_for_this_G4CR)
            tm_medians.append(tm_medians_for_this_G4CR)
            tm_mins.append(tm_mins_for_this_G4CR)
            ancestral_G4CR_regions.append(ancestral_regions_for_this_G4CR) 
            
            
    return gene_names, start_and_end_indices, real_match_rates, oldest_ancestral_G4CRs, oldest_real_match_rates, poly_vals, G4CR_lengths, g_percentages, ntots, ntands, tm_maxes, tm_medians, tm_mins, ancestral_G4CR_regions 

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
    
def mark_skipped_indices(sim_input_seqs, sim_output_seqs):
    '''
    get_new_match_rates(sim_input_seqs: List[str], sim_output_seqs: List[str]) -> List[float]
    IMPORTANT: len(sim_input_seqs) = len(sim_output_seqs)
    
    Given the simulation inputs and outputs, in list format (split by "3333333333"),
    Returns the match rate of each output sim compared to its corresponding input
    '''
    new_match_rates = []
    
    for i in range(len(sim_output_seqs)):

        sim_input_seq = sim_input_seqs[i]
        sim_output_seq = sim_output_seqs[i]
        
        # Get the match rate between the input and output seqs
        ungapped_length = 0
        k = 0
        matches = 0
        skip = False
        while k < len(sim_input_seq):
            try:
                if sim_output_seq[k] not in NTS:
                    k += 1
                    continue
                elif sim_output_seq[k] == sim_input_seq[k]:
                    matches += 1
                ungapped_length += 1
                k += 1
            except IndexError:
                skip = True
                break
        if skip:
            continue
            
        # Skip simulations that have zero ungapped_length, i.e., that consist of all "-" characters, these are useless
        if ungapped_length == 0:
            new_match_rates.append(-1.0) # Placeholder value, will skip when checking sim outputs in main()
            continue
        
        new_match_rate = matches/ungapped_length
        new_match_rates.append(new_match_rate)
    
    return new_match_rates

def write_G4CR_to_output_file(output_filepath, gene_name, poly_vals, this_G4CR_length, g_percentage, this_ntot, this_ntand, this_tm_max, this_tm_median, this_tm_min, G4CR_seq, \
    ancestral_seq_name, start_index, end_index, sim_output_seqs_for_this_G4CR, aligned_ancestors_for_this_G4CR, GQs, scores, real_match_rate, num_recursions, local_skipped_indices, oldest_match_rate_rel_to_human):

    chrom_filepath = f"{output_filepath}/{ancestral_seq_name}/{gene_name}_{start_index}_{end_index}.csv"
    with open(chrom_filepath, "w") as output_file:
        
        # Headers outlining the REGIONAL, ISOMER, and UNIVERSAL fields
        output_file.write("SIMULATION FIELDS,num_g4crs,sim_match_rate,aligned_ancestral_seq,aligned_output_seq,max_ntot,sum_ntot,max_ntand,sum_ntands,max_tm_max,max_tm_median,max_tm_min,sum_g4cr_length,max_g4cr_length\n")
        output_file.write("G4CR FIELDS,g4cr_seq,nt_polymorphism,start_index,end_index,g4cr_length,g_percentage,ntot,ntand,tm_max,tm_median,tm_min\n")
        
        sum_num_G4CRs = 0
        sum_sim_match_rate = 0.0
        sum_max_ntot = 0
        sum_sum_ntots = 0
        sum_max_ntand = 0
        sum_sum_ntands = 0
        sum_max_Tm_max = 0.0
        sum_max_Tm_median = 0.0
        sum_max_Tm_min = 0.0
        sum_sum_G4CR_length = 0
        sum_max_G4CR_length = 0
        sum_avg_g_percentage = 0.0
        sum_max_g_percentage = 0
        
        # print(f"outputting G4CR: {gene_name}_{start_index}_{end_index} for ancestor {ancestral_seq_name}")

        seq_counter = 0
        for z in range(len(sim_output_seqs_for_this_G4CR)):
            
            if z in local_skipped_indices: # Simulation with undefined output
                continue
            
            sim_input_seq = sim_output_seqs_for_this_G4CR[z] # Use input, not output, because its branch length is closer to the real one
            aligned_ancestor = aligned_ancestors_for_this_G4CR[z]
            
            # Split the lists to be able to replace mutations as needed
            aligned_ancestor = aligned_ancestors_for_this_G4CR[z]
            sim_output = sim_output_seqs_for_this_G4CR[z]
            
            matches = 0
            ungapped_length = 0
            for j in range(len(aligned_ancestor)):
                if aligned_ancestor[j] != "-" and aligned_ancestor[j] == sim_output[j]:
                    matches += 1
                if aligned_ancestor[j] != "-" and sim_output[j] != "-":
                    ungapped_length += 1
            
            # Skip simulations that have zero ungapped_length, i.e., that consist of all "-" characters, to avoid ZeroDivisionError below
            if ungapped_length == 0:
                continue

            sim_match_rate = matches/ungapped_length
            seq_counter += 1
            
            # Run GReg on the simulation output
            gapless_output = "".join([nt for nt in sim_input_seq if nt != "-"])
            G4CR_fields, num_G4CRs, max_ntot, sum_ntots, max_ntand, sum_ntands, max_Tm_max, max_Tm_median, max_Tm_min, sum_G4CR_length, max_G4CR_length, avg_g_percentage, max_g_percentage = greg_main(gapless_output, MAX_LOOP, GQs, scores)
            
            output_file.write(f"SIMULATION,{num_G4CRs},{sim_match_rate},{''.join(aligned_ancestor)},{''.join(sim_input_seq)},{max_ntot},{sum_ntots},{max_ntand},{sum_ntands},{max_Tm_max},{max_Tm_median},{max_Tm_min},{sum_G4CR_length},{max_G4CR_length}\n")
            output_file.write(G4CR_fields)
            
            sum_num_G4CRs += int(num_G4CRs)
            sum_sim_match_rate += float(sim_match_rate)
            sum_max_ntot += int(max_ntot)
            sum_sum_ntots += int(sum_ntots)
            sum_max_ntand += int(max_ntand)
            sum_sum_ntands += int(sum_ntands)
            sum_max_Tm_max += float(max_Tm_max)
            sum_max_Tm_median += float(max_Tm_median)
            sum_max_Tm_min += float(max_Tm_min)
            sum_sum_G4CR_length += int(sum_G4CR_length)
            sum_max_G4CR_length += int(max_G4CR_length)
            sum_avg_g_percentage += float(avg_g_percentage)
            sum_max_g_percentage += int(max_g_percentage)
            
        try:
            avg_num_G4CRs = sum_num_G4CRs/seq_counter
            avg_sim_match_rate = sum_sim_match_rate/seq_counter
            avg_max_ntot = sum_max_ntot/seq_counter
            avg_sum_ntots = sum_sum_ntots/seq_counter
            avg_max_ntand = sum_max_ntand/seq_counter
            avg_sum_ntands = sum_sum_ntands/seq_counter
            avg_max_Tm_max = sum_max_Tm_max/seq_counter
            avg_max_Tm_median = sum_max_Tm_median/seq_counter
            avg_max_Tm_min = sum_max_Tm_min/seq_counter
            avg_sum_G4CR_length = sum_sum_G4CR_length/seq_counter
            avg_max_G4CR_length = sum_max_G4CR_length/seq_counter
            avg_avg_g_percentage = sum_avg_g_percentage/seq_counter
            avg_max_g_percentage = sum_max_g_percentage/seq_counter
        except ZeroDivisionError:
            return
        
        output_file.write("\n")
        output_file.write("UNIVERSAL FIELD,ancestral_seq_name,num_recursions,avg_num_g4crs,avg_sim_match_rate,avg_max_ntot,avg_sum_ntots,avg_max_ntand,avg_sum_ntands,avg_max_tm_max,avg_max_tm_median,avg_max_tm_min,avg_sum_g4cr_length,avg_max_g4cr_length,avg_g_percentage,avg_max_g_percentage\n")
        output_file.write("REAL FIELD,ancestral_seq_name,real_match_rate,oldest_match_rate_relative_to_human,poly_vals,G4CR_length,g_percentage,ntot,ntand,tm_max,tm_median,tm_min,G4CR_seq\n")
        output_file.write(f"UNIVERSAL,{ancestral_seq_name},{num_recursions},{avg_num_G4CRs},{avg_sim_match_rate},{avg_max_ntot},{avg_sum_ntots},{avg_max_ntand},{avg_sum_ntands},{avg_max_Tm_max},{avg_max_Tm_median},\
        {avg_max_Tm_min},{avg_sum_G4CR_length},{avg_max_G4CR_length},{avg_avg_g_percentage},{avg_max_g_percentage}\n")
        output_file.write(f"REAL,{ancestral_seq_name},{real_match_rate},{oldest_match_rate_rel_to_human},{poly_vals},{this_G4CR_length},{g_percentage},{this_ntot},{this_ntand},{this_tm_max},{this_tm_median},{this_tm_min},{G4CR_seq}\n")
        
def calc_match_rate(aligned_ancestors_for_this_G4CR, sim_output_seqs_for_this_G4CR, local_skipped_indices):

    new_sim_match_rates = []
    for k in range(len(sim_output_seqs_for_this_G4CR)):
        if k in local_skipped_indices: # Simulation with undefined output
            continue
        
        # Split the lists to be able to replace mutations as needed
        aligned_ancestor = aligned_ancestors_for_this_G4CR[k]
        sim_output = sim_output_seqs_for_this_G4CR[k]
        
        matches = 0
        ungapped_count = 0
        for j in range(len(aligned_ancestor)):
            if aligned_ancestor[j] != "-" and aligned_ancestor[j] == sim_output[j]:
                matches += 1
            if aligned_ancestor[j] != "-" and sim_output[j] != "-":
                ungapped_count += 1
        
        if ungapped_count == 0:
            local_skipped_indices.add(k)
            continue
        
        new_sim_match_rate = matches / ungapped_count
        new_sim_match_rates.append(new_sim_match_rate)
        
    if len(new_sim_match_rates) == 0:
        return -1
    return sum(new_sim_match_rates) / len(new_sim_match_rates)
    
def desimulate_output_seq(aligned_ancestors_for_this_G4CR, sim_output_seqs_for_this_G4CR, desimulating_distance, local_skipped_indices):
    
    for k in range(len(sim_output_seqs_for_this_G4CR)): # iterate through each of the simulations
        
        if k in local_skipped_indices: # Simulation with undefined output
            continue
        
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
                
                if random.random() <= desimulating_distance:
                        reversing_indel = True

                # Assumes that there is never an index with both the ancestor and descendant containing a "-" gap
                if reversing_indel and ancestor_char == "-":
                    continue
                elif reversing_indel and sim_char == "-":
                    new_aligned_ancestor.append(ancestor_char)
                    new_sim_output.append(ancestor_char)
                else:
                    new_aligned_ancestor.append(ancestor_char)
                    new_sim_output.append(sim_char)
                    
            elif ancestor_char != sim_char: # substitution
                if random.random() <= desimulating_distance:
                    new_sim_output.append(ancestor_char)
                else:
                    new_sim_output.append(sim_char)
                new_aligned_ancestor.append(ancestor_char)
                
            elif ancestor_char == sim_char:
                new_aligned_ancestor.append(ancestor_char)
                new_sim_output.append(sim_char)
                
                if reversing_indel:
                    reversing_indel = False
        
        # Replace simulations for this G4CR after de-simulating them by desired branch length
        aligned_ancestors_for_this_G4CR[k] = new_aligned_ancestor
        sim_output_seqs_for_this_G4CR[k] = new_sim_output
  
def retain_undefined_sims(sim_input_seqs, sim_output_seqs):
    ''' Ensures that no inputs that were previously undefined (i.e., consisting of all gaps) has nts inserted in latest sim output
    '''
    for i in range(len(sim_input_seqs)):
        input_seq = sim_input_seqs[i]
        output_seq = sim_output_seqs[i]
        if input_seq.count("-") == len(input_seq):
            sim_output_seqs[i] = ''.join(['-' for _ in range(len(output_seq))])

def main():
    '''
    BASIC IDEA -> start by assembling the input string of num_sims * each G4CR input (from files)
    1. Run a single simulation
    2. Split output string, calculate whether each output G4CR is closer or further from real_match_rate than input
    3. Remove those that are FURTHER, run GReg on these, and output to files
    4. Build new string from those that are CLOSER, and repeat process from step 1 until string is empty
    '''
    
    # chrom = sys.argv[1]
    input_filepath = sys.argv[1]
    output_filepath = sys.argv[2]
    gpu = sys.argv[3]
    num_sims = int(sys.argv[4])
    
    GQs, scores = get_GQs(MAX_LOOP, MAX_BULGE, MIN_TEMP)
    
    # chrom_filepath = f"/home/mcb/users/lnelso12/evoGReg/outputs_filtered/chr{chrom}"
    
    # NOTE THAT EVERY OUTPUT OF THIS FUNCTION IS G4CR-INDEXED AND OF THE SAME LENGTH
    gene_names, start_and_end_indices, real_match_rates, oldest_ancestral_G4CRs, oldest_match_rates_rel_to_human, poly_vals, G4CR_lengths, g_percentages, ntots, ntands, tm_maxes, tm_medians, tm_mins, ancestral_G4CR_regions = read_input_files(input_filepath)
    # NOTE on the line above: the "oldest_ancestral_G4CRs" only consists of the furthest-back ancestral aligned sequence, which we use as input to simulate on.
    # Meanwhile, the "ancestral_G4CR_regions" is a dictionary containing the aligned sequences corresponding to EVERY ancestor for a given G4CR, to be written into the output file
    print("Successfully read all input files and stored their information to begin simulations")
    
    sim_input_seqs = [] # Will be passed in as input to simulation across every recursion
    for i in range(len(oldest_ancestral_G4CRs)): # num_sims simulation inputs for each G4CR
        sim_input_seqs.extend([oldest_ancestral_G4CRs[i] for _ in range(num_sims)])
        
    global_skipped_indices = [set() for _ in range(len(oldest_ancestral_G4CRs))]
    
    oldest_ancestral_G4CRs = sim_input_seqs.copy() # Will be progressively aligned to simulation outputs at each recursion 
    sim_input_seqs = "3333333333".join(sim_input_seqs) # length = num_G4CRs * num_sims
    
    num_recursions = 0
    while sim_input_seqs != "":
        
        print(f"Recursion number {num_recursions}")
        sim_output_seqs = simulate(gpu, sim_input_seqs) # TEST SIMULATION CALL
        oldest_ancestral_G4CRs, sim_output_seqs = align_simulated_outputs_with_ancestors("3333333333".join(oldest_ancestral_G4CRs), sim_output_seqs) # Always relative to prev sim
        
        oldest_ancestral_G4CRs = oldest_ancestral_G4CRs.split("3333333333") # LIST
        sim_input_seqs = sim_input_seqs.split("3333333333")
        sim_output_seqs = sim_output_seqs.split("3333333333") # LIST
        retain_undefined_sims(sim_input_seqs, sim_output_seqs)
        
        new_match_rates = mark_skipped_indices(oldest_ancestral_G4CRs, sim_output_seqs) # Always relative to real ancestors; LIST
            
        new_sim_output_seqs = []
        new_oldest_ancestral_G4CRs = []
        
        new_gene_names = []
        new_start_and_end_indices = []
        new_real_match_rates = []
        new_oldest_match_rates_rel_to_human = []
        new_poly_vals = []
        new_G4CR_lengths = []
        new_g_percentages = []
        new_ntots = []
        new_ntands = []
        new_tm_maxes = []
        new_tm_medians = []
        new_tm_mins = []
        new_ancestral_G4CR_regions = []
        
        for i in range(num_sims, len(new_match_rates)+num_sims, num_sims):
            
            # Note that certain new match rates are undefined because of length 0, so skip these and their corresponding attributes in other lists
            new_match_rates_for_this_G4CR = []
            num_new_match_rates = 0
            local_skipped_indices = global_skipped_indices[int(i/num_sims)-1]
            
            for j in range(i-num_sims,i):
                new_match_rates_for_this_G4CR.append(new_match_rates[j])
                if new_match_rates[j] != -1.0:
                    num_new_match_rates += 1
                else:
                    local_skipped_indices.add(j)
            
            if len(local_skipped_indices) == num_sims: # Case where every simulation results in an undefined output for a given G4CR
                continue # These will be set aside and skipped to be studied individually in closer detail
            
            # The attributes below are different for every simulation
            sim_output_seqs_for_this_G4CR = [sim_output_seqs[l] for l in range(i-num_sims,i)]
            aligned_ancestors_for_this_G4CR = [oldest_ancestral_G4CRs[b] for b in range(i-num_sims,i)]
            
            # The attributes below are the same for all simulations of a given G4CR
            real_match_rates_for_this_G4CR = real_match_rates[int(i/num_sims)-1]
            oldest_match_rate_rel_to_human = oldest_match_rates_rel_to_human[int(i/num_sims)-1]
            poly_vals_for_this_G4CR = poly_vals[int(i/num_sims)-1]
            this_G4CR_lengths = G4CR_lengths[int(i/num_sims)-1]
            g_percentages_for_this_G4CR = g_percentages[int(i/num_sims)-1]
            ntots_for_this_G4CR = ntots[int(i/num_sims)-1]
            ntands_for_this_G4CR = ntands[int(i/num_sims)-1]
            tm_maxes_for_this_G4CR = tm_maxes[int(i/num_sims)-1]
            tm_medians_for_this_G4CR = tm_medians[int(i/num_sims)-1]
            tm_mins_for_this_G4CR = tm_mins[int(i/num_sims)-1]
            ancestral_regions_for_this_G4CR = ancestral_G4CR_regions[int(i/num_sims)-1]
            gene_name = gene_names[int(i/num_sims)-1]
            start_index = int(start_and_end_indices[int(i/num_sims)-1][0])
            end_index = int(start_and_end_indices[int(i/num_sims)-1][1])
            
            avg_new_match_rate = calc_match_rate(aligned_ancestors_for_this_G4CR, sim_output_seqs_for_this_G4CR, local_skipped_indices)
            
            keys_to_remove = []
            for ancestral_seq_name in real_match_rates_for_this_G4CR.keys():
                
                real_match_rate = float(real_match_rates_for_this_G4CR[ancestral_seq_name])
                poly_val = poly_vals_for_this_G4CR[ancestral_seq_name]
                G4CR_length = this_G4CR_lengths[ancestral_seq_name]
                g_percentage = g_percentages_for_this_G4CR[ancestral_seq_name]
                ntot = ntots_for_this_G4CR[ancestral_seq_name]
                ntand = ntands_for_this_G4CR[ancestral_seq_name]
                tm_max = tm_maxes_for_this_G4CR[ancestral_seq_name]
                tm_median = tm_medians_for_this_G4CR[ancestral_seq_name]
                tm_min = tm_mins_for_this_G4CR[ancestral_seq_name]
                G4CR_seq = ancestral_regions_for_this_G4CR[ancestral_seq_name]
                    
                if real_match_rate == avg_new_match_rate: 
                    write_G4CR_to_output_file(output_filepath, gene_name, poly_val, G4CR_length, g_percentage, ntot, ntand, tm_max, tm_median, tm_min, G4CR_seq, \
                        ancestral_seq_name, start_index, end_index, sim_output_seqs_for_this_G4CR, aligned_ancestors_for_this_G4CR, GQs, scores, real_match_rate, num_recursions, local_skipped_indices, oldest_match_rate_rel_to_human)
                    keys_to_remove.append(ancestral_seq_name)
                    
                elif real_match_rate > avg_new_match_rate: # DE-SIMULATION
                    # print(f"Real match rate of {real_match_rate} calculated for ancestor {ancestral_seq_name}")
                    # print(f"Match rate of {avg_new_match_rate} calculated for ancestor {ancestral_seq_name} on iteration number 0")
                    
                    desimulated_aligned_ancestors_for_this_G4CR = copy.deepcopy(aligned_ancestors_for_this_G4CR)
                    desimulated_sim_output_seqs_for_this_G4CR = copy.deepcopy(sim_output_seqs_for_this_G4CR)
                    
                    avg_new_sim_match_rates = avg_new_match_rate
                    
                    desimulating_distance = 0.01
                    if avg_new_sim_match_rates == -1:
                        keys_to_remove.append(ancestral_seq_name)
                        continue
                    
                    counter = 0
                    while real_match_rate > avg_new_sim_match_rates: # THIS LOOP ENSURES THAT MATCH RATE AVERAGES ARE ALWAYS SLIGHTLY HIGHER FOR SIMS THAN REAL SEQS. Ensures we are conservative.
                        counter += 1
                        desimulate_output_seq(desimulated_aligned_ancestors_for_this_G4CR, desimulated_sim_output_seqs_for_this_G4CR, desimulating_distance, local_skipped_indices)
                        avg_new_sim_match_rates = calc_match_rate(desimulated_aligned_ancestors_for_this_G4CR, desimulated_sim_output_seqs_for_this_G4CR, local_skipped_indices)
                        # print(f"Match rate of {avg_new_sim_match_rates} calculated for ancestor {ancestral_seq_name} on iteration number {counter}")
                        if counter == 500 or avg_new_sim_match_rates == -1:
                            break
                        desimulating_distance += 0.01
                    write_G4CR_to_output_file(output_filepath, gene_name, poly_val, G4CR_length, g_percentage, ntot, ntand, tm_max, tm_median, tm_min, G4CR_seq, \
                        ancestral_seq_name, start_index, end_index, desimulated_sim_output_seqs_for_this_G4CR, desimulated_aligned_ancestors_for_this_G4CR, GQs, scores, real_match_rate, num_recursions, local_skipped_indices, oldest_match_rate_rel_to_human)
                    keys_to_remove.append(ancestral_seq_name)
            
            for ancestral_seq_name in keys_to_remove:
                real_match_rates_for_this_G4CR.pop(ancestral_seq_name)
            
            if bool(real_match_rates_for_this_G4CR):
                # BELOW: Preserve all the data for the G4CRs we will be analyzing in the next recursive call
                new_sim_output_seqs.extend(sim_output_seqs_for_this_G4CR)
                new_oldest_ancestral_G4CRs.extend(aligned_ancestors_for_this_G4CR)
                
                new_gene_names.append(gene_name)
                new_start_and_end_indices.append((start_index, end_index))
                new_real_match_rates.append(real_match_rates_for_this_G4CR)
                new_oldest_match_rates_rel_to_human.append(oldest_match_rate_rel_to_human)
                new_poly_vals.append(poly_vals_for_this_G4CR)
                new_G4CR_lengths.append(this_G4CR_lengths)
                new_g_percentages.append(g_percentages_for_this_G4CR)
                new_ntots.append(ntots_for_this_G4CR)
                new_ntands.append(ntands_for_this_G4CR)
                new_tm_maxes.append(tm_maxes_for_this_G4CR)
                new_tm_medians.append(tm_medians_for_this_G4CR)
                new_tm_mins.append(tm_mins_for_this_G4CR)
                new_ancestral_G4CR_regions.append(ancestral_regions_for_this_G4CR)
                
            
        num_recursions += 1
        
        sim_input_seqs = "3333333333".join(new_sim_output_seqs)
        oldest_ancestral_G4CRs = new_oldest_ancestral_G4CRs
        real_match_rates = new_real_match_rates
        print(real_match_rates)
        
        gene_names = new_gene_names
        start_and_end_indices = new_start_and_end_indices
        oldest_match_rates_rel_to_human = new_oldest_match_rates_rel_to_human
        poly_vals = new_poly_vals
        G4CR_lengths = new_G4CR_lengths
        g_percentages = new_g_percentages
        ntots = new_ntots
        ntands = new_ntands
        tm_maxes = new_tm_maxes
        tm_medians = new_tm_medians
        tm_mins = new_tm_mins
        ancestral_G4CR_regions = new_ancestral_G4CR_regions

    print(f"The main() function ended with num_recursions value {num_recursions}\n")    
    
     
main()

# COMMAND TO RUN SIMS ON chr_
'''
python3 run_X_simulations.py /home/mcb/users/lnelso12/evoGReg/outputs_ancestral_formatted/chr11 /home/mcb/users/lnelso12/G4_EvoLSTM/new_outputs/hundred_outputs/chr11 5 100
'''