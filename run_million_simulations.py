# AUTHOR: Lucas Nelson
import sys
import random
import numpy as np
from simulate import simulate
from G4_greg import main as greg_main
from get_GQs import get_GQs

INPUT_G4CRs = {
                '6': ['VEGFA_1_192_214', 'VEGFA_1_1912_1941'], 
               '7': ['HOXA3_1_comp_1994_2096'], 
               '8': ['MYC_5_comp_2485_2543', 'GPIHBP1_1_323_617'], 
               '9': ['TRMT10B_1_comp_322_429'], 
               '10': ['SEMA4G_1_comp_279_349'], 
               '11': ['TMEM216_2_846_937'],
               '12': ['KRAS_2_1548_1588', 'KRAS_2_2308_2359', 'KRT7_2_comp_77_188'], 
               '13': ['PDS5B_1_2650_2730'],
               '14': ['ACIN1_1_610_693'],
               '15': ['MEGF11_1_comp_443_518'],
               '16': ['NME4_1_2434_2569'],
               '17': ['C1QL1_1_comp_654_728'],
               '18': ['BCL2_4_3382_3461', 'CTDP1_2_comp_2447_2582'], 
               '19': ['USF2_2_comp_2482_2571', 'CAPN12_1_2003_2331'], 
               '20': ['RAE1_2_3582_3708', 'ZNF512B_1_2062_2141'],
               '21': ['AIRE_1_comp_3192_3284'],
               '22': ['TTLL12_1_comp_1997_2324', 'PIK3IP1_1_comp_945_1022']
               }
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

def read_input_file(chrom_filepath, start_index, end_index, input_file):
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
  
    input_filepath = f"{chrom_filepath}/{input_file}"

    with open(input_filepath, "r") as input_file:
        
        ancestral_aligned_sequences = {}
        ancestral_match_rates = {}
        
        input_lines = input_file.readlines()
        
        for line in input_lines:
            
            if line == "":
                continue
            line = line.split(",")
            
            if len(line) == 14:
                seq_name = line[0]
                match_rate = line[12]
                aligned_sequence = line[13]
                
                ancestral_aligned_sequences[seq_name] = aligned_sequence
                ancestral_match_rates[seq_name] = match_rate

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
            raise IndexError("The modern G4CR was inserted after the human-chimp ancestor (_HP)\n")
        
        # G4CR-INDEXED DATA
        ancestral_G4CR_region = "".join([ancestral_promoter[i] for i in range(start_index, end_index) if ancestral_promoter[i] in NTS])

        real_match_rate = ancestral_match_rates[ancestral_promoter_name]
        return ancestral_G4CR_region, start_index, end_index, real_match_rate, ancestral_promoter_name

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

def write_G4CR_to_output_file(chrom, gene_name, ancestral_seq_name, start_index, end_index, sim_input_seqs_for_this_G4CR, aligned_ancestors_for_this_G4CR, GQs, scores, real_match_rate, num_recursions, local_skipped_indices, sim_call_number, num_sims):
    
    output_filepath = f"/home/mcb/users/lnelso12/G4_EvoLSTM/million_outputs/chr{chrom}/{gene_name}_{start_index}_{end_index}_run{sim_call_number}.csv"
    with open(output_filepath, "w") as output_file:
        
        # Headers outlining the REGIONAL, ISOMER, and UNIVERSAL fields
        output_file.write("SIMULATION FIELDS,num_g4crs,sim_match_rate,aligned_ancestral_seq,aligned_output_seq,max_ntot,sum_ntot,max_ntand,sum_ntands,max_tm_max,max_tm_median,max_tm_min,sum_g4cr_length,max_g4cr_length\n")
        output_file.write("G4CR FIELDS,g4cr_seq,nt_polymorphism,start_index,end_index,g4cr_length,g_percentage,ntot,ntand,tm_max,tm_median,tm_min\n")
        output_file.write("UNIVERSAL FIELDS,ancestral_seq_name,real_match_rate,num_recursions,avg_num_g4crs,avg_sim_match_rate,avg_max_ntot,avg_sum_ntots,avg_max_ntand,avg_sum_ntands,avg_max_tm_max,avg_max_tm_median,avg_max_tm_min,avg_sum_g4cr_length,avg_max_g4cr_length,avg_g_percentage\n")
        output_file.write("\n")
        
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
        
        print(f"outputting G4CR: {gene_name}_{start_index}_{end_index}")
        
        for z in range(len(sim_input_seqs_for_this_G4CR)):
            
            if z in local_skipped_indices: # Simulation with undefined output
                continue
            
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
        
        output_file.write(f"UNIVERSAL,{ancestral_seq_name},{real_match_rate},{num_recursions},{sum_num_G4CRs/num_sims},{sum_sim_match_rate/num_sims},{sum_max_ntot/num_sims},{sum_sum_ntots/num_sims},{sum_max_ntand/num_sims},{sum_sum_ntands/num_sims},{sum_max_Tm_max/num_sims},{sum_max_Tm_median/num_sims},{sum_max_Tm_min/num_sims},{sum_sum_G4CR_length/num_sims},{sum_max_G4CR_length/num_sims},{sum_avg_g_percentage/num_sims},{sum_max_g_percentage/num_sims}\n")

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
            
            # Assumes that there is never an index with both the ancestor and descendant containing a "-" gap
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
    BASIC IDEA -> start by assembling the input string of 100 * each G4CR input (from files)
    1. Run a single simulation
    2. Split output string, calculate whether each output G4CR is closer or further from real_match_rate than input
    3. Remove those that are FURTHER, run GReg on these, and output to files
    4. Build new string from those that are CLOSER, and repeat process from step 1 until string is empty
    '''
    
    chrom = sys.argv[1]
    G4CR_info = sys.argv[2].split("_") # FORMAT: 'VEGFA_1_1912_1941' or 'VEGFA_1_comp_1912_1941'
    gpu = sys.argv[3]
    sim_call_number = sys.argv[4]
    num_sims = int(sys.argv[5])
    
    
    GQs, scores = get_GQs(MAX_LOOP, MAX_BULGE, MIN_TEMP)
    
    chrom_filepath = f"/home/mcb/users/lnelso12/evoGReg/outputs/chr{chrom}"
    if len(G4CR_info) == 4:
        gene_name = "_".join([G4CR_info[0], G4CR_info[1]])
        start_index = int(G4CR_info[2])
        end_index = int(G4CR_info[3])
        
    elif len(G4CR_info) == 5:
        gene_name = "_".join([G4CR_info[0], G4CR_info[1], G4CR_info[2]])
        start_index = int(G4CR_info[3])
        end_index = int(G4CR_info[4])
        
    input_file = gene_name + ".csv"
    aligned_ancestral_G4CR, start_index, end_index, real_match_rate, ancestral_seq_name = read_input_file(chrom_filepath, start_index, end_index, input_file)
    print("Successfully read the input file and stored its information to begin simulations")
    
    real_match_rate = float(real_match_rate)
    start_index = int(start_index)
    end_index = int(end_index)
    
    sim_input_seq = [aligned_ancestral_G4CR for _ in range(num_sims)] # Will be passed in as input to simulation across every recursion
        
    skipped_indices = set()
    
    aligned_ancestral_G4CR = sim_input_seq.copy() # Will be progressively aligned to simulation outputs at each recursion 
    sim_input_seq = "3333333333".join(sim_input_seq) # length = num_G4CRs * 100
    
    num_recursions = 0
    while sim_input_seq != "":
        
        print(f"Recursion number {num_recursions}")
        sim_output_seq = simulate(gpu, sim_input_seq) # TEST SIMULATION CALL
        aligned_ancestral_G4CR, sim_output_seq = align_simulated_outputs_with_ancestors("3333333333".join(aligned_ancestral_G4CR), sim_output_seq) # Always relative to prev sim
        
        aligned_ancestral_G4CR = aligned_ancestral_G4CR.split("3333333333") # LIST
        sim_input_seq = sim_input_seq.split("3333333333")
        sim_output_seq = sim_output_seq.split("3333333333") # LIST
        retain_undefined_sims(sim_input_seq, sim_output_seq)
        
        new_match_rates = get_new_match_rates(aligned_ancestral_G4CR, sim_output_seq) # Always relative to real ancestors; LIST
            
        new_sim_output_seq = []
        new_aligned_ancestral_G4CR = []
        
        # DELETED 100-INCREMENTED FOR-LOOP HERE
        
        # Note that certain new match rates are undefined because of length 0, so skip these and their corresponding attributes in other lists
        new_match_rates_for_this_G4CR = []
        num_new_match_rates = 0
        
        for j in range(num_sims):
            new_match_rates_for_this_G4CR.append(new_match_rates[j])
            if new_match_rates[j] != -1.0:
                num_new_match_rates += 1
            else:
                skipped_indices.add(j)
        
        if len(skipped_indices) == num_sims: # Case where every simulation results in an undefined output for a given G4CR
            raise IndexError("Every simulation consists entirely of gaps, match rate is completely undefined\n")
        
        # The attributes below are different for every simulation
        avg_new_match_rate = sum([x for x in new_match_rates_for_this_G4CR if x != -1.0])/num_new_match_rates
                    
        if real_match_rate == avg_new_match_rate: 
            write_G4CR_to_output_file(chrom, gene_name, ancestral_seq_name, start_index, end_index, sim_output_seq, aligned_ancestral_G4CR, GQs, scores, real_match_rate, num_recursions, skipped_indices, sim_call_number, num_sims)
            
        elif avg_new_match_rate > real_match_rate:
            # BELOW: Preserve all the data for the G4CRs we will be analyzing in the next recursive call
            new_sim_output_seq.extend(sim_output_seq)
            new_aligned_ancestral_G4CR.extend(aligned_ancestral_G4CR)
            
        elif real_match_rate > avg_new_match_rate: # DE-SIMULATION
            
            # For every single mutation (insertion/deletion/substitution), reject it with a probability of x
            # Where "x" is the difference between the real match rate and the average of the match rates of the simulation outputs
            real_sim_match_rate_diff = abs(avg_new_match_rate - real_match_rate) # use this to determine probability of reverting mutations
            # No output, modifies the aligned ancestors and sim output sequences in-place (mutable objects)
            
            desimulate_output_seq(aligned_ancestral_G4CR, sim_output_seq, real_sim_match_rate_diff, skipped_indices)
            avg_new_sim_match_rates = calc_match_rate(aligned_ancestral_G4CR, sim_output_seq, skipped_indices)
            if avg_new_sim_match_rates == -1:
                continue
            counter = 0
            while real_match_rate > avg_new_sim_match_rates: # THIS LOOP ENSURES THAT MATCH RATE AVERAGES ARE ALWAYS SLIGHTLY HIGHER FOR SIMS THAN REAL SEQS. Ensures we are conservative.
                counter += 1
                desimulate_output_seq(aligned_ancestral_G4CR, sim_output_seq, 0.01, skipped_indices)
                avg_new_sim_match_rates = calc_match_rate(aligned_ancestral_G4CR, sim_output_seq, skipped_indices)
                if counter == 500 or avg_new_sim_match_rates == -1:
                    break
                
            write_G4CR_to_output_file(chrom, gene_name, ancestral_seq_name, start_index, end_index, sim_output_seq, aligned_ancestral_G4CR, GQs, scores, real_match_rate, num_recursions, skipped_indices, sim_call_number, num_sims)
        
        num_recursions += 1
        
        sim_input_seq = "3333333333".join(new_sim_output_seq)
        aligned_ancestral_G4CR = new_aligned_ancestral_G4CR

    print(f"The main() function ended with num_recursions value {num_recursions}\n")     

main()
