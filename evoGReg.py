# AUTHOR: Lucas Nelson
import pandas as pd
from old_greg import main as greg_main
from get_GQs import get_GQs
import sys
import os

# CHROMOSOMES = ["X"] + ["Y"] + [str(i) for i in range(1,23)]
# CHROMOSOMES.reverse()
GENOMES = ['hg38', '_HP', '_HPG', '_HPGP', '_HPGPN', '_HPGPNRMPC', '_HPGPNRMPCCS', '_HPGPNRMPCCSO', '_HPGPNRMPCCSOT', '_HPGPNRMPCCSOTSJMCMMRHCCOOO', '_HPGPNRMPCCSOTSJMCMMRHCCOOOSVCTOPBOCECFCMAOLPPEMMESC', '_HPGPNRMPCCSOTSJMCMMRHCCOOOSVCTOPBOCECFCMAOLPPEMMESCLETCEOD']
MAX_LOOP = 7
MAX_BULGE = 3
MIN_TEMP = 50

def evoGReg():

    input_filepath = sys.argv[1] # Filepath to TSS file
    chrom_value = sys.argv[2]

    with open(input_filepath, "r") as tss_file:
        GQs, scores = get_GQs(MAX_LOOP, MAX_BULGE, MIN_TEMP)
        aligned_human_seq = ""
        chromosome = f"chr{chrom_value}"
        alignment = pd.read_pickle(f'/home/mcb/users/lnelso12/conservation/data/seqDictPad_{chromosome}.pkl')

        gene_infos = [gene_info.strip("\n").split("\t") for gene_info in tss_file.readlines()]

        # Minimize number of pandas initializations by loading in all G4CRs of a given chrosome genome-by-genome
        for genome in GENOMES:

            chrom_in_current_genome = alignment[genome]

            prev_gene_name = "placeholder_1"
            for gene_info in gene_infos:

                # gene_info = gene_info.split("\t")
                if len(gene_info) == 6:
                    gene_info.pop(4) # Remove the excess "1" present in certain rows

                chr_from_file = gene_info[0]
                tss = int(gene_info[1])

                gene_name = gene_info[3]
                gene_name_split = gene_name.split("_")
                prev_gene_name_split = prev_gene_name.split("_")

                if chr_from_file != chromosome:
                    break # Stop checking G4CRs if at end of current chromosome
                elif gene_name_split[0] == prev_gene_name_split[0]:
                    continue # Skip this line if the gene in question has already been analyzed

                # Make sure the ./outputs/chr? filepath exists prior to starting code, and that the file itself is empty/non-existent
                main_filepath = f"./outputs/{chromosome}/{gene_name}.csv"
                comp_filepath = f"./outputs/{chromosome}/{gene_name}_comp.csv"

                with open(main_filepath, "a") as main_file, open(comp_filepath, "a") as comp_file:
                    if os.path.getsize(main_filepath) == 0:
                        main_file.write("Name of Sequence,G4CR Sequence,NT Polymorphism Values,G4CR Start Index,G4CR End Index,G4CR Length in NTs,G-content Percentage,N-total,N-tandem,Tm-max,Tm-median,Tm-min,Match Rate,Aligned Promoter Sequence\n")
                    if os.path.getsize(comp_filepath) == 0:
                        comp_file.write("Name of Sequence,G4CR Sequence,NT Polymorphism Values,G4CR Start Index,G4CR End Index,G4CR Length in NTs,G-content Percentage,N-total,N-tandem,Tm-max,Tm-median,Tm-min,Match Rate,Aligned Promoter Sequence\n")

                    # aligned promoter regions extracted from sequence alignment
                    aligned_seq = [chrom_in_current_genome[i] for i in range((tss-2000),(tss+2000))]
                    condensed_seq = [] # promoter regions with gaps removed
                    condensed_indices = [] # array that, for each nt in the "condensed" sequence, records its index value in the original aligned sequence
                    nts = ("A","C","G","T")

                    for i in range(len(aligned_seq)):
                        if aligned_seq[i] in nts:
                            condensed_seq.append(aligned_seq[i])
                            condensed_indices.append(i)

                    # IMPORTANT: sequences that have too few nts (i.e., close to 4000 gaps) break GReg.
                    # We want to ignore such sequences and continue analyzing regions that are more relevant
                    if len(condensed_seq) < 15:
                        continue
                    aligned_seq_string = "".join(aligned_seq)
                    if genome == "hg38":
                        aligned_human_seq = aligned_seq_string # Keep a record of human genome to later calculate match rate
                        match_rate = 1.0
                    else: # Calculate match rate between the aligned ancestral sequence and the human sequence
                        num_matches = 0
                        sequence_length = 0
                        for i in range(4000):
                            if aligned_human_seq[i] == aligned_seq_string[i]:
                                num_matches += 1
                            if aligned_seq_string[i] in nts:
                                sequence_length += 1
                        match_rate = num_matches / sequence_length

                    # Generate complementary strand of the condensed sequence
                    comp_aligned_seq = []
                    for i in range(len(aligned_seq)):
                        if aligned_seq[i] == "A":
                            comp_aligned_seq.append("T")
                        elif aligned_seq[i] == "C":
                            comp_aligned_seq.append("G")
                        elif aligned_seq[i] == "G":
                            comp_aligned_seq.append("C")
                        elif aligned_seq[i] == "T":
                            comp_aligned_seq.append("A")
                        else:
                            comp_aligned_seq.append(aligned_seq[i])

                    comp_aligned_seq_string = "".join(comp_aligned_seq)
                    comp_condensed_seq = [nt for nt in comp_aligned_seq if nt in nts]
                    condensed_seq_string = "".join(condensed_seq)
                    comp_condensed_seq_string = "".join(comp_condensed_seq)

                    main_output = greg_main(condensed_seq_string, genome, MAX_LOOP, GQs, scores)
                    comp_output = greg_main(comp_condensed_seq_string, genome, MAX_LOOP, GQs, scores)

                    # Replace the start/end indices of G4CRs with the human-aligned values from condensed_indices
                    main_lines = []
                    for index,line in enumerate(main_output):

                        if index == 0:
                            line.append(str(match_rate))
                            line.append(aligned_seq_string)

                        start_index = line[3]
                        end_index = line[4]

                        line[3] = condensed_indices[int(start_index)-1]
                        line[4] = condensed_indices[int(end_index)-1]
                        line = [str(x) for x in line]

                        line.append("\n")
                        main_line = ",".join(line)
                        main_lines.append(main_line)

                    comp_lines = []
                    for index,line in enumerate(comp_output):

                        if index == 0:
                            line.append(str(match_rate))
                            line.append(comp_aligned_seq_string)

                        start_index = line[3]
                        end_index = line[4]

                        line[3] = condensed_indices[int(start_index)-1]
                        line[4] = condensed_indices[int(end_index)-1]
                        line = [str(x) for x in line]

                        line.append("\n")
                        comp_line = ",".join(line)
                        comp_lines.append(comp_line)

                    main_output = "".join(main_lines)
                    comp_output = "".join(comp_lines)

                    if main_output != "":
                        main_file.write(main_output)
                    else:
                        main_file.write(f"{genome},,,,,,,,,,,,{match_rate},{aligned_seq_string}\n")

                    if comp_output != "":
                        comp_file.write(comp_output)
                    else:
                        comp_file.write(f"{genome},,,,,,,,,,,,{match_rate},{comp_aligned_seq_string}\n")

                    # Write important info at the bottom of the file. Leave whitespace so that we can recognize it when reading for simulations
                    if genome == '_HPGPNRMPCCSOTSJMCMMRHCCOOOSVCTOPBOCECFCMAOLPPEMMESCLETCEOD':
                        main_file.write("\nTSS index,Chromosome Number,Gene Name/Number\n")
                        main_file.write(f"{tss},{chrom_value},{gene_name}\n")
                        comp_file.write("\nTSS index,Chromosome Number,Gene Name/Number\n")
                        comp_file.write(f"{tss},{chrom_value},{gene_name}\n")

                prev_gene_name = gene_name




evoGReg()
