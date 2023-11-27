# AUTHOR: Lucas Nelson
import os

CHROMOSOMES = ["22"]
NTS = ("A","C","G","T")

def main():
    for chrom in CHROMOSOMES:
        chrom_filepath = f"/home/mcb/users/lnelso12/evoGReg/outputs/chr{chrom}"

        # Iterate through every file in the directory
        for input_filename in os.listdir(chrom_filepath):

            output_lines = []
            input_filepath = f"/home/mcb/users/lnelso12/evoGReg/outputs/chr{chrom}/{input_filename}"
            
            with open(input_filepath, "r") as input_file:
                
                input_lines = input_file.readlines()
                human_promoter = ""

                first_line = True
                for line in input_lines:

                    if first_line:
                        first_line = False
                        output_lines.append(line)
                        continue
                    
                    line = line.strip("\n")
                    line = line.split(",")

                    if line[-1] == "":
                        line = line[:-1]
                        
                    if len(line) == 14 and line[0] == 'hg38':
                        human_promoter = line[13]
                    elif len(line) == 14:
                        ancestral_promoter = line[13]
                        num_matches = 0
                        sequence_length = 0

                        for i in range(4000):
                            if ancestral_promoter[i] == human_promoter[i]:
                                num_matches += 1
                            if ancestral_promoter[i] in NTS:
                                sequence_length += 1

                        match_rate = num_matches / sequence_length
                        line[12] = str(match_rate)
                    line = ",".join(line)
                    output_lines.append(f"{line}\n")

            with open(input_filepath, "w") as output_file:
                for line in output_lines:
                    output_file.write(line)

main()            
