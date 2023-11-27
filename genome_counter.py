# AUTHOR: Lucas Nelson
import os

def main():
    for chrom in range(3,23):
        
        chrom = str(chrom)
        chrom_filepath = f"/home/mcb/users/lnelso12/evoGReg/outputs/chr{chrom}"
        num_files_to_fix = 0
        
        # Iterate through every file in the directory
        for input_filename in os.listdir(chrom_filepath):
            input_filepath = f"/home/mcb/users/lnelso12/evoGReg/outputs/chr{chrom}/{input_filename}"
            
            with open(input_filepath, "r") as input_file:
                
                input_lines = input_file.readlines()
                genome_counter = 0
                first_line = True
                output_lines = []
                
                for line in input_lines:
                    if first_line:
                        first_line = False
                        continue
                    temp_line = line.strip("\n").split(",")

                    if len(temp_line) == 0:
                        continue
                    while len(temp_line) > 0 and temp_line[-1] == "":
                        temp_line = temp_line[:-1]
                    if len(temp_line) == 14:
                        genome_counter += 1
                    output_lines.append(temp_line)
                        
                if genome_counter != 12:
                    print(f"Gene file {input_filename} has {genome_counter} genomes")
                    num_files_to_fix += 1
                    
        print(f"Number of files in chromosome {chrom} with the wrong number of genomes: {num_files_to_fix} out of {len(os.listdir(chrom_filepath))} total gene files")
        print()
    
main()
