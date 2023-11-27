# AUTHOR: Lucas Nelson
import os

def main():
    
    evolstm_outputs_directory = '/home/mcb/users/lnelso12/G4_EvoLSTM/outputs'
    evogreg_outputs_directory = '/home/mcb/users/lnelso12/evoGReg/outputs'
    evogreg_filtered_outputs_directory = '/home/mcb/users/lnelso12/evoGReg/outputs_filtered'

    for chr_ in os.listdir(evolstm_outputs_directory):

        chrom_filepath = f'{evolstm_outputs_directory}/{chr_}'
        
        for input_file in os.listdir(chrom_filepath):
            
            higher_ntot_simulated = False
            file_name_elements = input_file.split(".")[0].split("_")
            
            if len(file_name_elements) == 4:
                gene_name = "_".join([file_name_elements[0], file_name_elements[1]])
                start_index = int(file_name_elements[2])
                end_index = int(file_name_elements[3])
            elif len(file_name_elements) == 5:
                gene_name = "_".join([file_name_elements[0], file_name_elements[1], file_name_elements[2]])
                start_index = int(file_name_elements[3])
                end_index = int(file_name_elements[4])
            else:
                continue
            
            sim_output_filepath = f'{chrom_filepath}/{input_file}'
            evogreg_output_filepath = f'{evogreg_outputs_directory}/{chr_}/{gene_name}.csv'
            evogreg_filtered_output_filepath = f'{evogreg_filtered_outputs_directory}/{chr_}/{gene_name}_{start_index}_{end_index}.csv'
            
            with open(sim_output_filepath, "r") as sim_outputs, open(evogreg_output_filepath, "r") as evogreg_output:
                
                sim_output_lines = sim_outputs.readlines()
                evogreg_output_lines = evogreg_output.readlines()
                
                real_ntot = -1.0
                for evogreg_output_line in evogreg_output_lines:
                    evogreg_line_list = evogreg_output_line.split(",")
                    if evogreg_line_list[0] == "hg38" and start_index <= int(evogreg_line_list[3]) and end_index >= int(evogreg_line_list[4]):
                        real_ntot = float(evogreg_line_list[7])
                        break
                
                for sim_output_line in sim_output_lines:
                    sim_line_list = sim_output_line.split(",")
                    if sim_line_list[0] == "G4CR" and float(sim_line_list[7]) >= real_ntot:
                        higher_ntot_simulated = True
                        break
                
                if not higher_ntot_simulated: # Copy evoGReg file to new directory only if no simulated ntot was higher than real ntot
                    with open(evogreg_filtered_output_filepath, "w") as evogreg_filtered_output:
                        for line in evogreg_output_lines:
                            evogreg_filtered_output.write(line)

                            
main()        
            
    