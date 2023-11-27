# AUTHOR: Lucas Nelson
import os

def return_sorted_filenames(chrom):
    
    chrom_filepath = f"/home/mcb/users/lnelso12/evoGReg/outputs/chr{chrom}"
    
    existing_genes = set() # Ensures only files that are already in the directory are checked from the TSS file
    for existing_gene in os.listdir(chrom_filepath):
        existing_genes.add(existing_gene)
        
    tss_filepath = f"chr{chrom}_tss_indices.bed"
    filenames = []

    with open(tss_filepath, "r") as tss_file:
        gene_infos = [gene_info.strip("\n").split("\t") for gene_info in tss_file.readlines()]
        
        for gene_info in gene_infos:
            
            gene_name = gene_info[3]
            main_filename = f"{gene_name}.csv"
            comp_filename = f"{gene_name}_comp.csv"
            
            if main_filename in existing_genes:
                filenames.append(main_filename)
            if comp_filename in existing_genes:
                filenames.append(comp_filename)
                
    
    sorted_filenames = sorted(filenames, key=lambda x: x.split("_")[0])
    return sorted_filenames


def main():
    
    for chrom in range(3,22):
        
        chrom = str(chrom)
        sorted_filenames = return_sorted_filenames(chrom)
        prev_filename = sorted_filenames[0]
        
        for filename in sorted_filenames[2::2]:
            split_filename = filename.split("_")[0]
            split_prev_filename = prev_filename.split("_")[0]
            csv_stripped_filename = filename.split(".")[0]
            
            main_filepath = f"/home/mcb/users/lnelso12/evoGReg/outputs/chr{chrom}/{filename}"
            comp_filepath = f"/home/mcb/users/lnelso12/evoGReg/outputs/chr{chrom}/{csv_stripped_filename}_comp.csv"
            
            if split_filename == split_prev_filename:
                os.remove(main_filepath)
                os.remove(comp_filepath)
                print(f"Removed file {filename}")
                print(f"Removed file {csv_stripped_filename}_comp.csv")
            
            prev_filename = filename
                



main()
