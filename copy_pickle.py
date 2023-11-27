# AUTHOR: Lucas Nelson
import pandas as pd

GENOMES = ['hg38', '_HP', '_HPG', '_HPGP', '_HPGPN', '_HPGPNRMPC', '_HPGPNRMPCCS', '_HPGPNRMPCCSO', '_HPGPNRMPCCSOT', '_HPGPNRMPCCSOTSJMCMMRHCCOOO', '_HPGPNRMPCCSOTSJMCMMRHCCOOOSVCTOPBOCECFCMAOLPPEMMESC', '_HPGPNRMPCCSOTSJMCMMRHCCOOOSVCTOPBOCECFCMAOLPPEMMESCLETCEOD']

chrom = input("Chromosome you wish to copy over to new user.\n")
old_filepath = f'/home/mcb/users/dlim63/conservation/data/seqDictPad_chr{chrom}.pkl'
old_alignment = pd.read_pickle(old_filepath)
new_alignment = {}

for genome in GENOMES:
    new_alignment[genome] = old_alignment[genome]

new_path = f'/home/mcb/users/lnelso12/conservation/data/seqDictPad_chr{chrom}.pkl'
alignment_dataframe = pd.DataFrame.from_dict(new_alignment)
alignment_dataframe.to_pickle(new_path)


    