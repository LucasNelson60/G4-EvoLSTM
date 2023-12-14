## G4-EvoLSTM
An application of EvoLSTM, a neutral evolutionary drift-simulating machine learning model developed by Dongjoon Lim and Mathieu Blanchette at McGill University, to a full-genome statistical analysis of G-quadruplex relevance.  

## INFO
`run_all_simulations.py` runs a first pass with 100 EvoLSTM simulations on every G-quadruplex-containing-region (G4CR) in the promoter regions of the human genome (hg38).  
`run_thousand_simulations.py` runs a second pass with 1000 EvoLSTM simulations on G4CRs that appear to be statistically (hence biologically) significant, with added interval-simulation outputs for each of the following ancestral sequences: _HP, _HPG, _HPGP, _HPGPN, _HPGPNRMPC, _HPGPNRMPCCS, _HPGPNRMPCCSO, _HPGPNRMPCCSOT, _HPGPNRMPCCSOTSJMCMMRHCCOOO, _HPGPNRMPCCSOTSJMCMMRHCCOOOSVCTOPBOCECFCMAOLPPEMMESC, _HPGPNRMPCCSOTSJMCMMRHCCOOOSVCTOPBOCECFCMAOLPPEMMESCLETCEOD 
`run_million_simulations.py` runs a third pass with one million EvoLSTM simulations on a handful of G4CRs (like C-MYC) that have been extensively studied.  
This is a work in progress, and although there are differences in logic between the three files (mostly pertaining to the processing of input directories/files), the goal is to merge all three into a single file.  

## Acknowledgements
evoGReg and G4_EvoLSTM (run_X_simulations): Lucas Nelson  
EvoLSTM: Dongjoon Lim and Mathieu Blanchette of the McGill Computational Genomics Lab. Link: https://github.com/DongjoonLim/EvoLSTM  
GReg: Christopher Hennecker, Lynn Yamout, Chuyang (Amos) Zhang, Chenzhi Zhao, David Hiraki, Nicolas Moitessier and Anthony Mittermaier. Link: https://github.com/Christopher-Hennecker/GReg  
Sequence alignment data: Mathieu Blanchette and the McGill Computational Genomics Lab  
