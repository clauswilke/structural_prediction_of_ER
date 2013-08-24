Written by Eleisha Jackson, August 3, 2013

Description:
Entropy for Designed Sequences 
Each PDB found in /structural_prediction_of_ER/molecular_dynamics/Amber/postproc/pdb was used a template to design 500 proteins using Rosetta
The proteins were designed using a flexible backbone approach "Backrub" at a temperature of 0.6. 
Due to size constraints the PDB files are not in this folder. If needed, please ask.
Note: To run get_entropies you will need NUMPY installed. 
Contents of Directories:

scripts/
	extract_sequences.py
		Script that extracts the sequences from the designed PDBS and creates the sequences found in designed_sequences/
	get_entropies.py
		Script that calculates entropies from the designed sequences.
	cMyPDBParser.py
		This is a helper script used to extract the sequences.
	rename_sequences.py 
		This is a script that just renamed the .pdb returned from Rosetta. 
	
designed_sequences/
	Contains the sequences for the designed proteins (500 for each of 10 PDBS = 5000 Sequences Total)
		
entropies/
	Contains all entropies calculated for the designed protein sequences  
	
