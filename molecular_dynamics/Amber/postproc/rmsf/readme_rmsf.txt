Amir Shahmoradi, Tuesday 10:51 PM January 21, 2014, Wilke Lab, ICMB, University of Texas Austin

	This directory contains the (Amber CPPTRAJ output) mean and standard_deviation of rmsf (root-mean-square-fluctuation per residue) and rmsd (root-mean-square-distance of all CA atoms from the original pdb structure) files. The convention for file naming is the following:
	
		<pdb name> _ <chain id> _ <reference structure to which fitting was made (Cpdb stands for Crystal pdb)> _ <the backbone atom used for fluctuation calculation (CA)> .rmsf
		<pdb name> _ <chain id> _ <reference structure to which fitting was made (Cpdb stands for Crystal pdb)> _ <the backbone atom used for fluctuation calculation (CA)> .rmsd
	
	The subdirectory 'src' contains the input scripts to Amber CPPTRAJ software used to calculate the rmsf values. First the mdcrd files from MD simulations for each pdb structure are corrected and combined into a single '*_image.mdrcd' files. The necessary scripts are in the subdirectory 'image'.

	Then the rmsf values are calculated using the CPPTRAJ scripts in subdirectory 'rmsf'.
