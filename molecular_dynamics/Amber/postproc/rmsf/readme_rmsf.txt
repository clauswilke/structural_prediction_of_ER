Amir Shahmoradi, Wilke Lab., University of Texas Austin, 3:22 PM Friday July 25, 2013.

This directory contains the (Amber cpptraj output) rmsf (root-mean-square-fluctuation per residue) files. Each file's name is detemined by the pdb_id + file_content (here, rmsf) + and the reference frame with respect to which the rmsf was calculated (inpcrd: the initial pdb coordinate file for MD simulation).

The subdirectory 'src' contains the input scripts to Amber cpptraj software that was used to calculate the rmsf values.
The rmsf values are calculated based on the positions of only CA atoms of the protein backbone. Otherwise, if all atoms in the amino acid are used, then the name of the corresponding rmsf and rmsd files have '_all' at the end.

My latest file naming convention is given in the example below:

	2FP7_B_rmsf_CA_ref_pdb.txt

		"2FP7_B" is the name of the protein.
		"rmsf" or "rmsd" indicates the type of the calculation.
		"CA" implies that the fitting and fluctuation calculations are all based on the positions of only CA atoms in the protein backbone.
		"ref_pdb" indicates that the original protein crystal structure is used as the reference for fitting of the MD trajectories.

Amir Shahmoradi, Wilke Lab., University of Texas Austin, 9:49 PM Saturday January 18, 2013.