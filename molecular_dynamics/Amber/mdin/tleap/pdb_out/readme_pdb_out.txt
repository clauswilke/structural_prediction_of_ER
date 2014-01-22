Amir Shahmoradi, Tuesday 7:45 PM January 21, 2014, Wilke Lab, ICMB, University of Texas Austin

	This folder contains the output pdb files from AmberTools12 TLEAP module.

	abbreviations:
	
	- 'amber' in pdb file names ending with '_amber.pdb':
		Indicates pdb files that have been corrected for their missing atoms in their original crystal pdb structures.
			ATTN:
				The original input crystal pdb files should NOT be used in any future post-processing analyses, since they lack the atoms that have been added to the structures by TLEAP for MD simulations. The MD trajectories correspond to '*_amber.pdb' atoms.
			ATTN:
				Amber TLEAP converts some of the amino acid three-letter names according to Amber forcefield conventions for internal use in MD simulations.
				It also removes the chain ID letters from the pdb files.
				Therefore, for post processing purposes, '*_amber.pdb' should be first converted to the orginal naming convention and ad hoc chain IDs should be added to these pdbs in order to make them compatible with external software such as DSSP or VMD. The corrected pdb files are named '*_bres.pdb' & '*_X.pdb' and are available in the directory 'postproc'.

	- 'I' in pdb file names ending with '_I.pdb':
		Indicates 'Ionized' pdb files with a specific number of NaCl ions, such that the solvated pdb structures would have an ionic strength of ~150 mM.
			ATTN:
				These files are not included here in Git repository since they are of no use in future analyses.
		
	- 'IS' in pdb file names ending with '_IS.pdb':
		Indicates 'Ionized_Solvated' pdb files that have specific number of water molecules as the solvent with a specific number of NaCl ions corresponding to an ionic strength of ~150 mM. The solvation box used here is either a truncated octahedron or rectangular cuboid box. The sides of the solvation boxes are set to be at least 10 Angstroms away from the solute (protein).
			ATTN:
				These files are not included here in Git repository since they are of no use in future analyses and take up a significant amount of memory space.