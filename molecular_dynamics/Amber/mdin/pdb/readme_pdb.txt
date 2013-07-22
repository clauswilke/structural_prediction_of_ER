Amir Shahmoradi, Wilke Lab., University of Texas Austin, 5:42 PM Friday July 19, 2013.

This folder contains the pdb files (originally taken from Daria's table) prepared by AmberTools12 tleap module in order to generate the required topology and coordinate files for MD simulations.

abbreviations:

- 'I' in pdb filenames ending with '_I.pdb':
	Indicates 'Ionized' pdb files with a specific number of NaCl ions, such that the solvated pdb structures would have an ionic strength of ~150 mM.
	
- 'IS' in pdb filenames ending with '_IS.pdb':
	Indicates 'Ionized_Solvated' pdb files that have specific number of water molecules as the solvent with a specific number of NaCl ions corresponding to an ionic strength of ~150 mM. The solvation box used here is either a truncated octahedron or rectangular cuboid box. The sides of the solvation boxes are set to be at least 10 Angstroms away from the solute (protein).
