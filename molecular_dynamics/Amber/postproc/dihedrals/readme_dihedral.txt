Amir Shahmoradi, Wednesday 2:44 AM January 23, 2014, Wilke Lab, ICMB, University of Texas Austin

	This folder contains dihedral angle (phi, psi, chi1) mean and variance values for each amino acid (AA) of each protein structure that have MD simulation trajectories. The total number of MD trajectories are 1500 frames for each structure. The first line of data, represented by frame 0, corresponds to the crystal structure taken from the protein data bank.
	The sites with non-defined dihedral angles (especially true about chi1 dihedral angles) are represented by NA as the their value.
	
	The code used to calculate dihedral angles is in subdirectory 'src'. The dihedral code, although working properly, is somewhat messy and needs improvements that I have left for future, after spending a few hours on it in vain.