Amir Shahmoradi, Wilke Lab., University of Texas Austin, Monday 5:47 PM Nov 18, 2013.

This folder contains Relative Solvent Accessibility (rsa) for each amino acid (aa) of each protein structure that have MD simulation trajectories. The total number of MD trajectories are 1500 frames for each structure.


Amir Shahmoradi, Monday 4:00 PM Jan 20, 2013, Wilke Lab., University of Texas Austin
	
	Each row in every *.rsa file represents the rsa of the protein amino acids for the MD trajectory frame which is determined by its number in the first column, followed by the rsa values for each amino acid.
	The rsa files are obtained  using the Fortran (2008) source codes that are in the subdirectory 'src'. This is done by normalizing the asa values to the corresponding maximum ASA values obtained in Tien et al (2013) paper.
	
	To compile the Fortran code and create the executable on TACC, use the following command:
	
		> ifort main_dynamic_rsa.f90 -o main_dynamic_rsa.o

	To run the code type:
	
		> main_dynamic_rsa.o
		
	The code will automatically inform the user about the required command-line input file names in order to generate the rsa values.
	This code also calculates the summary RSA files and outputs them in file names ending with '*_sum.rsa'.  The first two moments (mean and variance and median) of the rsa values from MD trajectories for each protein residue are summarized in the corresponding summary files in the subdirectory 'rsa_summaries'.
