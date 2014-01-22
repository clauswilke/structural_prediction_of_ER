Amir Shahmoradi, Tuesday 7:20 PM January 21, 2014, Wilke Lab, ICMB, University of Texas Austin

	This directory contains all input crystal structures modified to contain only the desired chains in the subdirectory 'pdb_in'.
	These input pdb files are first modified by Amber TLEAP preprocessing package in order to include all missing atoms in the crystal pdb files.
	The corrected pdb files are output by TLEAP with a file name of the format '*_amber.pdb'.
	ATTN: TLEAP removes the chain ID from the crystal pdb files and modifies the three-letter names of some amino acids to make them consistent with Amber forcefield.
		  Therefore, for post-processing purposes the TLEAP output pdb files must be first corrected for their lack of chain ID and different amino acid names.
		  The amino acid name corrections can be done using Amber software ambpdb on TACC Lonestar command line:
				
				Example ambpdb command format:
				
					ambpdb -bres -p 1RD8_I_cpptraj.prmtop <$WORK/setup/1RD8_I.inpcrd> 1RD8_I_amber.pdb
					ambpdb -bres -p 2FP7_B_cpptraj.prmtop <$WORK/setup/2FP7_B.inpcrd> 2FP7_B_bres.pdb

		  An ad hoc chain ID can be added to the resulting '*_bres.pdb' files using the codes that I have written available in directory 'postproc' in Amber directory.
		  
	In addition, TLEAP also generates pdb files that contain a specific number of water molecules such that each atom of the protein is at least 10 Angstroms away from the periodic box walls. The system is then neutralized by adding appropriate number of ions. Then given the number of water molecules added, a specific number of NaCl ions corresponding to an ionic strength of ~150 mMol is added to the system. The corresponding files are named '*_IS.pdb' in which 'I' stands for being ionized and 'S' stands for being solvated in water molecules. The pdb files with name '*_I.pdb' are the same as '*_IS.pdb' except that they don't contain water molecules.
	
	To run TLEAP software type the following on TACC command line:
	
		tleap -s -f $AMBERHOME/dat/leap/cmd/leaprc.ff12SB
		
	If TLEAP is not recognized, it should be in the following directory:
	
		$AMBERHOME/AmberTools/bin/tleap
		
	Once in TLEAP, the following commands will help modify the pdb files and generate the required files in order to be used in Amber MD simulation:

		source leaprc.ff12SB
			! load forcefield
			
		loadamberparams frcmod.ionsjc_tip3p
			! needed for adding counterions. Taken from : http://archive.ambermd.org/201105/0704.html
			
		com = loadpdb *.pdb
			!	load the pdb file
			
		comsol = copy com
			! This will be the solvated + ionized pdb
			
		solvateoct comsol TIP3PBOX 10.0
			! Add solvent (water) molecules using an octahedron cell shape. This will minimize the number of solvent molecules needed for globular proteins.
			
		solvatebox comsol TIP3PBOX 10.0
			! Add solvent (water) molecules using an cubic cell shape. This will minimize the number of solvent molecules needed for elongated proteins.
			
		check com
			! see what sign the total charge is. For example, 3GOL_A.pdb is 15. So, Cl- ions are needed to neutralize the system.
		
		check comsol
			! see what sign the total charge is. For example, 3GOL_A.pdb is 15. So, Cl- ions are needed to neutralize the system.
			
		addions com Cl- 0
		addions com Na+ 0
			! Neutralize the non-solvated system by adding Cl- or Na+ ions whichever needed.
			
		addionsrand comsol Cl- 0
		addionsrand comsol Na+ 0
			! Neutralize the solvated system by adding Cl- or Na+ ions whichever needed.
			
		charge com
		charge comsol
			! check the systems' total charge status (should give zero after adding ions).
			
		check com
		check comsol
			! check for any potential errors in the pdb files. The following error can be corrected by the command 'deletebond'.
				WARNING: There is a bond of 5.625991 angstroms between: 
				-------  .R<THR 97>.A<C 13> and .R<SER 98>.A<N 1>
				deletebond com.97.C com.98.N
		
		savepdb com *_I.pdb
		savepdb comsol *_IS.pdb
			! Save the modified pdb files.
		
		saveamberparm com *_I.prmtop _I.inpcrd
		saveamberparm comsol *_IS.prmtop *_IS.inpcrd
			! Save the parameter topology and input coordinate files for Amber MD simulations.
