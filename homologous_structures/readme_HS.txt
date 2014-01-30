Amir Shahmoradi, Wednesday 10:13 PM, Jan 28 2014, WilkeLab, ICMB, University of Texas Austin

This directory contains primarily the works of Daria Sydkova on homologous structures, such as common-residue rmsf among homologous structures in a protein family. All works of Daria are collected under directory '#Daria'. At present, due to the lack of a sufficient number of pdb crystal structures in the same protein families with high sequence similarity, Daria has considered only four protein families, for which there are also MD simulations available for one of the family's representative structures (1RD8_AB, 2FP7_B, 2Z83_A, 3GOL_A).

Other directories:
	
	alignments:
		contains the structural alignments done by Daria (or someone else?) that I arranged and renamed and restored in this separate directory.
		
	maps:
		contains mapping between the amino acids (AAs) of the homologous structures and the AAs of the representative structure for which MD simulations are available.
		
	rmsf:
		contains the rmsf values from structural alignments, corresponding to the existing residues in the representative structures.
		The naming convention '*_HS.rmsf' here means rmsf values from Homologous Structures, mapped to the representative structure. Here * indicates the name of the representative structure.
		First column of data in these files 'res_num' is the residue number in the representative structure.
		Second column of data in these files 'res_name' is the residue name in the representative structure.
		Third column of data in these files 'rmsf' is the site-specific rmsf value measured using homologous structures.

	src:
		should ideally contain the source codes that were used to generate these data. However, at the moment only Daria can help on how these data were generated.

