V.1: Amir Shahmoradi, Wilke Lab., University of Texas Austin, 11:47 Friday July 19, 2013.

This folder contains all input structures, input parameter files and some of the important (but small size) output data from "Amber" explicit solvent MD simulations.

The following pdb files are either taken directly from Daria's table, or are one of their close families:

- Hemaglutinin Precursor --> 1RD8_AB.pdb
- Crimean Congo Hemorrhagic --> 4AQF_B.pdb
- Dengue Protease Helicase --> 2JLY.pdb
- Hepatitis C Protease --> 3GOL.pdb, 3GSZ.pdb, 3I5K.pdb
- Influenza Nucleoprotein --> 4IRY.pdb
- Japanese Encephalitis Helicase --> 2Z83.pdb
- Marburg RNA Binding Domain --> 4GHA_A.pdb
- Rift Valley Fever Virus Nucleoprotein --> 3LYF.pdb

As of today, each structure has the following steps for its MS simulations by Amber, all available on Lonestar, TACC:

- A short 2000 step initial energy minimization of the structure.
- 100ps weakly-constrained heating to 300K.
- 100ps weakly-constrained density equilibration at temp~300K.
- 5ns unconstrained NPT equilibration.
- 15ns of MD trajectories.
