Written by Amir Shahmoradi, Saturday 3:18 PM, October 11, 2014, WilkeLab, iCMB, UT Austin

This directory contains the evolutionary rates calculated using the aaJC method of the software rate4site (the slow version of the code compiled through Makefile_slow). Example command line arguments that are used to generate the rates:

	rate4site -ma -t ~/git/structural_prediction_of_ER/evolutionary_rates/trees/tree_H1N1_HA.tre -s ~/git/structural_prediction_of_ER/evolutionary_rates/alignments/aminoacid/aln_aa_H1N1_HA.fasta -x ~/git/structural_prediction_of_ER/evolutionary_rates/rate4site_aaJC/HP.tree -y ~/git/structural_prediction_of_ER/evolutionary_rates/rate4site_aaJC/HP.rates

    rate4site -ma -t ~/git/structural_prediction_of_ER/evolutionary_rates/trees/tree_crimean_congo.Nucleoprotein.tre -s ~/git/structural_prediction_of_ER/evolutionary_rates/alignments/aminoacid/aln_aa_crimean_congo.Nucleoprotein.fasta -x ~/git/structural_prediction_of_ER/evolutionary_rates/rate4site_aaJC/CCHFN.tree -y ~/git/structural_prediction_of_ER/evolutionary_rates/rate4site_aaJC/CCHFN.rates
    
    rate4site -ma -t ~/git/structural_prediction_of_ER/evolutionary_rates/trees/tree_H1N1_NP.tre -s ~/git/structural_prediction_of_ER/evolutionary_rates/alignments/aminoacid/aln_aa_H1N1_NP.fasta -x ~/git/structural_prediction_of_ER/evolutionary_rates/rate4site_aaJC/INP.tree -y ~/git/structural_prediction_of_ER/evolutionary_rates/rate4site_aaJC/INP.rates
    
    rate4site -ma -t ~/git/structural_prediction_of_ER/evolutionary_rates/trees/tree_dengue_ns3.tre -s ~/git/structural_prediction_of_ER/evolutionary_rates/alignments/aminoacid/aln_aa_dengue_ns3.fasta -x ~/git/structural_prediction_of_ER/evolutionary_rates/rate4site_aaJC/DPH.tree -y ~/git/structural_prediction_of_ER/evolutionary_rates/rate4site_aaJC/DPH.rates
    
    rate4site -ma -t ~/git/structural_prediction_of_ER/evolutionary_rates/trees/tree_HCV_NS5B.tre -s ~/git/structural_prediction_of_ER/evolutionary_rates/alignments/aminoacid/aln_aa_HCV_NS5B.fasta -x ~/git/structural_prediction_of_ER/evolutionary_rates/rate4site_aaJC/HCP.tree -y ~/git/structural_prediction_of_ER/evolutionary_rates/rate4site_aaJC/HCP.rates
    
    rate4site -ma -t ~/git/structural_prediction_of_ER/evolutionary_rates/trees/tree_JEV.tre -s ~/git/structural_prediction_of_ER/evolutionary_rates/alignments/aminoacid/aln_aa_JEV.fasta -x ~/git/structural_prediction_of_ER/evolutionary_rates/rate4site_aaJC/JEHN.tree -y ~/git/structural_prediction_of_ER/evolutionary_rates/rate4site_aaJC/JEHN.rates
    
    rate4site -ma -t ~/git/structural_prediction_of_ER/evolutionary_rates/trees/tree_marburg.vp35.tre -s ~/git/structural_prediction_of_ER/evolutionary_rates/alignments/aminoacid/aln_aa_marburg.vp35.fasta -x ~/git/structural_prediction_of_ER/evolutionary_rates/rate4site_aaJC/MRNABD.tree -y ~/git/structural_prediction_of_ER/evolutionary_rates/rate4site_aaJC/MRNABD.rates
    
    rate4site -ma -t ~/git/structural_prediction_of_ER/evolutionary_rates/trees/tree_riftvalley.tre -s ~/git/structural_prediction_of_ER/evolutionary_rates/alignments/aminoacid/aln_aa_riftvalley.fasta -x ~/git/structural_prediction_of_ER/evolutionary_rates/rate4site_aaJC/RVFVNP.tree -y ~/git/structural_prediction_of_ER/evolutionary_rates/rate4site_aaJC/RVFVNP.rates
    
    rate4site -ma -t ~/git/structural_prediction_of_ER/evolutionary_rates/trees/tree_westnile_chainb.tre -s ~/git/structural_prediction_of_ER/evolutionary_rates/alignments/aminoacid/aln_aa_westnile_chainb.fasta -x ~/git/structural_prediction_of_ER/evolutionary_rates/rate4site_aaJC/WNPB.tree -y ~/git/structural_prediction_of_ER/evolutionary_rates/rate4site_aaJC/WNPB.rates
    
    
    
The naming convention for the output tree and rates files of each of the structures follows those of the README file in the root directory of the repository.
