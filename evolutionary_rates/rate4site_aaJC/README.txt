Last updated by Amir Shahmoradi, Monday 12:40 AM, October 13, 2014, WilkeLab, iCMB, UT Austin


For two of the above cases (HP & HCV) the software rate4site crashes with the following error:


===================================================================================
===================================================================================

amir@linux-5xxe:~/git/structural_prediction_of_ER/evolutionary_rates> rate4site -ma -t ~/git/structural_prediction_of_ER/evolutionary_rates/trees/tree_HCV_NS5B.tre -s ~/git/structural_prediction_of_ER/evolutionary_rates/alignments/aminoacid/aln_aa_HCV_NS5B.fasta -x ~/git/structural_prediction_of_ER/evolutionary_rates/HCP.tree -y ~/git/structural_prediction_of_ER/evolutionary_rates/HCP.rates
START OF LOG FILE

 =======================================================
 the rate for site project:                             
 Version: 3.0.0.                                                                                                                                                            
 Tal Pupko and his lab:     talp@post.tau.ac.il                                                                                                                             
 Nir Ben-Tal and his lab:   bental@ashtoret.tau.ac.il                                                                                                                       
 Itay Mayrose:  itayMay@post.tau.ac.il                                                                                                                                      
 For program support, please contact Itay Mayrose       
 =======================================================


 ---------------------- THE PARAMETERS ----------------------------
tree file is: /home/amir/git/structural_prediction_of_ER/evolutionary_rates/trees/tree_HCV_NS5B.tre
seq file is: /home/amir/git/structural_prediction_of_ER/evolutionary_rates/alignments/aminoacid/aln_aa_HCV_NS5B.fasta
output file is: r4s.res
rate inference method is: empirical Bayesian estimate
using a Gamma prior distribution with: 16 discrete categories
probablistic_model is: AAJC
branch lengths optimization is ML using a gamma model

 -----------------------------------------------------------------
Optimizing branch lengths and alpha...
The tree was written to a file name called /home/amir/git/structural_prediction_of_ER/evolutionary_rates/HCP.tree
Computing the rates...
rate4site: siteSpecificRate.cpp:202: void computeEB_EXP_siteSpecificRate(int, const sequenceContainer&, const stochasticProcess&, const computePijGam&, const tree&, double&, double&, double&, double&, double, VVdouble*, unObservableData*): Assertion `sum!=0' failed.
Aborted

===================================================================================
===================================================================================
    
I  have already tried twice to run the software for these two, all in vain. At this point I am leaving these two, hoping that someone will the know the cause of the error who could help on this issue.




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
