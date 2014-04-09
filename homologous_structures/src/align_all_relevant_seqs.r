##Rscript align_all_relevant_seqs.r <virus pdb directory> <${PDB_VIRUS_NAME}_cutoff_blast_hits.csv used as filtered pdb files>
##This script downloads all necessary pdb files from ${PDB_VIRUS_NAME}_cutoff_blast_hits.csv, aligns them, and outputs a fasta file in ../alignments/$PDB_VIRUS_NAME/$PDB_VIRUS_NAME_all_relevant_seq.fa. 

args <- commandArgs(trailingOnly = TRUE)
library(bio3d)

get_all_relevant_seqs <- function(cutoff_blast_hits)
{
	pdb_ids <- substr(cutoff_blast_hits$subjectids,nchar(as.character(cutoff_blast_hits$subjectids))-5, nchar(as.character(cutoff_blast_hits$subjectids))-2)
	pdb_files <- get.pdb(pdb_ids, path = "../raw_pdbs")
	split.pdb(pdb_files, path = "../raw_pdbs/split_chain/")
}

virus_pdb_dir <- as.character(args[1])
cutoff_blast_hits_file <- as.character(args[2])
cutoff_blast_hits <- read.csv(paste("../blast_hits/cutoff_blast_hits/",cutoff_blast_hits_file,sep=""))
 
##download all seqs
get_all_relevant_seqs(cutoff_blast_hits)

pdb_ids_chain <- paste(substr(cutoff_blast_hits$subjectids,nchar(as.character(cutoff_blast_hits$subjectids))-5,nchar(as.character(cutoff_blast_hits$subjectids))-2),"_", substr(cutoff_blast_hits$subjectids, nchar(as.character(cutoff_blast_hits$subjectids)), nchar(as.character(cutoff_blast_hits$subjectids))),sep = "")
if (virus_pdb_dir=="HP"){
	pdb_ids_chain <- append(pdb_ids_chain,"1RD8_AB")
}

##align all pdb structures from cutoff_blast_hits 
setwd("../raw_pdbs/split_chain/")
pdbs <- pdbaln(paste(pdb_ids_chain, ".pdb",sep = ""))
setwd("../src/")

write.fasta(pdbs, pdbs$id, pdbs$ali,file=paste("~/Desktop/Research/structural_prediction_of_ER_homologous_structures/alignments/",virus_pdb_dir,"/",virus_pdb_dir,"_all_relevant_seq.fa",sep=""))  



