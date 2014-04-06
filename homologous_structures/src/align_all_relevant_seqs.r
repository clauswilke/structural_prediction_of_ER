##get all the matching pdb structures
args <- commandArgs(trailingOnly = TRUE)
library(bio3d)

get_all_relevant_seqs <- function(cutoff_blast_hits)
{
	pdb_ids <- substr(cutoff_blast_hits$subjectids,nchar(as.character(cutoff_blast_hits$subjectids))-5, nchar(as.character(cutoff_blast_hits$subjectids))-2)
	pdb_files <- get.pdb(pdb_ids, path = "~/Desktop/Research/structural_prediction_of_ER_homologous_structures/raw_pdbs")
	split.pdb(pdb_files, path = "~/Desktop/Research/structural_prediction_of_ER_homologous_structures/raw_pdbs/split_chain/")
}

find_ref_seqs <- function(pdbs,virus_pdb_dir)
{
	## find reference seq or longest sequence
	min_num_gaps = length(pdbs$ali[1,])
	for (i in 1:length(pdbs$ali[,1]))
	{       
		ali_bool <- pdbs$ali[i,] == "-"
		num_gaps <- table(ali_bool)[2]
	
		if (is.na(num_gaps)) num_gaps = 0
	
		if (num_gaps < min_num_gaps){
			min_num_gaps = num_gaps
			ref_pdb_seq <- pdbs$ali[i,]
			ref_pdb_id <- pdbs$id[i]
			
		}
	}

write(paste(">",ref_pdb_id,sep=""),paste("~/Desktop/Research/structural_prediction_of_ER_homologous_structures/alignments/reference_seqs/",virus_pdb_dir,"_ref_seq.fa",sep=""))
write(paste(ref_pdb_seq,collapse="",sep=""),paste("~/Desktop/Research/structural_prediction_of_ER_homologous_structures/alignments/reference_seqs/",virus_pdb_dir,"_ref_seq.fa",sep=""),append=TRUE)
}



cutoff_blast_hits_file <- as.character(args[2])
cutoff_blast_hits <- read.csv(paste("~/Desktop/Research/structural_prediction_of_ER_homologous_structures/blast_hits/cutoff_blast_hits/",cutoff_blast_hits_file,sep=""))
virus_pdb_dir <- as.character(args[1])
 
##download all seqs
get_all_relevant_seqs(cutoff_blast_hits)

pdb_ids_chain <- paste(substr(cutoff_blast_hits$subjectids,nchar(as.character(cutoff_blast_hits$subjectids))-5,nchar(as.character(cutoff_blast_hits$subjectids))-2),"_", substr(cutoff_blast_hits$subjectids, nchar(as.character(cutoff_blast_hits$subjectids)), nchar(as.character(cutoff_blast_hits$subjectids))),sep = "")
setwd("~/Desktop/Research/structural_prediction_of_ER_homologous_structures/raw_pdbs/split_chain/")
pdbs <- pdbaln(paste(pdb_ids_chain, ".pdb",sep = ""))
setwd("~/Desktop/Research/structural_prediction_of_ER_homologous_structures/src/")
write.fasta(pdbs, pdbs$id, pdbs$ali,file=paste("~/Desktop/Research/structural_prediction_of_ER_homologous_structures/alignments/",virus_pdb_dir,"/",virus_pdb_dir,"_all_relevant_seq.fa",sep=""))  
#pdbs <- read.fasta(paste("~/Desktop/Research/structural_prediction_of_ER_homologous_structures/alignments/",virus_pdb_dir,"/",virus_pdb_dir,"_all_relevant_seq.fa",sep="")) 

#find reference sequence, used as the first sequence for comparison
find_ref_seqs(pdbs,virus_pdb_dir)


