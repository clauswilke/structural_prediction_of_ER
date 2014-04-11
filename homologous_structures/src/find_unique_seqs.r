##Rscript find_unique_seqs.r $PDB_ID $PDB_CHAIN $PDB_VIRUS_NAME ${PDB_VIRUS_NAME}_all_relevant_seq.fa 2
##This script finds unique sequence with a given percent difference in aligned amino acids or all unique sequence with option "all"
args <- commandArgs(trailingOnly = TRUE)
library(bio3d)

find_unique_seqs <- function(pdbs,ref_pdb,ref_seq,diff) {
	
	## initialize matrix for seq-s with enough difference
	final_aln <- matrix(" 0", nrow = 1 , ncol = length(ref_seq))
	final_aln[1,] <- ref_seq

	## creating vector with sorted pdb names
	final_pdb_id <- c()
	final_pdb_id <- append(final_pdb_id,ref_pdb)

	
	for (i in 1:length(pdbs$ali[,1])) {
		
		add_seq <- TRUE
		for (j in 1:length(final_aln[,1])) {
			
			
			n_diff <- sum(pdbs$ali[i,] != final_aln[j,] & pdbs$ali[i,] != "-" & final_aln[j,] != "-")
			n_comp <- sum(pdbs$ali[i,] != "-" & final_aln[j,] != "-")
			
			if (diff == "all") { percent_diff <- 1/n_comp
			} else { percent_diff <- as.numeric(diff)/100
			}	

			if (n_diff/n_comp < percent_diff | (pdbs$id[i] %in% final_pdb_id)) {
				add_seq <- FALSE
				break
			}
		}
		if ( add_seq ){
			final_pdb_id <- append(final_pdb_id,pdbs$id[i])
			final_aln <- rbind(final_aln,pdbs$ali[i,])
		}
	
	}
	return(final_pdb_id)
}
pdb_id <- args[1]
pdb_chain <- args[2]
ref_pdb <- paste(args[1],"_",args[2],".pdb", sep="")

virus_pdb_dir <- args[3]
	
all_relevant_seq_file <- args[4]

all_relevant_seq_aln <- read.fasta(paste("../alignments/",virus_pdb_dir,"/",all_relevant_seq_file,sep=""),)
ref_seq <- all_relevant_seq_aln$ali[all_relevant_seq_aln$id==ref_pdb,]

diff <- as.character(args[5])

final_pdb_id <- find_unique_seqs(all_relevant_seq_aln,ref_pdb,ref_seq,diff)

setwd("../raw_pdbs/split_chain/")
final_aln_fasta <- pdbaln(final_pdb_id)
setwd("../../src")


if (diff=='all') {
	write.fasta(final_aln_fasta,final_aln_fasta$id,final_aln_fasta$ali,file=paste("../alignments/",virus_pdb_dir,"/",virus_pdb_dir,"_all_unique_seq.fa", sep=""))  
} else {
	write.fasta(final_aln_fasta,final_aln_fasta$id,final_aln_fasta$ali,file=paste("../alignments/",virus_pdb_dir,"/",virus_pdb_dir,"_",as.character(diff),"%_unique_seq.fa", sep=""))  
}
