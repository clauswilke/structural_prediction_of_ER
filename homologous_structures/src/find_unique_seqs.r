## establish number of residues difference to match 10% difference
args <- commandArgs(trailingOnly = TRUE)
library(bio3d)

find_unique_seqs <- function(pdbs,ref_pdb,diff)
{
	
	## initialize matrix for seq-s with enough difference
	final_aln <- matrix(" 0", nrow = 1 , ncol = length(ref_pdb$ali))
	final_aln[1,] = ref_pdb$ali

	## creating vector with sorted pdb names
	final_pdb_id <- c()
	final_pdb_id <- append(final_pdb_id,ref_pdb$id)

	for (i in 1:length(pdbs$ali[,1]))
	{
		if (length(final_pdb_id) > 1) 
		{
			count_diff = 0
			for (j in 1:length(final_aln[,1])) 
			{
				diff_bool <- (pdbs$ali[i,] != final_aln[j,] & pdbs$ali[i,] != "-" & final_aln[j,] != "-")
				n_diff <- table(diff_bool)[2]
				if (is.na(n_diff)) n_diff = 0

				if (n_diff >= diff) count_diff = count_diff+1
                
                if (count_diff == length(final_pdb_id)) {
					final_pdb_id <- append(final_pdb_id,pdbs$id[i])
					final_aln <- rbind(final_aln,pdbs$ali[i,])
                }
			}
		}
        else 
		{
			##bool vector of all places where sequence in the original aln and refseq differ 
			diff_bool <- (pdbs$ali[i,] != ref_pdb$ali & pdbs$ali[i,] != "-" & ref_pdb$ali != "-")
			n_diff <- table(diff_bool)[2]

			if (is.na(n_diff)) n_diff = 0

			if (n_diff > diff) {
				final_pdb_id <- append(final_pdb_id,pdbs$id[i])
				final_aln <- rbind(final_aln,pdbs$ali[i,])
				next
			}
        }
	}
	return(final_pdb_id)
}

virus_pdb_dir <- args[1]
ref_pdb <- read.fasta(paste("~/Desktop/Research/structural_prediction_of_ER_homologous_structures/alignments/reference_seqs/",args[2],sep=""))
all_relevant_seq_aln <- read.fasta(paste("~/Desktop/Research/structural_prediction_of_ER_homologous_structures/alignments/",virus_pdb_dir,"/",args[3],sep=""),)

percent_diff <- as.character(args[4])
if (percent_diff=='all') { diff=1
} else { 
	len_ref_seq = length(ref_pdb$ali)
	diff <- len_ref_seq*as.numeric(percent_diff)/100
}
final_pdb_id <- find_unique_seqs(all_relevant_seq_aln,ref_pdb,diff)

setwd("~/Desktop/Research/structural_prediction_of_ER_homologous_structures/raw_pdbs/split_chain/")
final_aln_fasta <- pdbaln(final_pdb_id)
setwd("~/Desktop/Research/structural_prediction_of_ER_homologous_structures/src")

if (percent_diff=='all') {
	write.fasta(final_aln_fasta,final_aln_fasta$id,final_aln_fasta$ali,file=paste("~/Desktop/Research/structural_prediction_of_ER_homologous_structures/alignments/",virus_pdb_dir,"/",virus_pdb_dir,"_all_unique_seq.fa", sep=""))  
} else {
	write.fasta(final_aln_fasta,final_aln_fasta$id,final_aln_fasta$ali,file=paste("~/Desktop/Research/structural_prediction_of_ER_homologous_structures/alignments/",virus_pdb_dir,"/",virus_pdb_dir,"_",as.character(percent_diff),"%_unique_seq.fa", sep=""))  
}
