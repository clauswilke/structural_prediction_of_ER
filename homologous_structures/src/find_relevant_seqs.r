##Rscript find_relevant_seqs.r <pdb_id> <pdb_chaind> <all_blast_hits.csv to written>
##find all sequences similar to a given chain in a protein using BLAST.
args <- commandArgs(trailingOnly = TRUE)

get_blasthits <- function(pdb,chain) 
{	 
	library(bio3d)
	
##download the pdb file
	pdb <- as.character(pdb)	
	get.pdb(pdb,path="~/Desktop/Research/structural_prediction_of_ER_homologous_structures/raw_pdbs")
	
##get the given chain sequence
	pdb <- paste(substr(pdb,1,4),".pdb",sep="")
	
##split the pdb into chain structures
	split.pdb(paste("~/Desktop/Research/structural_prediction_of_ER_homologous_structures/raw_pdbs/",pdb,sep=""), "~/Desktop/Research/structural_prediction_of_ER_homologous_structures/raw_pdbs/")
	
##establish a given chain
	chain <- as.character(toupper(chain))
	
	read_chain <- paste(substr(pdb,1,4),"_", chain, ".pdb",sep="")
	
##Blast according to the given chain
	pdb_chain <- read.pdb(paste("~/Desktop/Research/structural_prediction_of_ER_homologous_structures/raw_pdbs/split_chain/",read_chain,sep="")) 
	
	chain_seq <- seq.pdb(pdb_chain)
	blast <- blast.pdb(chain_seq)
	
##plot the blast to find better cutoff. x-axis: similar sequence No; y-axis: e-value,bitscore,sequence similarity, etc. 
	#plot <- plot.blast(blast) 
	
	hit_table_all <- blast$hit.tbl
	
##throw away hits that are 100% matching
	mismatches <- hit_table_all[,6]
	
	zero_mismatches <- 0 != as.numeric(mismatches)
	zero_mismatches[1] = TRUE ##keep the original first sequence in BLAST
	
	hit_table <- hit_table_all[zero_mismatches,]
	
	return(hit_table_all)
}

##requesting pdb structure name and the chain to be blasted
pdb <- args[1]
chain <- args[2]

all_blast_hits <- get_blasthits(pdb, chain)
all_blast_hits_file <- as.character(args[3])

write.csv(all_blast_hits,paste("~/Desktop/Research/structural_prediction_of_ER_homologous_structures/all_blast_hits/",all_blast_hits_file,sep=""))


















































