##Rscript get_pdb_availability.r <pdb_id> <pdb_chain> <virus_pdb_dir>
##This script makes a pdb_availability table. The output is written to ../pdb_availability/pdb_availability.txt. 
args <- commandArgs(trailingOnly = TRUE)
library(bio3d)

pdb_availability <- c()

pdb_id <- args[1]
pdb_chain <- args[2]
pdb <- paste(pdb_id,"_",pdb_chain,sep="")
virus_pdb_dir <- args[3]

pdb_availability <- append(pdb_availability,virus_pdb_dir)
pdb_availability <- append(pdb_availability,pdb_id)

all_blast_hits <- read.csv(paste("../blast_hits/all_blast_hits/",virus_pdb_dir,"_all_blast_hits.csv",sep=""))
cutoff_blast_hits <- read.csv(paste("../blast_hits/cutoff_blast_hits/",virus_pdb_dir,"_cutoff_blast_hits.csv",sep=""))
all_unique_aln <- read.fasta(paste("../alignments/",virus_pdb_dir,"/",virus_pdb_dir,"_all_unique_seq.fa",sep=""))
ten_percent_unique_aln <- read.fasta(paste("../alignments/",virus_pdb_dir,"/",virus_pdb_dir,"_10%_unique_seq.fa",sep=""))
five_percent_unique_aln <- read.fasta(paste("../alignments/",virus_pdb_dir,"/",virus_pdb_dir,"_5%_unique_seq.fa",sep=""))  
two_percent_unique_aln <- read.fasta(paste("../alignments/",virus_pdb_dir,"/",virus_pdb_dir,"_2%_unique_seq.fa",sep=""))  

pdb_availability <- append(pdb_availability,length(all_blast_hits$subjectids))
pdb_availability <- append(pdb_availability,length(cutoff_blast_hits$subjectids))
pdb_availability <- append(pdb_availability,length(all_unique_aln$id))
pdb_availability <- append(pdb_availability,length(ten_percent_unique_aln$id))
pdb_availability <- append(pdb_availability,length(five_percent_unique_aln$id))
pdb_availability <- append(pdb_availability,length(two_percent_unique_aln$id))


header <- c("virus","pdb_id","all_available_seqs_blast","cutoff_seqs_blast","all_unique_seqs","seqs_with_10%_diff","seqs_with_5%_diff","seqs_with_2%_diff")
write(header,file=paste("../pdb_availability/",pdb,"_pdb_availability.txt",sep=""),ncolumns=8 ,sep="\t")
write(as.character(pdb_availability), file=paste("../pdb_availability/",pdb,"_pdb_availability.txt",sep=""),ncolumns=8 , sep="\t",append=TRUE)





