##Rscript get_rmsf.r <pdb_id> <pdb_chain> <virus_pdb_dir> <${PDB_VIRUS_NAME}_wtd_rmsf.txt file to be written>
##This script calculates rmsf from ${PDB_VIRUS_NAME}_5%_unique_seq.fa file. 
args <- commandArgs(trailingOnly = TRUE)
library(bio3d)

pdb_id <- args[1]
pdb_chain <- args[2]
virus_pdb_dir <- args[3]
wts_tbl <- read.csv(paste("../weights/",virus_pdb_dir,"_weights.csv",sep=""))

aln <- read.fasta(paste("../alignments/",virus_pdb_dir,"/",virus_pdb_dir,"_5%_unique_seq.fa", sep=""))

##match weight values order to pdbs in the alignment
wts_tbl <- wts_tbl[match(aln$id, wts_tbl$taxon),]
wts <- wts_tbl$weight

##realign the pdbs to get xyz coordinates of C alphas
setwd("../raw_pdbs/split_chain/")
aln <- pdbaln(as.character(wts_tbl$taxon))
setwd("../../src/")

gaps.xyz <- gap.inspect(aln$xyz)

##fit structures based on all non-gap positions
xyz <- fit.xyz( fixed = aln$xyz[1,],mobile = aln$xyz, fixed.inds = gaps.xyz$f.inds, mobile.inds = gaps.xyz$f.inds)

xyzmeans <- c()
xyzwtdmeans <- c()
xyzsd <- c()
xyzwtdsd <- c()

##the loop calculates weighted variance in the location of the new x,y,z positions produced by fit.xyz()
j = 1
for (i in 1:length(xyz[1,])){
	if (all(!is.na(xyz[,i]))) 
	{
		#xyzmeans[j] <- mean(xyz[,i])
		xyzwtdmeans[j] <- sum(wts*xyz[,i])
		xyzsd[j] <- sd(xyz[,i])
		xyzwtdsd[j] <- sqrt(sum(wts*(xyz[,i]-xyzwtdmeans[j])^2))
	}
	else
	{
		xyzwtdmeans[j] <- NA
		xyzwtdsd[j] <- NA
	}
	j = j+1
}

##calculate rmsf
m_mean <- matrix(xyzwtdmeans,ncol = 3, byrow = T) 
diswtdmean <- sqrt(rowSums(m_mean^2, na.rm = F))
m_sd <- matrix(xyzwtdsd, ncol = 3, byrow = T)
diswtdsd <- sqrt(rowSums(m_sd^2, na.rm = F))

##output a formated file
pdb <- paste(pdb_id,"_",pdb_chain,".pdb",sep="")

res_id <- aln$ali[aln$id==pdb,]
print(aln$id==pdb)
print(aln$ali)

rmsf <- diswtdsd[!(res_id=="-")]
res_id_no_gaps <- res_id[!(res_id=="-")]

rmsf_data <- data.frame(res_num=c(1:length(res_id_no_gaps)), res_id=res_id_no_gaps, rmsf=rmsf)

write.table(rmsf_data,paste("../rmsf/",args[4],sep=""),quote=F,sep="\t",row.names=F)


