args <- commandArgs(trailingOnly = TRUE)
library(bio3d)

wts_tbl <- read.csv(paste("~/Desktop/Research/structural_prediction_of_ER_homologous_structures/weights/",args[1],"_weights.csv",sep=""))
wts <- wts_tbl$weight

setwd("~/Desktop/Research/structural_prediction_of_ER_homologous_structures/raw_pdbs/split_chain/")
aln <- pdbaln(as.character(wts_tbl$taxon))
setwd("~/Desktop/Research/structural_prediction_of_ER_homologous_structures/src/")

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
		xyzmeans[j] <- mean(xyz[,i])
		xyzwtdmeans[j] <- sum(wts*xyz[,i])
		xyzsd[j] <- sd(xyz[,i])
		xyzwtdsd[j] <- sqrt(sum(wts*(xyz[,i]-xyzwtdmeans[j])^2))
	}
	else
	{
		xyzmeans[j] <- NA
		xyzsd[j] <- NA
	}
	j = j+1
}


m_mean <- matrix(xyzwtdmeans,ncol = 3, byrow = T) 
diswtdmean <- sqrt(rowSums(m_mean^2, na.rm = F))
m_sd <- matrix(xyzwtdsd, ncol = 3, byrow = T)
diswtdsd <- sqrt(rowSums(m_sd^2, na.rm = F))

write.csv(diswtdsd,paste("~/Desktop/Research/structural_prediction_of_ER_homologous_structures/rmsf/",args[2],sep=""))

##rmsf for non-weighted calculation
#m_mean <- matrix(xyzmeans,ncol = 3, byrow = T)
#dis_mean <- sqrt(rowSums(m_mean^2, na.rm = F))
#m_sd <- matrix(xyzsd, ncol = 3, byrow = T)
#dis_sd <- sqrt(rowSums(m_sd^2, na.rm = F))

