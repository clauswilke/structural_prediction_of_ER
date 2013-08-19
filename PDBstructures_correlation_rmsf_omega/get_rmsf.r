#wts <- read.csv("virus_weights.csv")[,2]

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
		xyzwtdmeans[j] <- NA
		xyzwtdsd[j] <- NA
	}
	j = j+1
}



m_mean <- matrix(xyzwtdmeans,ncol = 3, byrow = T) 
diswtdmean <- sqrt(rowSums(m_mean^2, na.rm = F))
m_sd <- matrix(xyzwtdsd, ncol = 3, byrow = T)
diswtdsd <- sqrt(rowSums(m_sd^2, na.rm = F))
