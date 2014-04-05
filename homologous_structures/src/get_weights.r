args <- commandArgs(trailingOnly = TRUE)
library(bio3d)
library(ape)

BM.executable = "./BranchManager-highres.jar" # executable for BranchManager (BM) program
tmp.tree = "tree.FHDieck239swl94mv.txt" # temporary tree file for use with external BM program
tmp.output = "BM.wle238LKswo209sal.out" # temporary output file for use with external BM program

BM.calculate.weights <- function( phylo ) # invokes branch manager to calculate weights for the given tree. Works only on unix.
{
	phylo$node.label<-NULL # remove node labels since BM cannot deal with them
	write.tree( phylo, tmp.tree )
	result <- system( paste( "java -cp", BM.executable, "BM", tmp.tree, ">", tmp.output ), intern=F )
	result <- read.table( tmp.output, header=F )
	colnames( result )<-c("taxon","weight")
	unlink( tmp.tree )
	unlink( tmp.output)
	result
}

tree <- read.tree(paste("~/Desktop/Research/structural_prediction_of_ER_homologous_structures/trees/best_tree/RAxML_bestTree.",args[1],sep=""))

weights <- BM.calculate.weights(tree)

write.csv(weights,file=paste("~/Desktop/Research/structural_prediction_of_ER_homologous_structures/weights/",args[1],"_weights.csv",sep=""))
