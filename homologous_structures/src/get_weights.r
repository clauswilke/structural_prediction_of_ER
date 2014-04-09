##Rscript get_weights.r <virus_pdb_dir>
##This script calculates weights from raxml tree 

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

tree <- read.tree(paste("../trees/",args[1],"/RAxML_bestTree.",args[1],sep=""))

weights <- BM.calculate.weights(tree)

write.csv(weights,file=paste("../weights/",args[1],"_weights.csv",sep=""))
