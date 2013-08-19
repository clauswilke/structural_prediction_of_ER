#go through each column of the alignment 
#calculate the occurence of each residue at this column 
#add up and calculate all the occurences 
seq_ent <- c()

for (j in 1:length(aln$ali[1,])){
  
  if (all(aln$ali[,j] != "-")){ 
    t = table(aln$ali[,j])
    total = sum(t)
    rel_occurence = t/total 
    ent = -sum(rel_occurence*log(rel_occurence))
    seq_ent[j] <- ent
  }
  else seq_ent[j] <- NA
}
