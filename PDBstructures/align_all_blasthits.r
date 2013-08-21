setwd("~/raw_pdbs/split_chain")
pdb_chains <- paste(substr(blasthits[,2],nchar(blasthits[,2])-5,nchar(blasthits[,2])-2),"_", substr(blasthits[,2], nchar(blasthits[,2]), nchar(blasthits[,2])), sep = "")
pdbs <- pdbaln(paste(pdb_chains, ".pdb",sep = ""))

## find reference seq or longest sequence

min_num_gaps = length(pdbs$ali[1,])
for (i in 1:length(pdbs$ali[,1]))
{       ali_bool <- pdbs$ali[i,] == "-"
        num_gaps <- table(ali_bool)[2]

        if (is.na(num_gaps)) num_gaps = 0

        if (num_gaps < min_num_gaps){
                min_num_gaps = num_gaps
                refseq <- pdbs$ali[i,]
                refpdb <- pdbs$id[i]
        }
}
