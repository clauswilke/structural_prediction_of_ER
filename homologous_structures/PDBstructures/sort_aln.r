#get pdbs$ali,refseq,refpdb from align_all_blasthits.r

## establish number of residues difference to match 10% difference

len_total = length(refseq)

diff = len_total*0.02
#diff = 1
diff = round(diff, digits = 0)

## initialize matrix for seq-s with enough difference
final_aln <- matrix(" 0", nrow = 1 , ncol = length(refseq))
final_aln[1,] = refseq

## creating vector with sorted pdb names
pdbid = c()
pdbid[1] = refpdb

for (i in 1:length(pdbs$ali[,1])){
        if (length(pdbid) > 1) {
                count_diff = 0
                for (j in 1:length(final_aln[,1])) {
                        diff_bool <- (pdbs$ali[i,] != final_aln[j,] & pdbs$ali[i,] != "-" & final_aln[j,] != "-")
                        n_diff <- table(diff_bool)[2]
                        if (is.na(n_diff)) n_diff = 0

                        if (n_diff >= diff) {
                                count_diff = count_diff+1
                        }
                }
                if (count_diff == length(pdbid)){
                        pdbid[length(pdbid)+1] = pdbs$id[i]
                        final_aln <- rbind(final_aln,pdbs$ali[i,])
                }
        }
        else {
                ##bool vector of all places where sequence in the original aln and refseq differ 
                diff_bool <- (pdbs$ali[i,] != refseq & pdbs$ali[i,] != "-" & refseq != "-")
                n_diff <- table(diff_bool)[2]

                if (is.na(n_diff)) n_diff = 0

                if (n_diff > diff) {
                        pdbid[length(pdbid)+1] = pdbs$id[i]
                        final_aln <- rbind(final_aln,pdbs$ali[i,])
                        next
                }
        }

}

setwd("../raw_pdbs/split_chain/")
final_aln <- pdbaln(pdbid)
setwd("../../main_analysis_r_scripts")
