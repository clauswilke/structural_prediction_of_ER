##plug in needed values
pdb <- read.pdb("3GOL_X_md.renumbered.pdb")
pdb_rmsf_map <- read.csv("HepC_rmsf_map.txt", sep = "\t")
md_rmsf_map <- read.table("3GOL_rmsf_ref_inpcrd_md.txt")
md_rmsf <- as.numeric(md_rmsf_map[,3])
pdb_rmsf <- as.numeric(read.csv("HepC_wtdRMSF.csv")[,2]) 

##convert both aln positions and pdb position vectors to numerical values
num_rmsf_pos <- as.numeric(pdb_rmsf_map[,1])
num_pdb_pos_rmsf <- as.numeric(pdb_rmsf_map[,3])
num_md_pos_rmsf <- as.numeric(md_rmsf_map[,1])

##get AA for rmsf for verification with real pdb AA
pdb_aa_rmsf <- as.vector(pdb_rmsf_map[,4])

##make data frame
seq <- seq.pdb(pdb)

na_vector <- rep(NA, length(seq))
cor_vals <- data.frame("pdb_aa" = seq, "pdb_pos" = 1:length(seq), "md_rmsf_pos" = na_vector,"md_rmsf" = na_vector, "pdb_aa_rmsf" = na_vector, "pdb_pos_rmsf" = na_vector, "pdb_rmsf" = na_vector) 

##map pdb rmsf to md rmsf

cor_vals$md_rmsf_pos[na.omit(num_md_pos_rmsf)] <- md_rmsf_map[,1]
cor_vals$md_rmsf[na.omit(num_md_pos_rmsf)] <- md_rmsf
cor_vals$pdb_aa_rmsf[na.omit(num_pdb_pos_rmsf)] <- pdb_aa_rmsf[num_rmsf_pos[!is.na(num_pdb_pos_rmsf)]]
cor_vals$pdb_pos_rmsf[na.omit(num_pdb_pos_rmsf)] <- num_pdb_pos_rmsf[num_rmsf_pos[!is.na(num_pdb_pos_rmsf)]]
cor_vals$pdb_rmsf[na.omit(num_pdb_pos_rmsf)] <- pdb_rmsf[num_rmsf_pos[!is.na(num_pdb_pos_rmsf)]]

write.csv(cor_vals, "HepC_rmsf_pdb_rmsf_md_corr_vals.csv")

md_rmsf_no_na <- cor_vals$md_rmsf[!is.na(cor_vals$pdb_rmsf)]
pdb_rmsf_no_na <- cor_vals$pdb_rmsf[!is.na(cor_vals$pdb_rmsf)]
cor.test(md_rmsf_no_na, pdb_rmsf_no_na)
cor.test(md_rmsf_no_na, pdb_rmsf_no_na, method = "spearman")
