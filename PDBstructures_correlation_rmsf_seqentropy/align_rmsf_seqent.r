##plug in needed values
#rmsf <- pdb_wtdRMSF

##convert both aln positions and pdb positions vectors to numerical values
num_seq_ent_pos <- as.numeric(pdb_seqent_map[,1])
num_rmsf_pos <- as.numeric(pdb_rmsf_map[,1])

num_pdb_pos_seq_ent <- as.numeric(pdb_seqent_map[,3])
num_pdb_pos_rmsf <- as.numeric(pdb_rmsf_map[,3])

##get AA for both maps for verification with real pdb AA
pdb_aa_seq_ent <- as.vector(pdb_seqent_map[,4])
pdb_aa_rmsf <- as.vector(pdb_rmsf_map[,4])

##make data frame
seq <- seq.pdb(pdb)

na_vector <- rep(NA, length(seq))
cor_vals <- data.frame("pdb_aa" = seq, "pdb_pos" = 1:length(seq),"pdb_aa_seq_ent" = na_vector, "pdb_pos_seq_ent" = na_vector, "seq_ent" = na_vector, "pdb_aa_rmsf" = na_vector, "pdb_pos_rmsf" = na_vector, "rmsf" = na_vector) 

##map rmsf to seq ent

cor_vals$pdb_aa_seq_ent[na.omit(num_pdb_pos_seq_ent)] <- pdb_aa_seq_ent[num_seq_ent_pos[!is.na(num_pdb_pos_seq_ent)]]
cor_vals$pdb_pos_seq_ent[na.omit(num_pdb_pos_seq_ent)] <- pdb_evrate_map[,3][num_seq_ent_pos[!is.na(num_pdb_pos_seq_ent)]]
cor_vals$seq_ent[na.omit(num_pdb_pos_seq_ent)] <- seq_ent[num_seq_ent_pos[!is.na(num_pdb_pos_seq_ent)]]
cor_vals$pdb_aa_rmsf[na.omit(num_pdb_pos_rmsf)] <- pdb_aa_rmsf[num_rmsf_pos[!is.na(num_pdb_pos_rmsf)]]
cor_vals$pdb_pos_rmsf[na.omit(num_pdb_pos_rmsf)] <- pdb_aa_rmsf[num_rmsf_pos[!is.na(num_pdb_pos_rmsf)]]
cor_vals$rmsf[na.omit(num_pdb_pos_rmsf)] <- rmsf[num_rmsf_pos[!is.na(num_pdb_pos_rmsf)]]

