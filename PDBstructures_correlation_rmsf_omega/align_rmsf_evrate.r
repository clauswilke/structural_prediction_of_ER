##plug in needed values
#ev_rate <- pdb_alpha_beta_table[,3]/pdb_alpha_beta_table[,2]
#rmsf <- pdbRMSFwtd

##convert both aln positions and pdb positions vectors to numerical values
num_ev_rate_pos <- as.numeric(pdb_evrate_map[,1])
num_rmsf_pos <- as.numeric(pdb_rmsf_map[,1])

num_pdb_pos_ev_rate <- as.numeric(pdb_evrate_map[,3])
num_pdb_pos_rmsf <- as.numeric(pdb_rmsf_map[,3])

##get AA for both maps for verification with real pdb AA
pdb_aa_ev_rate <- as.vector(pdb_evrate_map[,4])
pdb_aa_rmsf <- as.vector(pdb_rmsf_map[,4])

##make data frame
seq <- seq.pdb(pdb)

na_vector <- rep(NA, length(seq))
cor_vals <- data.frame("pdb_aa" = seq, "pdb_pos" = 1:length(seq),"pdb_aa_ev_rate" = na_vector, "pdb_pos_ev_rate" = na_vector, "ev_rate" = na_vector, "pdb_aa_rmsf" = na_vector, "pdb_pos_rmsf" = na_vector, "rmsf" = na_vector) 

##map rmsf to ev rate

cor_vals$pdb_aa_ev_rate[na.omit(num_pdb_pos_ev_rate)] <- pdb_aa_ev_rate[num_ev_rate_pos[!is.na(num_pdb_pos_ev_rate)]]
cor_vals$pdb_pos_ev_rate[na.omit(num_pdb_pos_ev_rate)] <- pdb_evrate_map[,3][num_ev_rate_pos[!is.na(num_pdb_pos_ev_rate)]]
cor_vals$ev_rate[na.omit(num_pdb_pos_ev_rate)] <- ev_rate[num_ev_rate_pos[!is.na(num_pdb_pos_ev_rate)]]
cor_vals$pdb_aa_rmsf[na.omit(num_pdb_pos_rmsf)] <- pdb_aa_rmsf[num_rmsf_pos[!is.na(num_pdb_pos_rmsf)]]
cor_vals$pdb_pos_rmsf[na.omit(num_pdb_pos_rmsf)] <- pdb_aa_rmsf[num_rmsf_pos[!is.na(num_pdb_pos_rmsf)]]
cor_vals$rmsf[na.omit(num_pdb_pos_rmsf)] <- rmsf[num_rmsf_pos[!is.na(num_pdb_pos_rmsf)]]

