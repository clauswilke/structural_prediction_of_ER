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

seq <- seq.pdb(pdb)
correlated_vals <- matrix(data = NA, nrow = length(seq), ncol = 8) 

na_vector <- rep(NA, length(seq))
correlated_vals <- cbind("pdb_aa" = seq, "pdb_pos" = 1:length(seq),"pdb_aa_ev_rate" = na_vector, "pdb_pos_ev_rate" = na_vector, "ev_rate" = na_vector, "pdb_aa_rmsf" = na_vector, "pdb_pos_rmsf" = na_vector, "rmsf" = na_vector) 

##loop through each vector

#evolutionary rates vector

for (i in 1:length(num_ev_rate_pos)){
  ind_pdb <- num_pdb_pos_ev_rate[i]
  ind_ev_rate <- num_ev_rate_pos[i]

  ##check if pdb_aa from ev_rates and rmsf match the seq
  if (!(is.na(ind_pdb)) && pdb_aa_ev_rate[i] == seq[ind_pdb]) {
    ##assign values from ev_rates vector using indecies in aln positions vector to the aa in pdb using indecies from pdb positions vector
    correlated_vals[ind_pdb,3] <- pdb_aa_ev_rate[i]
    correlated_vals[ind_pdb,4] <- num_pdb_pos_ev_rate[i]
    correlated_vals[ind_pdb,5] <- ev_rate[ind_ev_rate]
  }
}

#rmsf vector

for (i in 1:length(num_rmsf_pos)){
  ind_pdb <- num_pdb_pos_rmsf[i]
  ind_rmsf <- num_rmsf_pos[i]

  if (!(is.na(ind_pdb)) && pdb_aa_rmsf[i] == seq[ind_pdb]) {
    correlated_vals[ind_pdb,6] <- pdb_aa_rmsf[i]
    correlated_vals[ind_pdb,7] <- num_pdb_pos_rmsf[i]
    correlated_vals[ind_pdb,8] <- rmsf[ind_rmsf]

  }
}

#if (!all(correlated_vals[,4] == correlated_vals[,7])) {
#        print("pdb_pos do not match")
#        print(correlated_vals[correlated_vals[,4] == correlated_vals[,7],])
#}

rmsf_ev_rate_cor <- correlated_vals
tf_rmsf <- !is.na(correlated_vals[,8])
tf_evrate <- !is.na(correlated_vals[,5])
tf <- tf_rmsf & tf_evrate
plot(rmsf_ev_rate_cor[,8][tf],rmsf_ev_rate_cor[,5][tf])

num_final_rmsf_vals <- as.numeric(correlated_vals[,8][tf])
num_final_ev_rate_vals <- as.numeric(correlated_vals[,5][tf])
r_sqrd <- cor(num_final_rmsf_vals,num_final_ev_rate_vals)
print("Correlation RMSF ev_rate:")
print(r_sqrd)
