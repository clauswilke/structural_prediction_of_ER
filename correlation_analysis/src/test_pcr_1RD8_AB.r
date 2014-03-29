library( pls )

# Replace NAs with mean values in data:
  pcr.desent      = data_1RD8_AB$desent
 #pcr.rsa_cr      = data_1RD8_AB$rsa_cr
  pcr.bfca		  = data_1RD8_AB$bfca
  pcr.rsa_avg_md  = data_1RD8_AB$rsa_avg_md
  pcr.rsa_var_md  = data_1RD8_AB$rsa_var_md
  pcr.phi_var_md  = data_1RD8_AB$phi_var_md
  pcr.psi_var_md  = data_1RD8_AB$psi_var_md
  pcr.chi1_var_md = data_1RD8_AB$chi1_var_md
  pcr.rmsf_avg_md = data_1RD8_AB$rmsf_avg_md
  pcr.rmsf_std_md = data_1RD8_AB$rmsf_std_md
  pcr.iwcn_avg_md = 1./data_1RD8_AB$wcn_avg_md
  
  pcr.bfca[is.na(pcr.bfca)]               = mean(na.omit(pcr.bfca))
  pcr.rsa_avg_md[is.na(pcr.rsa_avg_md)]   = mean(na.omit(pcr.rsa_avg_md))
  pcr.rsa_var_md[is.na(pcr.rsa_var_md)]   = mean(na.omit(pcr.rsa_var_md))
  pcr.phi_var_md[is.na(pcr.phi_var_md)]   = mean(na.omit(pcr.phi_var_md))
  pcr.psi_var_md[is.na(pcr.psi_var_md)]   = mean(na.omit(pcr.psi_var_md))
  pcr.chi1_var_md[is.na(pcr.chi1_var_md)] = mean(na.omit(pcr.chi1_var_md))

data_mna_1RD8_AB = data_1RD8_AB[is.na(data_1RD8_AB)]=mean(na.omit(data_1RD8_AB))
g <- pcr( data_1RD8_AB$entropy ~ scale() + scale(log(cai)) + scale(log(abund)) + scale(log(len)) + scale(disp) + scale(log(deg)) + scale(log(BC+BCadd)), data=d, y=T)
bars <- mysum( g )
if (file_out) pdf(paste(dir_out, "dn-pcr.pdf", sep=''), family="Times", 
	width=10, height=7) else par(font=1, family="serif")
barplot(bars, col=seven_colors)
if (file_out) dev.off() else par(def.par)
print( summary( lm( g$y ~ g$scores ) ) )