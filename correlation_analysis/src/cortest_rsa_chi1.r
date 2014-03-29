#Calculate the correlation between rmsf and dihedral angle chi1:
print(cor.test(data_1RD8_AB$rsa_cr,data_1RD8_AB$rmsf_avg_md, method='spearman'))
print( cor.test(data_2JLY_A$rsa_cr, data_2JLY_A$rmsf_avg_md, method='spearman'))
print( cor.test(data_2Z83_A$rsa_cr, data_2Z83_A$rmsf_avg_md, method='spearman'))
print( cor.test(data_3GOL_A$rsa_cr, data_3GOL_A$rmsf_avg_md, method='spearman'))
print( cor.test(data_3GSZ_A$rsa_cr, data_3GSZ_A$rmsf_avg_md, method='spearman'))
print( cor.test(data_3I5K_A$rsa_cr, data_3I5K_A$rmsf_avg_md, method='spearman'))
print( cor.test(data_3LYF_A$rsa_cr, data_3LYF_A$rmsf_avg_md, method='spearman'))
print( cor.test(data_4AQF_B$rsa_cr, data_4AQF_B$rmsf_avg_md, method='spearman'))
print( cor.test(data_4GHA_A$rsa_cr, data_4GHA_A$rmsf_avg_md, method='spearman'))
print( cor.test(data_4IRY_A$rsa_cr, data_4IRY_A$chi1_var_md, method='spearman'))
#print(cor.test(data_1AJ8_A$AvgRMSD, dihedral_1AJ8_A$var_chi, method='spearman'))
#print(cor.test(data_1AOR_A$AvgRMSD, dihedral_1AOR_A$var_chi, method='spearman'))
#print(cor.test(data_1CTS_A$AvgRMSD, dihedral_1CTS_A$var_chi, method='spearman'))
#print(cor.test(data_1MP9_A$AvgRMSD, dihedral_1MP9_A$var_chi, method='spearman'))