# This function calculates the correlation matrix of the variables considered for each protein structure in the study.
# Amir Shahmoradi, Wednesday 1:30 PM, Jan 15, 2015, Wilke Lab, ICMB, UT Austin

data = rbind(data_1RD8, data_2JLY, data_2Z83, data_3GOL, data_3GSZ, data_3I5K, data_4AQF_B, data_4GHA_A, data_4IRY,
             data_2JLY_temp_50, data_2JLY_temp_100, data_2JLY_temp_200, data_2JLY_temp_450)
data$protein = factor(data$protein)

result_rsa = data.frame()
result_rmsf = data.frame()
result_dihedral = data.frame()
result_cn13 = data.frame()
result_wcn = data.frame()

for(protein in levels(data$protein))
{
    d = data[data$protein==protein,]
    
	#x = cor.test( d$phi_var_md[-1], d$entropy[-1] )
	x = cor.test( d$phi_var_md, d$entropy, na.action="na.omit" )
    r.phi_var_md = x$estimate
    p.phi_var_md = x$p.value
    
	
	#x = cor.test( d$psi_var_md[-nrow(d)], d$entropy[-nrow(d)] )
    r.psi_var_md = x$estimate
    p.psi_var_md = x$p.value
    row = data.frame( protein=protein, r.phi_var_md, p.phi_var_md, r.psi_var_md, p.psi_var_md )
    result = rbind( result, row )
    print( protein )
}
row.names(result) = c()
write.csv( result, "correlations.csv", row.names=F )

##########################################################

#analyze.cor = function( data, columns, ... )
#{
#    matrix = data[,colnames(data) %in% columns]
#    cor(na.omit(matrix), ...)
#}
#
#
#
#for(protein in levels(data$protein))
#{
#    d = data[data$protein==protein,]
#    columns = c( "phi_var_md", "psi_var_md", "entropy" )
#    m.pears = analyze.cor( d, columns )
#    m.spear = analyze.cor( d, columns, method="sp")
#    print( protein )
#    print( m.pears )
#    print( m.spear )
#    fname = paste( "cor_matrix_", protein, ".csv", sep='' )
#    write.csv( m.pears, fname )
#    
#}

