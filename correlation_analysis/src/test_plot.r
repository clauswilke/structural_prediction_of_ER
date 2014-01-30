# This function calculates the correlation between the two variables entropy and RSA for each protein structure in the study and outputs the result into a CSV file.
# Amir Shahmoradi, Wednesday 3:06 PM, Jan 15, 2015, Wilke Lab, ICMB, UT Austin

# source("input_data.r")

setwd('C:/Users/Amir/Documents/GitHub/structural_prediction_of_ER/')

data = rbind(data_1RD8_AB, data_2FP7_B, data_2JLY_A, data_2Z83_A, data_3GOL_A, data_3GSZ_A, data_3I5K_A, data_3LYF_A, data_4AQF_B, data_4GHA_A, data_4IRY_A)
            #,data_2JLY_A_temp_50, data_2JLY_A_temp_100, data_2JLY_A_temp_200, data_2JLY_A_temp_450)

data$protein = factor(data$protein)

result = data.frame()

for(protein in levels(data$protein))
{
    d = data[data$protein==protein,]
    
	x = cor.test( d$cn13_cr, d$entropy, method="spearman", na.action="na.omit" )
    r.cn13_cr = x$estimate
    p.cn13_cr = x$p.value
    
	x = cor.test( d$cn13_avg_md, d$entropy, method="spearman", na.action="na.omit" )
    r.cn13_avg_md = x$estimate
    p.cn13_avg_md = x$p.value

	x = cor.test( d$cn13_var_md, d$entropy, method="spearman", na.action="na.omit" )
    r.cn13_var_md = x$estimate
    p.cn13_var_md = x$p.value

    row = data.frame( protein=protein, cr_rho = -r.cn13_cr, cr_P = p.cn13_cr, avg_md_rho = -r.cn13_avg_md, avg_md_P = p.cn13_avg_md, var_md_rho = -r.cn13_var_md, var_md_P = p.cn13_var_md )
    result = rbind( result, row )
    print( protein )
}
row.names(result) = c()
write.csv( result, "correlation_analysis/cor_entropy_cn13.csv", row.names=F )



colors = c('red', 'blue', 'green', 'purple', 'orange', 'yellow', 'black', 'gray', 'red', 'blue', 'green', 'purple')
proteins = c('1RD8_AB', '2FP7_B', '2JLY_A', '2Z83_A', '3GOL_A', '3GSZ_A', '3I5K_A', '3LYF_A', '4AQF_B', '4GHA_A', '4IRY_A')
index = match(result$protein,proteins)

pdf( "correlation_analysis/figures/test.pdf", width=4.5, height=4, useDingbats=FALSE )
par( mai=c(0.65, 0.65, 0.1, 0.05), mgp=c(2, 0.5, 0), tck=-0.03 )
plot(0,xaxt='n',yaxt='n',bty='n',pch='',ylab='Correlation',xlab='Variable', xlim=c(-1,1),ylim=c(.5,.5+length(index)))
points( result$cr_rho, index, pch=19, col=colors[i])
lines( result$cr_rho, index, col=colors[i] )
points( result$avg_md_rho, index, pch=19, col=colors[i])
lines( result$avg_md_rho, index, col=colors[i] )
points( result$var_md_rho, index, pch=19, col=colors[i])
lines( result$var_md_rho, index, col=colors[i] )

axis( 1 ) # y axis
axis( 2,
      at=index,
      proteins[index])

legend( 1.5, -.2, c("cr","avg_md","var_md"), pch=19, col=colors[1:3], bty='n' )
dev.off()