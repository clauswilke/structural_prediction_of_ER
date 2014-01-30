# This function calculates the correlation between the two variables entropy and rmsf for each protein structure in the study and outputs the result into a CSV file. It also generates the corresponding graphs of entropy vs. rmsf correlations.

# Amir Shahmoradi, Saturday 1:51 AM, Jan 25, 2014, Wilke Lab, ICMB, UT Austin

# source("input_data.r")

setwd('C:/Users/Amir/Documents/GitHub/structural_prediction_of_ER/')

data = rbind(data_1RD8_AB, data_2FP7_B, data_2JLY_A, data_2Z83_A, data_3GOL_A, data_3LYF_A, data_4AQF_B, data_4GHA_A, data_4IRY_A,
             data_3GSZ_A, data_3I5K_A, data_2JLY_A_temp_50, data_2JLY_A_temp_100, data_2JLY_A_temp_200, data_2JLY_A_temp_450)

data$protein = factor(data$protein)

result = data.frame()

for(protein in levels(data$protein))
{
    d = data[data$protein==protein,]
    
	x = cor.test( d$rmsf_avg_md, d$entropy, method="spearman", na.action="na.omit" )
    r.rmsf_cr = x$estimate
    p.rmsf_cr = x$p.value
    
	x = cor.test( d$rmsf_avg_md, d$entropy, method="spearman", na.action="na.omit" )
    r.rmsf_avg_md = x$estimate
    p.rmsf_avg_md = x$p.value

	x = cor.test( d$rmsf_std_md, d$entropy, method="spearman", na.action="na.omit" )
    r.rmsf_std_md = x$estimate
    p.rmsf_std_md = x$p.value

    row = data.frame( protein=protein, cr_rho = r.rmsf_avg_md, cr_P = p.rmsf_cr, avg_md_rho = r.rmsf_avg_md, avg_md_P = p.rmsf_avg_md, std_md_rho = r.rmsf_std_md, std_md_P = p.rmsf_std_md )
    result = rbind( result, row )
    print( protein )
}
row.names(result) = c()
write.csv( result, "correlation_analysis/cor_tables/cor_entropy_cr_md.csv", row.names=F )


index = names(result) %in% c("cr_rho", "avg_md_rho", "std_md_rho") # columns we want to plot
index.p = names(result) %in% c("cr_P", "avg_md_P", "std_md_P") # columns that store significance

proteins = c('1RD8_AB', '2FP7_B', '2JLY_A', '2Z83_A', '3GOL_A', '3LYF_A', '4AQF_B', '4GHA_A', '4IRY_A') #, '3GSZ_A', '3I5K_A', '2JLY_A_temp_50', '2JLY_A_temp_100', '2JLY_A_temp_200', '2JLY_A_temp_450')

pdf( "correlation_analysis/figures/cor_entropy_cr_md.pdf", width=9, height=8, useDingbats=FALSE )

split.screen(c(2,2))
screen(1)
par( mai=c(0.65, 0.65, 0.1, 0.05), mgp=c(2, 0.5, 0), tck=-0.03 )
plot(0,xaxt='n',yaxt='n',bty='n',pch='',xlab="entropy vs RSA (crystal)",ylab="entropy vs RSA (MD)", xlim=c(-.4,.4),ylim=c(-.4,.4))
axis( 1,  # x axis
	  at=c(-.4,-.3,-.2,-.1,0,.1,.2,.3,.4),
	  c("-0.4","","-0.2","","0.0","","0.2","","0.4"))
axis( 2,  # y axis
	  at=c(-.4,-.3,-.2,-.1,0,.1,.2,.3,.4),
	  c("-0.4","","-0.2","","0.0","","0.2","","0.4"))

	  
points( result$cr_rho, result$avg_md_rho, pch=19 )
abline(0,1)

screen(3)
par( mai=c(0.65, 0.65, 0.1, 0.05), mgp=c(2, 0.5, 0), tck=-0.03 )
plot(0,xaxt='n',yaxt='n',bty='n',pch='',xlab="entropy vs RMSF (crystal)",ylab="entropy vs RMSF (MD)", xlim=c(-.4,.4),ylim=c(-.4,.4))
axis( 1,  # x axis
	  at=c(-.4,-.3,-.2,-.1,0,.1,.2,.3,.4),
	  c("-0.4","","-0.2","","0.0","","0.2","","0.4"))
axis( 2,  # y axis
	  at=c(-.4,-.3,-.2,-.1,0,.1,.2,.3,.4),
	  c("-0.4","","-0.2","","0.0","","0.2","","0.4"))

	  
points( result$cr_rho, result$avg_md_rho, pch=19 )
abline(0,1)



close.screen(all = TRUE)

dev.off()