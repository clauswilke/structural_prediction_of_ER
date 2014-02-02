# This function calculates the correlation between the two variables entropy and contact number (CN13) for each protein structure in the study and outputs the result into a CSV file. It also generates the corresponding graph of entropy vs. CN13 correlations.

# Amir Shahmoradi, Friday 11:02 PM, Jan 24, 2014, Wilke Lab, ICMB, UT Austin

# source("input_data.r")

# install.packages("Hmisc", dependencies=TRUE) If Hmisc is not already installed, uncomment this command and let the code install it.
# library(Hmisc)	# adds minor ticks to plots

setwd('C:/Users/Amir/Documents/GitHub/structural_prediction_of_ER/')

data = rbind(data_1RD8_AB, data_2FP7_B, data_2JLY_A, data_2Z83_A, data_3GOL_A, data_3LYF_A, data_4AQF_B, data_4GHA_A, data_4IRY_A,
             data_3GSZ_A, data_3I5K_A, data_2JLY_A_temp_50, data_2JLY_A_temp_100, data_2JLY_A_temp_200, data_2JLY_A_temp_450)

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
write.csv( result, "correlation_analysis/cor_tables/cor_entropy_cn13.csv", row.names=F )


index = names(result) %in% c("cr_rho", "avg_md_rho", "var_md_rho") # columns we want to plot

colors = c('red', 'blue', 'green', 'purple', 'orange3', 'darkgreen', 'black', 'gray', 'cyan2') #, 'darkred', 'darkgreen', 'bisque2')
proteins = c('1RD8_AB', '2FP7_B', '2JLY_A', '2Z83_A', '3GOL_A', '3LYF_A', '4AQF_B', '4GHA_A', '4IRY_A') #, '3GSZ_A', '3I5K_A', '2JLY_A_temp_50', '2JLY_A_temp_100', '2JLY_A_temp_200', '2JLY_A_temp_450')

pdf( "correlation_analysis/figures/cor_entropy_cn13.pdf", width=4.5, height=4, useDingbats=FALSE )
par( mai=c(0.65, 0.65, 0.1, 0.05), mgp=c(2, 0.5, 0), tck=-0.03 )
plot(0,xaxt='n',yaxt='n',bty='n',pch='',ylab="Spearman's rank correlation coefficient",xlab='Contact Number: CN13', xlim=c(1,3.2),ylim=c(-.5,.5)) #, main = 'sequence entropy - contact number association', cex.main=0.8)
#minor.tick(nx=0, ny=4, tick.ratio=2)
axis( 2 ) # y axis
axis( 1,
      at=c(1, 2, 3),
      c("Crystal structure", "Avg MD", "Variance MD"))


for( i in 1:nrow(result) )
{
    row = unlist( result[i,index] )
    x = 1:sum(index)
    points( x, row, pch=19, col=colors[i])
    lines( x, row, col=colors[i] )
}

legend( 1.0, -.28, proteins[1:3], pch=19, col=colors[1:3], bty='n', cex=0.9)
legend( 1.7, -.28, proteins[4:6], pch=19, col=colors[4:6], bty='n', cex=0.9)
legend( 2.4, -.28, proteins[7:9], pch=19, col=colors[7:9], bty='n', cex=0.9)
dev.off()