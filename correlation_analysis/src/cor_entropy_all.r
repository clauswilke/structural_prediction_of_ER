# This function compares correlations of sequence entropy with RSA, WCN, Chi1 & RMSF from MD simulations. All measures of the variables are chosen from MD, since either the values from crystal structures are highly similar to MD results, or the statistic from crystal structures is not significant enough. It also creates the corresponding plot.

# Amir Shahmoradi, Friday 1:14 PM, Jan 31, 2014, Wilke Lab, ICMB, UT Austin

# source("input_data.r")

setwd('C:/Users/Amir/Documents/GitHub/structural_prediction_of_ER/')

data = rbind(data_1RD8_AB, data_2FP7_B, data_2JLY_A, data_2Z83_A, data_3GOL_A, data_3LYF_A, data_4AQF_B, data_4GHA_A, data_4IRY_A,
             data_3GSZ_A, data_3I5K_A, data_2JLY_A_temp_50, data_2JLY_A_temp_100, data_2JLY_A_temp_200, data_2JLY_A_temp_450)

data$protein = factor(data$protein)

result = data.frame()

for(protein in levels(data$protein))
{
    d = data[data$protein==protein,]
    
	x = cor.test( d$rsa_avg_md, d$entropy, method="spearman", na.action="na.omit" )
    r.rsa_avg_md = x$estimate
    p.rsa_avg_md = x$p.value
    
	x = cor.test( d$wcn_avg_md, d$entropy, method="spearman", na.action="na.omit" )
    r.wcn_avg_md = x$estimate
    p.wcn_avg_md = x$p.value

	x = cor.test( d$chi1_var_md, d$entropy, method="spearman", na.action="na.omit" )
    r.chi1_var_md = x$estimate
    p.chi1_var_md = x$p.value
	
	x = cor.test( d$rmsf_avg_md, d$entropy, method="spearman", na.action="na.omit" )
    r.rmsf_avg_md = x$estimate
    p.rmsf_avg_md = x$p.value
	
	x = cor.test( d$bfca, d$entropy, method="spearman", na.action="na.omit" )
    r.bfca = x$estimate
    p.bfca = x$p.value
	
	x = cor.test( d$desent, d$entropy, method="spearman", na.action="na.omit" )
    r.desent = x$estimate
    p.desent = x$p.value

#	if (length(na.omit(d$rmsfHS))>=3)
#	{
#		x = cor.test( d$rmsfHS, d$entropy, method="spearman", na.action="na.omit" )
#		r.rmsf_cr = x$estimate
#		p.rmsf_cr = x$p.value
#	}
#	else
#	{
#		r.rmsf_cr = NA
#		p.rmsf_cr = NA
#	}
    
    row = data.frame( protein=protein,
                      r.rsa_avg_md = r.rsa_avg_md, p.rsa_avg_md = p.rsa_avg_md, r.wcn_avg_md = -r.wcn_avg_md, p.wcn_avg_md = p.wcn_avg_md,
					  r.chi1_var_md = r.chi1_var_md, p.chi1_var_md = p.chi1_var_md, r.rmsf_avg_md = r.rmsf_avg_md, p.rmsf_avg_md = p.rmsf_avg_md,
                      r.bfca = r.bfca, p.bfca = p.bfca, r.desent = r.desent, p.desent = p.desent) # avg_cr_rho = r.rmsf_cr, avg_cr_P = p.rmsf_cr ) # , std_md_rho = r.rmsf_std_md, std_md_P = p.rmsf_std_md
    result = rbind( result, row )
    print( protein )
}
row.names(result) = c()
write.csv( result, "correlation_analysis/cor_tables/cor_entropy_all.csv", row.names=F )

index = names(result) %in% c("r.rsa_avg_md", "r.wcn_avg_md", "r.chi1_var_md", "r.rmsf_avg_md", "r.bfca", "r.desent") #, "avg_cr_rho") # columns we want to plot
index.p = names(result) %in% c("p.rsa_avg_md", "p.wcn_avg_md", "p.chi1_var_md", "p.rmsf_avg_md", "p.bfca", "p.desent") # columns that store significance

colors = c('red', 'blue', 'green', 'purple', 'orange3', 'darkgreen', 'black', 'gray', 'cyan2', 'darkred') #, 'darkgreen', 'bisque2')
proteins = c('1RD8_AB', '2FP7_B', '2JLY_A', '2Z83_A', '3GOL_A', '3LYF_A', '4AQF_B', '4GHA_A', '4IRY_A') #, '3GSZ_A', '3I5K_A', '2JLY_A_temp_50', '2JLY_A_temp_100', '2JLY_A_temp_200', '2JLY_A_temp_450')
labels = c('1RD8', '2FP7', '2JLY', '2Z83', '3GOL', '3LYF', '4AQF', '4GHA', '4IRY') #, '3GSZ_A', '3I5K_A', '2JLY_A_temp_50', '2JLY_A_temp_100', '2JLY_A_temp_200', '2JLY_A_temp_450')

pdf( "correlation_analysis/figures/cor_entropy_all.pdf", width=7, height=4, useDingbats=FALSE )
par( mai=c(0.65, 0.65, 0.0, 0.05), mgp=c(2, 0.5, 0), tck=-0.03 )
plot(0,xaxt='n',yaxt='n',bty='n',pch='',ylab=expression(paste("Correlation (", rho, ") with entropy")),xlab='Correlating Variable', xlim=c(1,6.2),ylim=c(-.5,.45)) #, main = 'sequence entropy - contact number association', cex.main=0.8)
#minor.tick(nx=0, ny=4, tick.ratio=2)
axis( 2,  # y axis
	  at=c(-.4,-.3,-.2,-.1,0,.1,.2,.3,.4),
	  c("-0.4","","-0.2","","0.0","","0.2","","0.4"))
axis( 1,
      at=c(1, 2, 3, 4, 5, 6),
      padj=c(0,0,0.2,0,0,1.5),
	  #hadj=0,
      c("MD RSA", "MD iWCN", expression(paste("MD Var(", chi[1],")")), "MD RMSF", "B factor", "designed\n entropy")) #,  "CS rmsf"))


for( i in 1:nrow(result) )
{
	if ( i <= length(proteins) )
	{
		if ( result$protein[i] == proteins[i] )
		{
			row = unlist( result[i,index] )
			x = 1:sum(index)
			print( i )
			lines( x, row, col=colors[i] )
			p = unlist( result[i,index.p] )
			sign = p < 0.05
			print(row)
			print(p)
			print(sign)
			points( x[sign], row[sign], pch=19, col=colors[i])
			points( x[!sign], row[!sign], pch=1, col=colors[i])
		}
		else
		{
			cat( "Error: data mismatch!\n" ) 
		}
	}
	else
	{
		cat( "skipping ", as.character(result$protein[i]), "\n" ) 
	}
}

legend( 1.0, -.28, labels[1:3], pch=19, col=colors[1:3], bty='n', cex=0.9)
legend( 1.7, -.28, labels[4:6], pch=19, col=colors[4:6], bty='n', cex=0.9)
legend( 2.4, -.28, labels[7:9], pch=19, col=colors[7:9], bty='n', cex=0.9)
dev.off()