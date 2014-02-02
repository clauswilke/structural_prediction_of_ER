# This function compares the bfactor with rmsf calculated from MD and Homologous crystal Structures (HS) for each protein structure in the study and outputs the result into a CSV file. It also generates the corresponding graph of bfactor (bfca) vs. rmsf (MD) and bfactor (bfca) vs. rmsf (HS).


# Amir Shahmoradi, Saturday 4:07 PM, Feb 1 2014, Wilke Lab, ICMB, UT Austin

# source("input_data.r")

setwd('C:/Users/Amir/Documents/GitHub/structural_prediction_of_ER/')

data = rbind(data_1RD8_AB, data_2FP7_B, data_2JLY_A, data_2Z83_A, data_3GOL_A, data_3LYF_A, data_4AQF_B, data_4GHA_A, data_4IRY_A,
             data_3GSZ_A, data_3I5K_A, data_2JLY_A_temp_50, data_2JLY_A_temp_100, data_2JLY_A_temp_200, data_2JLY_A_temp_450)

data$protein = factor(data$protein)

result = data.frame()

for(protein in levels(data$protein))
{
    d = data[data$protein==protein,]
    
	# bfactor - rmsf (MD):
	
		x = cor.test( d$bfca, d$rmsf_avg_md, method="spearman", na.action="na.omit" )
		r.bfca_rmsfMD = x$estimate
		p.bfca_rmsfMD = x$p.value
		
	# bfactor - rmsf (HS):
	
		if (length(na.omit(d$rmsfHS))>=3)
		{
			x = cor.test( d$bfca, d$rmsfHS, method="spearman", na.action="na.omit" )
			r.bfca_rmsfHS = x$estimate
			p.bfca_rmsfHS = x$p.value
		}
		else
		{
			r.bfca_rmsfHS = NA
			p.bfca_rmsfHS = NA
		}	
	
    row = data.frame( protein=protein,r.bfca_rmsfHS = r.bfca_rmsfHS, p.bfca_rmsfHS = p.bfca_rmsfHS, r.bfca_rmsfMD = r.bfca_rmsfMD, p.bfca_rmsfMD = p.bfca_rmsfMD)
    result = rbind( result, row )
    print( protein )
}
row.names(result) = c()
write.csv( result, "correlation_analysis/cor_tables/cor_bfca_rmsf.csv", row.names=F )


proteins = c('1RD8_AB', '2FP7_B', '2JLY_A', '2Z83_A', '3GOL_A', '3LYF_A', '4AQF_B', '4GHA_A', '4IRY_A') #, '3GSZ_A', '3I5K_A', '2JLY_A_temp_50', '2JLY_A_temp_100', '2JLY_A_temp_200', '2JLY_A_temp_450')

pdf( "correlation_analysis/figures/cor_bfca_rmsf.pdf", width=4.5, height=4, useDingbats=FALSE )

split.screen(c(1,1))

screen(1)
	par( mai=c(0.65, 0.65, 0.1, 0.05), mgp=c(2, 0.5, 0), tck=-0.03 )
	plot(0,xaxt='n',yaxt='n',bty='n',pch='',xlab="B factor - RMSF (Homologous Structures)",ylab="B factor - RMSF (MD simulation)", xlim=c(0,.8),ylim=c(0,.8))
	axis( 1,  # x axis
		at=c(0.0,0.1,.2,0.3,0.4,0.5,0.6,0.7,0.8),
		c("0.0","","0.2","","0.4","","0.6","","0.8"))
	axis( 2,  # y axis
		at=c(0.0,0.1,.2,0.3,0.4,0.5,0.6,0.7,0.8),
		c("0.0","","0.2","","0.4","","0.6","","0.8"))
	
		
	points( result$r.bfca_rmsfHS, result$r.bfca_rmsfMD, pch=19 )
	abline(0,1)

#screen(2)
#	par( mai=c(0.65, 0.65, 0.1, 0.05), mgp=c(2, 0.5, 0), tck=-0.03 )
#	plot(0,xaxt='n',yaxt='n',bty='n',pch='',xlab="B factor",ylab="RMSF (Homologous Structures)", xlim=c(-.4,.4),ylim=c(-.4,.4))
#	axis( 1,  # x axis
#		at=c(-.4,-.3,-.2,-.1,0,.1,.2,.3,.4),
#		c("-0.4","","-0.2","","0.0","","0.2","","0.4"))
#	axis( 2,  # y axis
#		at=c(-.4,-.3,-.2,-.1,0,.1,.2,.3,.4),
#		c("-0.4","","-0.2","","0.0","","0.2","","0.4"))
#	
#		
#	points( result$r.rmsf_cr, result$r.rmsf_avg_md, pch=19 )
#	abline(0,1)

close.screen(all = TRUE)

dev.off()