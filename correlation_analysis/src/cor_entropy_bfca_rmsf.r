# This function calculates the correlation between entropy and RMSF_CR, RMSF_MD and CA B-factors (bfca), then compares the correlations of rmsf-entropy with bfca-entropy, and outputs the data ina CSV file and the plots in a pdf file.

# Amir Shahmoradi, Monday 5:44 PM, Feb 17, 2014, Wilke Lab, ICMB, UT Austin

# source("input_data.r")

setwd('C:/Users/Amir/Documents/GitHub/structural_prediction_of_ER/')

data = rbind(data_1RD8_AB, data_2FP7_B, data_2JLY_A, data_2Z83_A, data_3GOL_A, data_3LYF_A, data_4AQF_B, data_4GHA_A, data_4IRY_A)
            #,data_3GSZ_A, data_3I5K_A, data_2JLY_A_temp_50, data_2JLY_A_temp_100, data_2JLY_A_temp_200, data_2JLY_A_temp_450)

data$protein = factor(data$protein)

result = data.frame()

for(protein in levels(data$protein))
{
    d = data[data$protein==protein,]
    
	# entropy - RMSF:
	
		x = cor.test( d$rmsf_avg_md, d$entropy, method="spearman", na.action="na.omit" )
		r.rmsf_avg_md = x$estimate
		p.rmsf_avg_md = x$p.value
	
		if (length(na.omit(d$rmsfHS))>=3)
		{
			x = cor.test( d$rmsfHS, d$entropy, method="spearman", na.action="na.omit" )
			r.rmsf_cr = x$estimate
			p.rmsf_cr = x$p.value
		}
		else
		{
			r.rmsf_cr = NA
			p.rmsf_cr = NA
		}	

	# entropy - B factor:
		
		x = cor.test( d$bfca, d$entropy, method="spearman", na.action="na.omit" )
		r.bfca = x$estimate
		p.bfca = x$p.value
	
	
    row = data.frame( protein=protein,
                      r.rmsf_cr =  r.rmsf_cr, p.rmsf_cr = p.rmsf_cr, r.rmsf_avg_md = r.rmsf_avg_md,  p.rmsf_avg_md = p.rmsf_avg_md,
                      r.bfca = r.bfca, p.bfca = p.bfca)
    result = rbind( result, row )
    print( protein )
}
row.names(result) = c()
#write.csv( result, "correlation_analysis/cor_tables/cor_entropy_bfca_rmsf.csv", row.names=F )


proteins = c('1RD8_AB', '2FP7_B', '2JLY_A', '2Z83_A', '3GOL_A', '3LYF_A', '4AQF_B', '4GHA_A', '4IRY_A') #, '3GSZ_A', '3I5K_A', '2JLY_A_temp_50', '2JLY_A_temp_100', '2JLY_A_temp_200', '2JLY_A_temp_450')

pdf( "correlation_analysis/figures/cor_entropy_bfca_rmsf.pdf", width=9.3, height=4.4, useDingbats=FALSE )

split.screen(c(1,2))

screen(1)
	par( mai=c(0.65, 0.65, 0.4, 0.2), mgp=c(2, 0.5, 0), tck=-0.03 )
	plot(0,xaxt='n',yaxt='n',bty='n',pch='',xlab="entropy - B factor",ylab="entropy - MD RMSF", xlim=c(-.1,.4),ylim=c(-.1,.4))
	axis( 1,  # x axis
		at=c(-.1,0,.1,.2,.3,.4),
		c("","0.0","","0.2","","0.4"))
	axis( 2,  # y axis
		at=c(-.1,0,.1,.2,.3,.4),
		c("","0.0","","0.2","","0.4"))
	points( result$r.bfca, result$r.rmsf_avg_md, pch=19 )
	abline(0,1)
	mtext('A', side = 3, at=-.2, font=2, cex=1.2)

screen(2)
	par( mai=c(0.65, 0.65, 0.4, 0.2), mgp=c(2, 0.5, 0), tck=-0.03 )
	plot(0,xaxt='n',yaxt='n',bty='n',pch='',xlab="entropy - B factor",ylab="entropy - CS RMSF", xlim=c(-.1,.4),ylim=c(-.1,.4))
	axis( 1,  # x axis
		at=c(-.1,0,.1,.2,.3,.4),
		c("","0.0","","0.2","","0.4"))
	axis( 2,  # y axis
		at=c(-.1,0,.1,.2,.3,.4),
		c("","0.0","","0.2","","0.4"))
	points( result$r.bfca, result$r.rmsf_cr, pch=19 )
	abline(0,1)
	mtext('B', side = 3, at=-.2, font=2, cex=1.2)

close.screen(all = TRUE)

dev.off()
