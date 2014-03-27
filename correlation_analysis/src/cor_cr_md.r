# This function calculates the correlation between the corresponding variables calculated from MD and Crystal Structures (CS) for each protein structure in the study and outputs the result into a CSV file. It also generates the corresponding graphs of MD_variable vs. CS_variable correlations.
# The variables considered are: RSA_CS-RSA_MD, RMSF_CS-RMSF_MD, CN13_CS-CN13_MD, WCN_CS-WCN_MD

# Amir Shahmoradi, Friday 7:57 PM, Jan 31, 2014, Wilke Lab, ICMB, UT Austin

# source("input_data.r")

setwd('C:/Users/Amir/Documents/GitHub/structural_prediction_of_ER/')

data = rbind(data_1RD8_AB, data_2FP7_B, data_2JLY_A, data_2Z83_A, data_3GOL_A, data_3LYF_A, data_4AQF_B, data_4GHA_A, data_4IRY_A)
            #,data_3GSZ_A, data_3I5K_A, data_2JLY_A_temp_50, data_2JLY_A_temp_100, data_2JLY_A_temp_200, data_2JLY_A_temp_450)

data$protein = factor(data$protein)

result = data.frame()

for(protein in levels(data$protein))
{
    d = data[data$protein==protein,]
    
	# entropy - RSA:
	
		x = cor.test( d$rsa_cr, d$entropy, method="spearman", na.action="na.omit" )
		r.rsa_cr = x$estimate
		p.rsa_cr = x$p.value
		
		x = cor.test( d$rsa_avg_md, d$entropy, method="spearman", na.action="na.omit" )
		r.rsa_avg_md = x$estimate
		p.rsa_avg_md = x$p.value
    
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
	
	# entropy - CN13:
	
		x = cor.test( d$cn13_cr, d$entropy, method="spearman", na.action="na.omit" )
		r.cn13_cr = x$estimate
		p.cn13_cr = x$p.value
		
		x = cor.test( d$cn13_avg_md, d$entropy, method="spearman", na.action="na.omit" )
		r.cn13_avg_md = x$estimate
		p.cn13_avg_md = x$p.value
	
	# entropy - WCN:
	
		x = cor.test( d$wcn_cr, d$entropy, method="spearman", na.action="na.omit" )
		r.wcn_cr = x$estimate
		p.wcn_cr = x$p.value
		
		x = cor.test( d$wcn_avg_md, d$entropy, method="spearman", na.action="na.omit" )
		r.wcn_avg_md = x$estimate
		p.wcn_avg_md = x$p.value
	# entropy - B factor:
		
		x = cor.test( d$bfca, d$entropy, method="spearman", na.action="na.omit" )
		r.bfca = x$estimate
		p.bfca = x$p.value
	
	
    row = data.frame( protein=protein,
                      r.rsa_cr  =  r.rsa_cr,  p.rsa_cr  = p.rsa_cr,  r.rsa_avg_md  = r.rsa_avg_md,   p.rsa_avg_md = p.rsa_avg_md,
                      r.rmsf_cr =  r.rmsf_cr, p.rmsf_cr = p.rmsf_cr, r.rmsf_avg_md = r.rmsf_avg_md,  p.rmsf_avg_md = p.rmsf_avg_md,
                      r.cn13_cr = -r.cn13_cr, p.cn13_cr = p.cn13_cr, r.cn13_avg_md = -r.cn13_avg_md, p.cn13_avg_md = p.cn13_avg_md,
                      r.wcn_cr  = -r.wcn_cr,  p.wcn_cr  = p.wcn_cr,  r.wcn_avg_md  = -r.wcn_avg_md,  p.wcn_avg_md = p.wcn_avg_md,
                      r.bfca = r.bfca, p.bfca = p.bfca)
    result = rbind( result, row )
    print( protein )
}
row.names(result) = c()
#write.csv( result, "correlation_analysis/cor_tables/cor_entropy_cr_md.csv", row.names=F )


proteins = c('1RD8_AB', '2FP7_B', '2JLY_A', '2Z83_A', '3GOL_A', '3LYF_A', '4AQF_B', '4GHA_A', '4IRY_A') #, '3GSZ_A', '3I5K_A', '2JLY_A_temp_50', '2JLY_A_temp_100', '2JLY_A_temp_200', '2JLY_A_temp_450')

pdf( "correlation_analysis/figures/cor_cr_md.pdf", width=9.3, height=8.8, useDingbats=FALSE )

split.screen(c(2,2))

screen(3)
	par( mai=c(0.65, 0.65, 0.4, 0.2), mgp=c(2, 0.5, 0), tck=-0.03 )
	plot(0,xaxt='n',yaxt='n',bty='n',pch='',xlab="entropy - RSA (Crystal Structure)",ylab="entropy - RSA (MD simulation)", xlim=c(-.1,.4),ylim=c(-.1,.4))
	axis( 1,  # x axis
		at=c(-.1,0,.1,.2,.3,.4),
		c("","0.0","","0.2","","0.4"))
	axis( 2,  # y axis
		at=c(-.1,0,.1,.2,.3,.4),
		c("","0.0","","0.2","","0.4"))
	points( result$r.rsa_cr, result$r.rsa_avg_md, pch=19 )
	abline(0,1)
  mtext('C', side = 3, at=-.2, font=2, cex=1.2)

screen(4)
	par( mai=c(0.65, 0.65, 0.4, 0.2), mgp=c(2, 0.5, 0), tck=-0.03 )
	plot(0,xaxt='n',yaxt='n',bty='n',pch='',xlab="entropy - RMSF (Crystal Structure)",ylab="entropy - RMSF (MD simulation)", xlim=c(-.1,.4),ylim=c(-.1,.4))
	axis( 1,  # x axis
		at=c(-.1,0,.1,.2,.3,.4),
		c("","0.0","","0.2","","0.4"))
	axis( 2,  # y axis
		at=c(-.1,0,.1,.2,.3,.4),
		c("","0.0","","0.2","","0.4"))
	points( result$r.rmsf_cr, result$r.rmsf_avg_md, pch=19 )
	abline(0,1)
	mtext('D', side = 3, at=-.2, font=2, cex=1.2)
	

screen(1)
	par( mai=c(0.65, 0.65, 0.4, 0.2), mgp=c(2, 0.5, 0), tck=-0.03 )
	plot(0,xaxt='n',yaxt='n',bty='n',pch='',xlab="entropy - iCN (Crystal Structure)",ylab="entropy - iCN (MD simulation)", xlim=c(-.1,.4),ylim=c(-.1,.4))
	axis( 1,  # x axis
		at=c(-.1,0,.1,.2,.3,.4),
		c("","0.0","","0.2","","0.4"))
	axis( 2,  # y axis
		at=c(-.1,0,.1,.2,.3,.4),
		c("","0.0","","0.2","","0.4"))
	points( result$r.cn13_cr, result$r.cn13_avg_md, pch=19 )
	abline(0,1)
	mtext('A', side = 3, at=-.2, font=2, cex=1.2)
	
screen(2)
	par( mai=c(0.65, 0.65, 0.4, 0.2), mgp=c(2, 0.5, 0), tck=-0.03 )
	plot(0,xaxt='n',yaxt='n',bty='n',pch='',xlab="entropy - iWCN (Crystal Structure)",ylab="entropy - iWCN (MD simulation)", xlim=c(-.1,.4),ylim=c(-.1,.4))
	axis( 1,  # x axis
		at=c(-.1,0,.1,.2,.3,.4),
		c("","0.0","","0.2","","0.4"))
	axis( 2,  # y axis
		at=c(-.1,0,.1,.2,.3,.4),
		c("","0.0","","0.2","","0.4"))	
	points( result$r.wcn_cr, result$r.wcn_avg_md, pch=19 )
	abline(0,1)
	mtext('B', side = 3, at=-.2, font=2, cex=1.2)
	
#screen(5)
#	par( mai=c(0.65, 0.65, 0.1, 0.05), mgp=c(2, 0.5, 0), tck=-0.03 )
#	plot(0,xaxt='n',yaxt='n',bty='n',pch='',xlab="entropy - B factor (Crystal Structure)",ylab="entropy - RMSF (MD simulation)", xlim=c(-.1,.4),ylim=c(-.1,.4))
#	axis( 1,  # x axis
#		at=c(-.4,-.3,-.2,-.1,0,.1,.2,.3,.4),
#		c("-0.4","","-0.2","","0.0","","0.2","","0.4"))
#	axis( 2,  # y axis
#		at=c(-.4,-.3,-.2,-.1,0,.1,.2,.3,.4),
#		c("-0.4","","-0.2","","0.0","","0.2","","0.4"))
#	
#		
#	points( result$r.bfca, result$r.rmsf_avg_md, pch=19 )
#	abline(0,1)
#
#screen(6)
#	par( mai=c(0.65, 0.65, 0.1, 0.05), mgp=c(2, 0.5, 0), tck=-0.03 )
#	plot(0,xaxt='n',yaxt='n',bty='n',pch='',xlab="entropy - B factor (Crystal Structure)",ylab="entropy - RMSF (Homologous Structures)", xlim=c(-.4,.4),ylim=c(-.4,.4))
#	axis( 1,  # x axis
#		at=c(-.4,-.3,-.2,-.1,0,.1,.2,.3,.4),
#		c("-0.4","","-0.2","","0.0","","0.2","","0.4"))
#	axis( 2,  # y axis
#		at=c(-.4,-.3,-.2,-.1,0,.1,.2,.3,.4),
#		c("-0.4","","-0.2","","0.0","","0.2","","0.4"))
#	
#		
#	points( result$r.bfca, result$r.rmsf_cr, pch=19 )
#	abline(0,1)


close.screen(all = TRUE)

dev.off()