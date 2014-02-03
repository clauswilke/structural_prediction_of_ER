# This function calculates the correlation between the two variables entropy and omega (site-specific evolutionary rates) for each protein structure in the study and outputs the result into a CSV file. It also generates the corresponding graphs of correlations between entropy and the evolutionary rates (omega).

# Amir Shahmoradi, Saturday 11:14 PM, Feb  1, 2014, ICMB, UT Austin

# source("input_data.r")

# install.packages("Hmisc", dependencies=TRUE) If Hmisc is not already installed, uncomment this command and let the code install it.
#library(Hmisc)	# adds minor ticks to plots

setwd('C:/Users/Amir/Documents/GitHub/structural_prediction_of_ER/')

data = rbind(data_1RD8_AB, data_2FP7_B, data_2JLY_A, data_2Z83_A, data_3GOL_A, data_3LYF_A, data_4AQF_B, data_4GHA_A, data_4IRY_A,
             data_3GSZ_A, data_3I5K_A, data_2JLY_A_temp_50, data_2JLY_A_temp_100, data_2JLY_A_temp_200, data_2JLY_A_temp_450)

data$protein = factor(data$protein)

result = data.frame()

for(protein in levels(data$protein))
{
    d = data[data$protein==protein,]
    
	# ENTROPY - VARIABLES :
	
		x = cor.test( d$rsa_avg_md, d$entropy, method="spearman", na.action="na.omit" )
		r.entropy_rsa_avg_md = x$estimate
		p.entropy_rsa_avg_md = x$p.value
		
		x = cor.test( d$wcn_avg_md, d$entropy, method="spearman", na.action="na.omit" )
		r.entropy_wcn_avg_md = x$estimate
		p.entropy_wcn_avg_md = x$p.value
	
		x = cor.test( d$chi1_var_md, d$entropy, method="spearman", na.action="na.omit" )
		r.entropy_chi1_var_md = x$estimate
		p.entropy_chi1_var_md = x$p.value
		
		x = cor.test( d$rmsf_avg_md, d$entropy, method="spearman", na.action="na.omit" )
		r.entropy_rmsf_avg_md = x$estimate
		p.entropy_rmsf_avg_md = x$p.value
		
		x = cor.test( d$bfca, d$entropy, method="spearman", na.action="na.omit" )
		r.entropy_bfca = x$estimate
		p.entropy_bfca = x$p.value
		
		x = cor.test( d$cn13_avg_md, d$entropy, method="spearman", na.action="na.omit" )
		r.entropy_cn13_avg_md = x$estimate
		p.entropy_cn13_avg_md = x$p.value
		
	# OMEGA - VARIABLES :
	
		x = cor.test( d$rsa_avg_md, d$omega, method="spearman", na.action="na.omit" )
		r.omega_rsa_avg_md = x$estimate
		p.omega_rsa_avg_md = x$p.value
		
		x = cor.test( d$wcn_avg_md, d$omega, method="spearman", na.action="na.omit" )
		r.omega_wcn_avg_md = x$estimate
		p.omega_wcn_avg_md = x$p.value
	
		x = cor.test( d$chi1_var_md, d$omega, method="spearman", na.action="na.omit" )
		r.omega_chi1_var_md = x$estimate
		p.omega_chi1_var_md = x$p.value
		
		x = cor.test( d$rmsf_avg_md, d$omega, method="spearman", na.action="na.omit" )
		r.omega_rmsf_avg_md = x$estimate
		p.omega_rmsf_avg_md = x$p.value
		
		x = cor.test( d$bfca, d$omega, method="spearman", na.action="na.omit" )
		r.omega_bfca = x$estimate
		p.omega_bfca = x$p.value
		
		x = cor.test( d$cn13_avg_md, d$omega, method="spearman", na.action="na.omit" )
		r.omega_cn13_avg_md = x$estimate
		p.omega_cn13_avg_md = x$p.value
		
    row = data.frame( protein=protein,
                      r.entropy_rsa_avg_md  =  r.entropy_rsa_avg_md,  p.entropy_rsa_avg_md  = p.entropy_rsa_avg_md,
                      r.entropy_wcn_avg_md  = -r.entropy_wcn_avg_md,  p.entropy_wcn_avg_md  = p.entropy_wcn_avg_md,
                      r.entropy_cn13_avg_md =  r.entropy_cn13_avg_md, p.entropy_cn13_avg_md = p.entropy_cn13_avg_md,
                      r.entropy_chi1_var_md =  r.entropy_chi1_var_md, p.entropy_chi1_var_md = p.entropy_chi1_var_md,
                      r.entropy_rmsf_avg_md =  r.entropy_rmsf_avg_md, p.entropy_rmsf_avg_md = p.entropy_rmsf_avg_md,
                      r.entropy_bfca        =  r.entropy_bfca,        p.entropy_bfca        = p.entropy_bfca,
					  r.omega_rsa_avg_md    =  r.omega_rsa_avg_md,    p.omega_rsa_avg_md    = p.omega_rsa_avg_md,
                      r.omega_wcn_avg_md    = -r.omega_wcn_avg_md,    p.omega_wcn_avg_md    = p.omega_wcn_avg_md,
                      r.omega_cn13_avg_md   =  r.omega_cn13_avg_md,   p.omega_cn13_avg_md   = p.omega_cn13_avg_md,
                      r.omega_chi1_var_md   =  r.omega_chi1_var_md,   p.omega_chi1_var_md   = p.omega_chi1_var_md,
                      r.omega_rmsf_avg_md   =  r.omega_rmsf_avg_md,   p.omega_rmsf_avg_md   = p.omega_rmsf_avg_md,
                      r.omega_bfca          =  r.omega_bfca,          p.omega_bfca          = p.omega_bfca )
	
    result = rbind( result, row )
    print( protein )
}
row.names(result) = c()
write.csv( result, "correlation_analysis/cor_tables/cor_entropy_omega.csv", row.names=F )


#index   = names(result) %in% c("r.entropy_rsa_avg_md", "r.entropy_wcn_avg_md", "r.entropy_chi1_var_md", "r.entropy_rmsf_avg_md", "r.entropy_bfca",
#                               "r.omega_rsa_avg_md", "r.omega_wcn_avg_md", "r.omega_chi1_var_md", "r.omega_rmsf_avg_md", "r.omega_bfca") # columns we want to plot
#index.p = names(result) %in% c("p.entropy_rsa_avg_md", "p.entropy_wcn_avg_md", "p.entropy_chi1_var_md", "p.entropy_rmsf_avg_md", "p.entropy_bfca",
#                               "p.omega_rsa_avg_md", "p.omega_wcn_avg_md", "p.omega_chi1_var_md", "p.omega_rmsf_avg_md", "p.omega_bfca") # columns that store significance

colors = c('red', 'blue', 'green', 'bisque3', 'black', 'gray') #, 'gray', 'cyan2') #, 'darkred', 'darkgreen', 'bisque2')
variables = c('MD RSA', 'MD WCN', expression(paste("MD var(", chi[1], ")")), 'MD RMSF', 'B factor')
labels = c('1RD8', '2FP7', '2JLY', '2Z83', '3GOL', '3LYF', '4AQF', '4GHA', '4IRY')

pdf( "correlation_analysis/figures/cor_entropy_omega.pdf", width=4.5, height=4, useDingbats=FALSE )
par( mai=c(0.65, 0.65, 0.1, 0.05), mgp=c(2, 0.5, 0), tck=-0.03 )
plot(0,xaxt='n',yaxt='n',bty='n',pch='',ylab=expression(paste("Correlation (", rho, ") with omega")),xlab=expression(paste("Correlation (", rho, ") with entropy")), xlim=c(-.5,.5),ylim=c(-.5,.5)) #, main = 'sequence entropy - contact number association', cex.main=0.8)
#minor.tick(nx=0, ny=4, tick.ratio=2)
#axis( 2 ) # y axis
#axis( 1 ) # , at=c(1, 2, 3, 4, 5), c("Designed Entropy" )) #, "Avg MD", "Variance MD"))


	axis( 1, at=c(-.4,-.3,-.2,-.1,0,.1,.2,.3,.4), c("-0.4","","-0.2","","0.0","","0.2","","0.4"))   # x axis
	axis( 2, at=c(-.4,-.3,-.2,-.1,0,.1,.2,.3,.4), c("-0.4","","-0.2","","0.0","","0.2","","0.4"))   # y axis
		
	points( result$r.entropy_rsa_avg_md , result$r.omega_rsa_avg_md , pch=19, col = colors[1] )
	points( result$r.entropy_wcn_avg_md , result$r.omega_wcn_avg_md , pch=19, col = colors[2] )
	points( result$r.entropy_chi1_var_md, result$r.omega_chi1_var_md, pch=19, col = colors[3] )
	points( result$r.entropy_rmsf_avg_md, result$r.omega_rmsf_avg_md, pch=19, col = colors[4] )
	points( result$r.entropy_bfca       , result$r.omega_bfca       , pch=19, col = colors[5] )
	#points( result$r.entropy_cn13_avg_md, result$r.omega_cn13_avg_md, pch=19, col = colors[6] )
	abline(0,1)
	#abline(-0.2,1)

legend( -.5, 0.5, variables[1:5], pch=19, col=colors[1:5], bty='n', cex=0.9)
#legend( 0.4, -.28, variables[4:5], pch=19, col=colors[4:5], bty='n', cex=0.9)
dev.off()