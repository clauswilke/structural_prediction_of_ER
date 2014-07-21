# This R code performs Principal Component ANALYSIS on different variables of THE COMBINED PROTEINS' DATA potentially correlating with OMEGA, then finds the rotation matrix, then outputs the component loadings and the percentage of the omega's explained variance by each principal component in a CVS file.

# Amir Shahmoradi, Wednesday 7:51 PM, March 26 2014, Wilke Lab, ICMB, UT Austin

# source("input_data.r")

out_path = 'C:/Users/Amir/Documents/GitHub/structural_prediction_of_ER/correlation_analysis/pcr_combined/pcr_combined_cors_omega.csv'

# Replace NAs with mean values in data:

data = rbind(data_1RD8_AB, data_2FP7_B, data_2JLY_A, data_2Z83_A, data_3GOL_A, data_3LYF_A, data_4AQF_B, data_4GHA_A, data_4IRY_A) #, data_3GSZ_A, data_3I5K_A) #, data_2JLY_A_temp_50, data_2JLY_A_temp_100, data_2JLY_A_temp_200, data_2JLY_A_temp_450)

pca_data = data[,names(data) %in% c('protein','omega','desent','rsa_avg_md','rmsf_avg_md','chi1_var_md','wcn_avg_md','bfca')] 
pca_data$chi1_var_md[is.na(pca_data$chi1_var_md)] = mean(na.omit(pca_data$chi1_var_md))

	pca_data$rsa_avg_md  = scale(pca_data$rsa_avg_md  , center = TRUE, scale = TRUE)
	pca_data$wcn_avg_md  = scale(1/pca_data$wcn_avg_md, center = TRUE, scale = TRUE)
	pca_data$chi1_var_md = scale(pca_data$chi1_var_md , center = TRUE, scale = TRUE)
	pca_data$rmsf_avg_md = scale(pca_data$rmsf_avg_md , center = TRUE, scale = TRUE)
	pca_data$bfca        = scale(pca_data$bfca        , center = TRUE, scale = TRUE)
	pca_data$desent      = scale(pca_data$desent      , center = TRUE, scale = TRUE)

	# Rename the column headers
	#	names(pca_data) <- c('protein','omega','designed_entropy','RSA','RMSF','VAR(chi1)','iWCN','Bfactor')
	
pca=prcomp(pca_data[,c(-1,-2)], scale = T, center=F)
#biplot(pca, xlabs=rep('.',nrow(pca_data)))

rotmat = pca$rotation

pca_data$protein = factor(pca_data$protein)

result = data.frame()

for(protein in levels(pca_data$protein))
{	
	#protein = '1RD8_AB'
	rotated_pca_data = pca_data[pca_data$protein==protein,c(1,2)]
	rotated_pca_data = cbind(rotated_pca_data,as.matrix(pca_data[pca_data$protein==protein,c(-1,-2)]) %*% rotmat)
	
	# CALCULATE CORRELATIONS WITH omega
		
		x = cor.test( rotated_pca_data$omega, rotated_pca_data$PC1 )
		r.PC1 = x$estimate
		p.PC1 = x$p.value

		x = cor.test( rotated_pca_data$omega, rotated_pca_data$PC2 )
		r.PC2 = x$estimate
		p.PC2 = x$p.value

		x = cor.test( rotated_pca_data$omega, rotated_pca_data$PC3 )
		r.PC3 = x$estimate
		p.PC3 = x$p.value

		x = cor.test( rotated_pca_data$omega, rotated_pca_data$PC4 )
		r.PC4 = x$estimate
		p.PC4 = x$p.value
		
		x = cor.test( rotated_pca_data$omega, rotated_pca_data$PC5 )
		r.PC5 = x$estimate
		p.PC5 = x$p.value
		
		x = cor.test( rotated_pca_data$omega, rotated_pca_data$PC6 )
		r.PC6 = x$estimate
		p.PC6 = x$p.value

		x = lm( omega~PC1+PC2+PC3+PC4+PC5+PC6, data=rotated_pca_data )
		x2 = summary(x)
		
	# Combine correlation values in row for each protein and dump it in RESULT data.frame
		
		row = data.frame( protein=protein,
						  rsq.PC1 = r.PC1^2, p.PC1 = p.PC1, 
						  rsq.PC2 = r.PC2^2, p.PC2 = p.PC2, 
						  rsq.PC3 = r.PC3^2, p.PC3 = p.PC3, 
						  rsq.PC4 = r.PC4^2, p.PC4 = p.PC4, 
						  rsq.PC5 = r.PC5^2, p.PC5 = p.PC5, 
						  rsq.PC6 = r.PC6^2, p.PC6 = p.PC6, 
						  rsq.total = x2$adj.r.squared )
	  sum = sum(unlist(row[1,c(2,4,6,8,10,12)]))
		diff = row$rsq.total[1]-sum
		row = data.frame( row, sum, diff )
		result = rbind( result, row )
}
#row.names(result) = c()
write.csv( result, out_path, row.names=F )

index = names(result) %in% c('rsq.PC1', 'rsq.PC2', 'rsq.PC3', 'rsq.PC4', 'rsq.PC5', 'rsq.PC6') # columns we want to plot
index.p = names(result) %in% c('p.PC1', 'p.PC2', 'p.PC3', 'p.PC4', 'p.PC5', 'p.PC6') # columns that store significance

colors = c('red', 'blue', 'green', 'purple', 'orange3', 'darkgreen', 'black', 'gray', 'cyan2') #, 'darkred', 'darkgreen', 'bisque2')
proteins = c('1RD8_AB', '2FP7_B', '2JLY_A', '2Z83_A', '3GOL_A', '3LYF_A', '4AQF_B', '4GHA_A', '4IRY_A') #, '3GSZ_A', '3I5K_A', '2JLY_A_temp_50', '2JLY_A_temp_100', '2JLY_A_temp_200', '2JLY_A_temp_450')
labels = c('1RD8', '2FP7', '2JLY', '2Z83', '3GOL', '3LYF', '4AQF', '4GHA', '4IRY') #, '3GSZ_A', '3I5K_A', '2JLY_A_temp_50', '2JLY_A_temp_100', '2JLY_A_temp_200', '2JLY_A_temp_450')

pdf( "correlation_analysis/figures/PC_screen_omega.pdf", width=9.3, height=9.2, useDingbats=FALSE )

split.screen(c(2,1))

screen(1)
	par( mai=c(0.65, 0.65, 0.4, 0.05), mgp=c(2, 0.5, 0), tck=-0.03 )
	#par( mai=c(0.65, 0.65, 0.4, 0.2), mgp=c(2, 0.5, 0), tck=-0.03 )
	plot(0,xaxt='n',yaxt='n',bty='n',pch='',ylab=expression(paste(omega," explained variance (", r^2, ")")),xlab='Principal Component (PC)', xlim=c(1,6.2),ylim=c(-.1,.2))
	#minor.tick(nx=0, ny=4, tick.ratio=2)
	axis( 1,  # x axis
		at=c(1, 2, 3, 4, 5, 6),
		padj=c(0,0,0,0,0,0),
		c('PC1', 'PC2', 'PC3', 'PC4', 'PC5', 'PC6')
		)
	axis( 2,  # y axis
		at=c(-.1,-.05,.0,.05,0.1,0.15,0.2),
		c("-0.1","","0.0","","0.1","","0.2")
		)

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
	
	legend( 1.0, -0.03, labels[1:3], pch=19, col=colors[1:3], bty='n', cex=0.9)
	legend( 1.7, -0.03, labels[4:6], pch=19, col=colors[4:6], bty='n', cex=0.9)
	legend( 2.4, -0.03, labels[7:9], pch=19, col=colors[7:9], bty='n', cex=0.9)
	mtext('A', side = 3, at = 0.5, font=2, cex=1.2)

split.screen(c(1, 2), 2)    # split bottom half in two

screen(3)
	par( mai=c(0.65, 0.65, 0.4, 0.2), mgp=c(2, 0.5, 0), tck=-0.03 )
	par( mai=c(0.65, 0.65, 0.65, 0.65), mgp=c(2, 0.5, 0), tck=-0.03 )
	mybiplot(pca,
	       choices = c(1,2),
	       xlabs=rep('.',nrow(pca_data)),
           ylabs = c('designed entropy','RSA','MD RMSF','VAR(chi1)','iWCN','B factor'),
           xlab="PC 1",
           ylab="PC 2",
           xlim=c(-.09,.03),
           ylim=c(-.07,.07),
           col=c('grey','red')
		   )
	mtext('B', side = 3, at=-80, font=2, cex=1.2)
screen(4)
	par( mai=c(0.65, 0.65, 0.65, 0.65), mgp=c(2, 0.5, 0), tck=-0.03 )
	mybiplot(pca,
	       choices = c(2,3),
	       xlabs=rep('.',nrow(pca_data)),
           ylabs = c('designed entropy','RSA','MD RMSF','VAR(chi1)','iWCN','B factor'),
           xlab="PC 2",
           ylab="PC 3",
           xlim=c(-.08,.10),
		   ylim=c(-.07,.07),
		   col=c('grey','red')
		   )
	mtext('C', side = 3, at=-80, font=2, cex=1.2)

close.screen(all = TRUE)

dev.off()