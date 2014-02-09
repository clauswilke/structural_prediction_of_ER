# This R code performs Principal Component Regression on different variables correlating with entropy and outputs the component loadings and the percentage of the entropy's explained variance by each principal component in a CVS file.

# Amir Shahmoradi, Wednesday 10:07 PM, February 5 2014, Wilke Lab, ICMB, UT Austin

library( pls )

# source("input_data.r")

file_out <- T

out_path = 'C:/Users/Amir/Documents/GitHub/structural_prediction_of_ER/correlation_analysis/pcr/'

# modified summary function for class mvr

#	mysum <- function( g ){
#	  # function to normalize the components
#	  # prints actually the squared values, to get percentages
#	  f <- function(i, g) {round( g$projection[,i]*g$projection[,i]/(g$projection[,i]%*% g$projection[,i]), 3)}
#	  # print the summary
#	  h <- summary( g )
#	  htmp <- h
#	  # number of components
#	  comps = length( attr(h, "dimnames")[[2]] )
#	  # print summary with percentage differences
#	  for ( i in 2:comps ){
#	    h[[1,i]] <- htmp[[1,i]]-htmp[[1,i-1]]
#	    h[[2,i]] <- htmp[[2,i]]-htmp[[2,i-1]]
#	  }
#	  cat( "\nPercentage differences:\n" )
#	  print( round( h, 2 ), print.gap=2 )
#	
#	  cat( "\nPercentage contributions to components:\n" )
#	
#	  # print the normalized projection
#	  p = g$projection
#	  for ( i in 1:comps ){
#	    x <- p[,i]
#	    p[,i] <- x^2/t(x) %*% x
#	  }
#	  print( round( p, 3 ) )
#	
#	 	# print total variance contributions
#		cat( "\nTotal variance explained by each variable\n")
#		total_var_explained = as.array(0*(1:comps))
#		#print(dimnames(total_var_explained))
#		#dimnames(total_var_explained) <- paste('comp', 1:comps)
#		for (i in 1:comps) {
#			total_var_explained[i] <- sum(h[2,] %*% p[i,])
#		}
#		print( round(total_var_explained, 3) )
#		
#		#bars <- as.vector(h[[2,]]) %*% p
#		#bars
#		bars <- p
#		for (i in 1:comps) {
#			bars[,i] <- h[2,i] * p[,i]
#		}
#		for (i in 1:comps) {
#			bars[,i] <- sort(bars[,i])
#		}
#		bars
#	}

# Replace NAs with mean values in data:

data = rbind(data_1RD8_AB, data_2FP7_B, data_2JLY_A, data_2Z83_A, data_3GOL_A, data_3LYF_A, data_4AQF_B, data_4GHA_A, data_4IRY_A,
             data_3GSZ_A, data_3I5K_A) #, data_2JLY_A_temp_50, data_2JLY_A_temp_100, data_2JLY_A_temp_200, data_2JLY_A_temp_450)

data$protein = factor(data$protein)

for(protein in levels(data$protein))
{	
	pcr_data = data[data$protein==protein,]
	
	pcr_data$phi_var_md[is.na(pcr_data$phi_var_md)] = mean(na.omit(pcr_data$phi_var_md))
	pcr_data$psi_var_md[is.na(pcr_data$psi_var_md)] = mean(na.omit(pcr_data$psi_var_md))
	pcr_data$chi1_var_md[is.na(pcr_data$chi1_var_md)] = mean(na.omit(pcr_data$chi1_var_md))

	# This part is for the sake of better output naming of the variables
	
	rsa_avg_md  = scale(pcr_data$rsa_avg_md)
	wcn_avg_md  = scale(pcr_data$wcn_avg_md)
	phi_var_md  = scale(pcr_data$phi_var_md)
	psi_var_md  = scale(pcr_data$psi_var_md)
	chi1_var_md = scale(pcr_data$chi1_var_md)
	rmsf_avg_md = scale(pcr_data$rmsf_avg_md)
	bfca        = scale(pcr_data$bfca)
	desent      = scale(pcr_data$desent)
	
	g <- pcr( pcr_data$entropy ~ rsa_avg_md
							   + wcn_avg_md
							   + phi_var_md
							   + psi_var_md
							   + chi1_var_md
							   + rmsf_avg_md
							   + bfca
							   + desent
							   , y=T)
							   
	#bars <- mysum( g )

	p = g$projection	# projection matrix
	
	PEV = c()	# percentage of explained variance
	
	for ( i in 1:g$ncomp ){
						x <- p[,i]
						p[,i] <- x^2/t(x) %*% x
						PEV = cbind(PEV,100*cor(g$scores[,i],g$y)^2)
						#print (i)
						}
	
	row.names(PEV) = 'Explained Variance'
	
	row.names(p) = row.names(g$projection)
	
	output_file = paste(out_path,'pcr_',protein,'.csv',sep='')
	
	if (file_out) write.csv( rbind(round(g$projection,3),PEV) , output_file, row.names=T )
	
	#seven_colors = rev(c("red", "orange","yellow", "green", "blue", "darkblue", "violet", "chocolate1", 'black', 'gray', 'purple'))
	#
	#if (file_out) pdf(paste(dir_out, "pcr_test.pdf", sep=''), family="Times", 
	#	width=10, height=7) else par(font=1, family="serif")
	#barplot(bars, col=seven_colors)
	#if (file_out) dev.off() else par(def.par)
	
	# print( summary( lm( g$y ~ g$scores ) ) )
}