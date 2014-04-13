# This code calculates the Spearman correlation coefficient between RSA as measured from Crystal Structures and the average of MD simulation trajectories.
# Amir Shahmoradi, Friday 6:53 PM, April 11, 2014, ICMB, UT Austin

setwd('C:/Users/Amir/Documents/GitHub/structural_prediction_of_ER/')

# source('~/GitHub/structural_prediction_of_ER/correlation_analysis/src/input_data.r')

data = rbind(data_1RD8_AB, data_2FP7_B, data_2JLY_A, data_2Z83_A, data_3GOL_A, data_3LYF_A, data_4AQF_B, data_4GHA_A, data_4IRY_A) #,data_3GSZ_A, data_3I5K_A, data_2JLY_A_temp_50, data_2JLY_A_temp_100, data_2JLY_A_temp_200, data_2JLY_A_temp_450)

data$protein = factor(data$protein)

result = data.frame()

for(protein in levels(data$protein))
{
    d = data[data$protein==protein,]
    
    # RSA:
    x = cor.test( d$rsa_cr, d$rsa_avg_md, method="spearman", na.action="na.omit" )
    r.rsa = x$estimate
    p.rsa = x$p.value

    # CN13:
    x = cor.test( d$cn13_cr, d$cn13_avg_md, method="spearman", na.action="na.omit" )
    r.cn13 = x$estimate
    p.cn13 = x$p.value

    # WCN:
    x = cor.test( d$wcn_cr, d$wcn_avg_md, method="spearman", na.action="na.omit" )
    r.wcn = x$estimate
    p.wcn = x$p.value

    # RMSF:
    if (length(na.omit(d$rmsfHS))>=3)
	{
		x = cor.test( d$rmsf_avg_md, d$rmsfHS, method="spearman", na.action="na.omit" )
		r.rmsf = x$estimate
		p.rmsf = x$p.value
	}
	else
	{
		r.rmsf = NA
		p.rmsf = NA
	}	
	
    row = data.frame( protein = protein, r.rsa = r.rsa, p.rsa = p.rsa, r.cn13 = r.cn13, p.cn13 = p.cn13, r.wcn = r.wcn, p.wcn = p.wcn, r.rmsf = r.rmsf, p.rmsf = p.rmsf)
    result = rbind( result, row )
    print( protein )
}

row.names(result) = c()

row = data.frame( protein = 'MEAN',
                  r.rsa  = mean(result$r.rsa, na.rm = TRUE), p.rsa = 'NA',
				  r.cn13 = mean(result$r.cn13, na.rm = TRUE), p.cn13 = 'NA',
				  r.wcn  = mean(result$r.wcn, na.rm = TRUE), p.wcn = 'NA',
				  r.rmsf = mean(result$r.rmsf, na.rm = TRUE), p.rmsf = 'NA' )
result = rbind( result, row )

row = data.frame( protein = 'STDEV',
                  r.rsa  = sd(result$r.rsa, na.rm = TRUE), p.rsa = 'NA',
				  r.cn13 = sd(result$r.cn13, na.rm = TRUE), p.cn13 = 'NA',
				  r.wcn  = sd(result$r.wcn, na.rm = TRUE), p.wcn = 'NA',
				  r.rmsf = sd(result$r.rmsf, na.rm = TRUE), p.rmsf = 'NA' )
result = rbind( result, row )
				  
				  
write.csv( result, "correlation_analysis/cor_tables/cor_MD_CR.csv", row.names=F )