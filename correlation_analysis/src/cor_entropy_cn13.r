# This function calculates the correlation between the two variables entropy and RSA for each protein structure in the study and outputs the result into a CSV file.
# Amir Shahmoradi, Wednesday 3:06 PM, Jan 15, 2015, Wilke Lab, ICMB, UT Austin

# source("input_data.r")

data = rbind(data_1RD8, data_2Z83, data_3GOL, data_3GSZ, data_3I5K, data_4AQF_B, data_4GHA_A, data_4IRY,
             data_2JLY, data_2JLY_temp_50, data_2JLY_temp_100, data_2JLY_temp_200, data_2JLY_temp_450)
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

    row = data.frame( protein=protein, -r.cn13_cr, p.cn13_cr, -r.cn13_avg_md, p.cn13_avg_md, -r.cn13_var_md, p.cn13_var_md )
    result = rbind( result, row )
    print( protein )
}
row.names(result) = c()
write.csv( result, "cor_entropy_cn13.csv", row.names=F )


index = names(result) %in% c("X.r.cn13_cr", "X.r.cn13_avg_md", "X.r.cn13_var_md") # columns we want to plot

colors = c('red', 'blue', 'green', 'purple', 'orange', 'yellow', 'black', 'gray', 'red', 'blue', 'green', 'purple')
proteins = c('p1', 'p2', 'p3', 'p1', 'p2', 'p3', 'p1', 'p2', 'p3')

pdf( "test.pdf", width=4.5, height=4, useDingbats=FALSE )
par( mai=c(0.65, 0.65, 0.1, 0.05), mgp=c(2, 0.5, 0), tck=-0.03 )
plot(0,xaxt='n',yaxt='n',bty='n',pch='',ylab='Correlation',xlab='Variable', xlim=c(1,3),ylim=c(-.5,.5))
axis( 2 ) # y axis
axis( 1,
      at=c(1, 2, 3),
      c("CR", "Avg MD", "Var MD"))

for( i in 1:nrow(result) )
{
    row = unlist( result[i,index] )
    x = 1:sum(index)
    points( x, row, pch=19, col=colors[i])
    lines( x, row, col=colors[i] )
}

legend( 1.5, -.2, proteins[1:5], pch=19, col=colors[1:5], bty='n' )
legend( 2.2, -.2, proteins[6:10], pch=19, col=colors[6:10], bty='n' )
dev.off()