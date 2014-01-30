setwd('~/Desktop/Amir_Project/correlation_testing/')

get.rsas <- function(name) {
  return(read.table(name, sep=','))
}

get.omega <- function(fubar, this.map) {
  evo <- read.table(fubar, head=T, sep=',')
  map <- read.table(this.map, sep='\t', head=T)
  
  omega <- evo$beta[!is.na(map$pdb_pos)]/evo$alpha[!is.na(map$pdb_pos)]
  aas <- map$pdb_aa[!is.na(map$pdb_pos)]
  return.values <- data.frame(aa=c(aas), omega=c(omega))
  return(return.values)
}

rsa.1rd8 <- get.rsas('1rd8/1rd8.rsa')
rsa.2bmf <- get.rsas('2bmf/2bmf.rsa')
rsa.2fp7 <- get.rsas('2fp7/2fp7.rsa')
rsa.2hmi <- get.rsas('2hmi/2hmi.rsa')
rsa.2z83 <- get.rsas('2z83/2z83.rsa')
rsa.3tg6 <- get.rsas('3tg6/3rg6.rsa')
rsa.4akl <- get.rsas('4akl/4akl.rsa')
rsa.4gh9 <- get.rsas('4gh9/4gh9.rsa')


omega.1rd8 <- get.omega('1rd8/fubar_H1N1_HA.csv', '1rd8/1rd8A.map')
omega.2bmf <- get.omega('2bmf/fubar_dengue_ns3.csv', '2bmf/2bmf.map')

omega.2fp7A <- get.omega('2fp7/fubar_WestNile_ns2b.csv', '2fp7/2fp7A.map')
omega.2fp7B <- get.omega('2fp7/fubar_WestNile_ns2b.csv', '2fp7/2fp7B.map')
omega.2fp7 <- c(omega.2fp7A, omega.2fp7B)

omega.2hmi <- get.omega('2hmi/fubar_HIV_REV.csv', '2hmi/2hmiB.map')
omega.2z83 <- get.omega('2z83/fubar_JEV.csv', '2z83/2z83A.map')
omega.3tg6 <- get.omega('3tg6/fubar_H1N1_NP.csv', '3tg6/3tgcA.map')
omega.4akl <- get.omega('4akl/fubar_crimean_congo.Nucleoprotein.csv', '4akl/4aklA.map')
omega.4gh9 <- get.omega('4gh9/fubar_marburg.vp35.csv', '4gh9/4gh9A.map')

cor.1rd8 <- cor.test(rsa.1rd8$V2, omega.1rd8$omega)
cor.2bmf <- cor.test(rsa.2bmf$V2, omega.2bmf$omega)
#cor.2fp7 <- cor.test(rsa.2fp7, omega.2fp7)
#cor.2hmi <- cor.test(rsa.2hmi, omega.2hmi)
cor.2z83 <- cor.test(rsa.2z83$V2, omega.2z83$omega)
cor.3tg6 <- cor.test(rsa.3tg6$V2, omega.3tg6$omega)
cor.4akl <- cor.test(rsa.4akl$V2, omega.4akl$omega)
cor.4gh9 <- cor.test(rsa.4gh9$V2, omega.4gh9$omega)