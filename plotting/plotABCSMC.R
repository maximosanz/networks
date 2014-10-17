
library(locfit)
library(RColorBrewer)

d <- read.table('', header=F)
w <- read.table('', header=F)

nPar = 2
parNames = c('','')



a = c(0,0) # Lower limit of each parameter
b = c(1,1) # Upper limit


par(mfrow=c(nPar,nPar), mar=c(2,2,1,1))

for(j in 1:nPar){
  for(i in 1:nPar){	
    if(i==j){ 
      plot(density(d[,i], weights=w[,1]), main=parNames[i], xlab="",ylab="", col="red",xlim=c(a[i],b[i]))
    }
    else{
      fit <- locfit( ~ lp(d[,i],d[,j], nn=0.7,h=0.4,scale=T),data=d, weights=w[,1], xlim = c(c(0,0),c(1,1)), ev = lfgrid(10,c(a[i],a[j]),c(b[i],b[j])) )
      plot(fit,type="image",col=brewer.pal(9,"YlOrBr")[1:9],axes=F)
    }
  }
}