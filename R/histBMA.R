histBMA <- function(x, covariate, n.breaks=100, text = TRUE, gray.0 = 0.6, gray.no0=0.8){
  
  #x must be an object of class bma.coeffs
  if (!inherits(x, "bma.coeffs")){stop("x must be an object of class bma.coeffs")}
  #covariate must be a variable in the analysis
  i <- which(colnames(x)==covariate)
  if (sum(colnames(x)==covariate)==0){stop(paste("covariate ",covariate," is not in the analysis",sep=""))}
  
  v<- x[,covariate]
  	
  igual <- sum(v==0)
  menor <- sum(v<0)
  mayor <- sum(v>0)
  total <- length(v)
  dist <-(max(v)-min(v))/(n.breaks)
  if(menor!=0 & mayor!=0&igual!=0){
    cortes <- c(seq(from=min(v),to = 0-dist/2,by=dist),0-dist/2,seq(0+dist/2,max(v)+dist,by=dist))
    cortes <- unique(cortes)
    h <- hist(v,breaks=cortes,plot=F)
    cutoff <- cut(h$breaks,breaks = c(-Inf,0-dist/2,0+dist/2,+Inf))
    color <- c(rep(gray(gray.no0),sum(as.numeric(cutoff)==1)-1),rep(gray(gray.0),1), rep(gray(gray.no0),sum(as.numeric(cutoff)==3)))
    plot(h,col=color,freq=FALSE,yaxt="n",main=covariate,xlab="",ylab="")
    if(text) text(x=0,y=h$density[h$mids==0]+ 0.025*h$density[h$mids==0],labels = round(igual/total,3))
  }
  if(menor==0&igual!=0&mayor!=0){
    cortes <- c(0-dist/2,0+dist/2,seq(from=0+dist/2,to = max(v)+dist,by=dist))
    cortes <- unique(cortes)
    h <- hist(v,breaks=cortes,plot=F)
    cutoff <- cut(h$breaks,breaks = c(-Inf,0-dist/2,0+dist/2,+Inf))
    color <- c(rep(gray(gray.no0),sum(as.numeric(cutoff)==1)-1),rep(gray(gray.0),1), rep(gray(gray.no0),sum(as.numeric(cutoff)==3)))
    plot(h,col=color,freq=FALSE,yaxt="n",xlab="",ylab="",main=covariate)
    if(text) text(x=0,y=h$density[h$mids==0]+ 0.025*h$density[h$mids==0],labels = round(igual/total,3))
  }
  if(mayor==0&igual!=0&menor!=0){
    cortes <- c(seq(from=min(v),to = 0-dist/2,by=dist),0-dist/2,0+dist/2)
    cortes <- unique(cortes)
    h <- hist(v,breaks=cortes,plot=F)
    cutoff <- cut(h$breaks,breaks = c(-Inf,0-dist/2,0+dist/2,+Inf))
    color <- c(rep(gray(gray.no0),sum(as.numeric(cutoff)==1)-1),rep(gray(gray.0),1), rep(gray(gray.no0),sum(as.numeric(cutoff)==3)))
    plot(h,col=color,freq=FALSE,yaxt="n",xlab="",ylab="",main=covariate)
    if(text) text(x=0,y=h$density[h$mids==0]+ 0.025*h$density[h$mids==0],labels = round(igual/total,3))
  }
  
  if(igual==0){
    hist(v,breaks=100,freq=FALSE,col=gray(gray.no0),yaxt="n",xlab="",ylab="",main=covariate)
  }
  if(igual==total){
    print(paste("The regression coefficient associated to",covariate,"is 0 in all the models used for BMA"))
  }
}





