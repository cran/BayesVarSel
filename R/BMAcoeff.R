BMAcoeff<- function(x, n.sim=10000, method="svd"){
  #x is an object of class Bvs
  
  #simulations of the model averaging mixture of the objective posterior distributions of the regression
  #coefficients. Weights are taken from the previous Bvs or GibbsBvs run.
  
  #(next is for compatibility between notations for Intercept in lm and Bvs)
  if (!is.null(x$lmnull)) {if(colnames(x$lmnull$x)[1]=="(Intercept)"){colnames(x$lmnull$x)[1]<- "Intercept"}}
  if(colnames(x$lmfull$x)[1]=="(Intercept)"){colnames(x$lmfull$x)[1]<- "Intercept"}
  #name of dep var
  name.y<- colnames(x$lmfull$model[1])
  
  
  #results are given as a matrix bma.coeffs
  bma.coeffs<- matrix(0, nrow=n.sim, ncol=length(x$lmfull$coefficients))
  colnames(bma.coeffs)<- colnames(x$lmfull$x)
  
  #differentiate if method="full" (enumeration) or method="Gibbs"
  if (x$method=="full" | x$method=="parallel"){
    
    cat("\n")
    cat("Simulations obtained using the best", dim(x$modelsprob)[1], "models\n")
    cat("that accumulate",round(sum(x$modelsprob[,"prob"]),2), "of the total posterior probability\n")
    
    #draw n.sim models with replacement and with weights proportional
    #to their posterior prob:
    models<- sample(x=dim(x$modelsprob)[1], size=n.sim, replace=TRUE, prob=x$modelsprob[,"prob"])
    #table of these models
    t.models<- table(models)
    cs.tmodels<- cumsum(t.models)
    
    X<- x$lmfull$x
    
    for (iter in 1:length(t.models)){
      #rMD is model drawn (a number between 1 and n.keep)
      rMD<- as.numeric(names(t.models)[iter]); howmany<- t.models[iter]
      
      #covs in that model rMD (apart from the fixed ones)
      covsrMD<- names(x$modelsprob[rMD,])[x$modelsprob[rMD,]=="*"]
      #the data in model drawn (dependent variable in first column)
      datarMD<- as.data.frame(cbind(x$lmfull$model[,name.y],X[,covsrMD]))
      
      colnames(datarMD)<- c(name.y,covsrMD)
      #now add the fixed.cov if any
      if (!is.null(x$lmnull)) {datarMD<- cbind(datarMD,x$lmnull$x)}
				
			#remove rare characters because a double application of lm command
			#(happens for instance in covariates with ":")
			colnames(datarMD)<- gsub("`","",colnames(datarMD))
      
			#formula for that model
      formMD<- as.formula(paste(name.y,"~.-1",sep=""))
      
      #fit
      fitrMD<- lm(formula=formMD, data=as.data.frame(datarMD),qr=TRUE)
      
      #simulated value for the beta:
      #simple version
      #rcoeff<- rmvnorm(n=howmany, mean=fitrMD$coefficients, sigma=vcov(fitrMD))
      #exact version
      Rinv <- qr.solve(qr.R(fitrMD$qr))
      iXtX<- Rinv%*%t(Rinv)				
      Sigma<- sum(fitrMD$residuals*datarMD[,name.y])*iXtX/fitrMD$df
      rcoeff<- rmvt(n=howmany, sigma=Sigma, df=fitrMD$df, delta=fitrMD$coefficients, type="shifted", method=method)				
#      if(sum(names(fitrMD$coefficients)%in%colnames(bma.coeffs))!= length(names(fitrMD$coefficients))){
#        stop("The names of your covariates may contain a non-standar charater ($,',:,ect.). Please check.")
#      }
        
      bma.coeffs[(max(cs.tmodels[iter-1],0)+1):cs.tmodels[iter],names(fitrMD$coefficients)]<- rcoeff
      
    }
    
  }
  if (x$method=="gibbs"){
    
    cat("\n")
    cat("Simulations obtained using the ",dim(x$modelslogBF)[1]," sampled models.\n") 
    cat("Their frequencies are taken as the true posterior probabilities\n")
    
    #draw n.sim models with replacement (weights are implicitly proportional
    #to their posterior prob since repetitions are included in the sample):
    models<- sample(x=dim(x$modelslogBF)[1], size=n.sim, replace=TRUE)
    #table of these models
    t.models<- table(models)
    cs.tmodels<- cumsum(t.models)
    
    X<- x$lmfull$x
    
    for (iter in 1:length(t.models)){
      #rMD is model drawn (a number between 1 and n.keep)
      rMD<- as.numeric(names(t.models)[iter]); howmany<- t.models[iter]
      
      #covs in that model rMD (apart from the fixed ones)
      covsrMD<- names(x$modelslogBF[rMD,])[x$modelslogBF[rMD,]=="1"]
      #the data in model drawn (dependent variable in first column)
      datarMD<- as.data.frame(cbind(x$lmfull$model[,name.y],X[,covsrMD]))
      
      colnames(datarMD)<- c(name.y,covsrMD)
      #now add the fixed.cov if any
      if (!is.null(x$lmnull)) {datarMD<- cbind(datarMD,x$lmnull$x)}
		
			#remove rare characters because a double application of lm command
			#(happens for instance in covariates with ":")
			colnames(datarMD)<- gsub("`","",colnames(datarMD))
      #formula for that model
      formMD<- as.formula(paste(name.y,"~.-1",sep=""))
      
      #fit
      fitrMD<- lm(formula=formMD, data=as.data.frame(datarMD), qr=TRUE)
      
      #simulated value for the beta:
      #simple version
      #rcoeff<- rmvnorm(n=howmany, mean=fitrMD$coefficients, sigma=vcov(fitrMD))
      #exact version
      Rinv <- qr.solve(qr.R(fitrMD$qr))
      iXtX<- Rinv%*%t(Rinv)				
      Sigma<- sum(fitrMD$residuals*datarMD[,name.y])*iXtX/fitrMD$df
      rcoeff<- rmvt(n=howmany, sigma=Sigma, df=fitrMD$df, delta=fitrMD$coefficients, type="shifted", method=method)				
#      if(sum(names(fitrMD$coefficients)%in%colnames(bma.coeffs))!= length(names(fitrMD$coefficients))){
#        stop("The names of your covariates may contain a non-standar charater ($,',:,ect.). Please check and replace.")
#      }
      
      bma.coeffs[(max(cs.tmodels[iter-1],0)+1):cs.tmodels[iter],names(fitrMD$coefficients)]<- rcoeff
      
    }
  }
  #A plot to show the users potential multimodalities
  
  bma.density <- function(x,name="x"){
    x.no0 <- x[x!=0]; igual <- sum(x==0); total=length(x)
    if(igual==total) print(paste(name,"is not present in any of the considered models"))
    if(igual!=total){
      dens <- density(x.no0)
      dens.mod <- dens$y*(1-(igual/total))
      plot(dens$x,dens.mod,type="l",yaxt="n",xlab="",ylab="",main=name,col="blue")
      text(x=quantile(dens$x,probs = 0.15),y=mean(dens.mod),round(1-(igual/total),3))
      #print("The numeric value in the picture, if any, indicates the probability of that variable to be in the considered models")
    }
    
  }
  class(bma.coeffs) <- "bma.coeffs"
  dx<- dim(bma.coeffs)[2]
  nc<- round(sqrt(dx)); nr<- ceiling(sqrt(dx))
  par(mfrow=c(nr,nc))
  for (i in 1:dx){
    par(mar=c(2,1.5,1.5,1))
   histBMA(bma.coeffs,colnames(bma.coeffs)[i],text=F)		
    #plot(density(bma.coeffs[,i]), main=colnames(bma.coeffs)[i],yaxt="n",col=4,xlab="",ylab="")
  }
  
  par(mfrow=c(1,1))
  par(mar=c(5, 4, 4, 2) + 0.1)
  return(bma.coeffs)
  
}
