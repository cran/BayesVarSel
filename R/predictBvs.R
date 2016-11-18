predictBvs<- function(x, newdata, n.sim=10000){
	#x is an object of class Bvs

	#simulations of the model averaging mixture of the objective posterior predictive distributions. 
	#Weights are taken from the previous Bvs or GibbsBvs run.
	
	#(next is for compatibility between notations for Intercept in lm and Bvs)
	if (!is.null(x$lmnull)) {if(colnames(x$lmnull$x)[1]=="(Intercept)"){colnames(x$lmnull$x)[1]<- "Intercept"}}
	if (colnames(x$lmfull$x)[1]=="(Intercept)"){colnames(x$lmfull$x)[1]<- "Intercept"}
  #name of dep var
  name.y<- colnames(x$lmfull$model[1])
	
	if(!is.data.frame(newdata)){stop("newdata must be a data.frame\n")}
	
	#Now add the intercept if needed
	if (colnames(x$lmfull$x)[1]=="Intercept"){newdata$Intercept<- 1}
	
	#results are given as a matrix rpredictions
	rpredictions<- matrix(0, nrow=n.sim, ncol=dim(newdata)[1])

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
				#formula for that model
				formMD<- as.formula(paste(name.y,"~.-1",sep=""))

				#fit
				fitrMD<- lm(formula=formMD, data=as.data.frame(datarMD), qr=TRUE, x=TRUE)
				
				#corresponding new design for that model:
				newdataMD<- newdata[ ,names(fitrMD$coefficients)]
				
				fnx<- apply(X=as.matrix(newdataMD), MARGIN=1,
						  FUN=function(x, DM){
								Xstar<- rbind(DM, x)
								qrXstar<- qr(Xstar); Rinv <- qr.solve(qr.R(qrXstar))
								iXstartXstar<- Rinv%*%t(Rinv)
								1-t(x)%*%iXstartXstar%*%x},
							DM=as.matrix(fitrMD$x))
					
				sigmas<- sum(fitrMD$residuals*datarMD[,name.y])/(fnx*fitrMD$df)
				#compute means:
				means<- as.matrix(newdataMD)%*%fitrMD$coefficients
				
				#simulated value for the new y's:
				if (dim(newdata)[1]==1){sigmaM<- as.matrix(sigmas, nr=1, nc=1)}
				if (dim(newdata)[1]>1){sigmaM<- diag(sigmas)}
				rpredictions[(max(cs.tmodels[iter-1],0)+1):cs.tmodels[iter],]<- 
				   rmvt(n=howmany, delta=means, sigma=sigmaM, df=fitrMD$df, type="shifted")

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
				#formula for that model
				formMD<- as.formula(paste(name.y,"~.-1",sep=""))

				#fit
				fitrMD<- lm(formula=formMD, data=as.data.frame(datarMD), qr=TRUE, x=TRUE)
				
				#corresponding new design for that model:
				newdataMD<- newdata[ ,names(fitrMD$coefficients)]
				
				#compute scales (one for each configuration)
				#(notation fnx from Bernardo's book)
				fnx<- apply(X=as.matrix(newdataMD), MARGIN=1,
						  FUN=function(x, DM){
								Xstar<- rbind(DM, x)
								qrXstar<- qr(Xstar); Rinv <- qr.solve(qr.R(qrXstar))
								iXstartXstar<- Rinv%*%t(Rinv)
								1-t(x)%*%iXstartXstar%*%x},
							DM=as.matrix(fitrMD$x))
				
				sigmas<- sum(fitrMD$residuals*datarMD[,name.y])/(fnx*fitrMD$df)
				#compute means:
				means<- as.matrix(newdataMD)%*%fitrMD$coefficients
				
				#simulated value for the new y's:
				if (dim(newdata)[1]==1){sigmaM<- as.matrix(sigmas, nr=1, nc=1)}
				if (dim(newdata)[1]>1){sigmaM<- diag(sigmas)}
				rpredictions[(max(cs.tmodels[iter-1],0)+1):cs.tmodels[iter],]<- 
				   rmvt(n=howmany, delta=means, sigma=sigmaM, df=fitrMD$df, type="shifted")
			
		}
	}

return(rpredictions)
	
}
