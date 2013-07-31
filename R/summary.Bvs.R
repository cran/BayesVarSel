summary.Bvs <-
function(object,...){
  z <- object
  p <- z$p
  if (!inherits(object, "Bvs")) 
    warning("calling summary.Bvs(<fake-Bvs-object>) ...")
ans<-list()
ans$coefficients <- z$betahat
dimnames(ans$coefficients) <- list(names(z$lm$coefficients),"Estimate")

HPM <- z$HPMbin
MPM <- as.numeric(z$inclprob>=0.5)
astHPM <- matrix(" ",ncol=1,nrow=(p-1))
astMPM <- matrix(" ",ncol=1,nrow=(p-1))
astHPM[HPM==1] <- "*"
astMPM[MPM==1] <- "*"
astHPM <- rbind("*",astHPM)
astMPM <- rbind("*",astMPM)

incl.prob<-z$inclprob
incl.prob <- rbind(1,incl.prob)
summ.Bvs <- cbind(incl.prob,astHPM,astMPM)
dimnames(summ.Bvs)<- list(names(z$lm$coefficients),c("Incl.prob.","HPM","MPM"))

ans$summary<-summ.Bvs  
ans$method <- z$method
ans$call<-z$call
class(ans) <- "summary.Bvs"
ans

}
