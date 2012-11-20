print.Bvs <-
function(x, ...){
  cat("Call:\n")
  print(x$call)
  cat("\nThis is the result for a model selection problem with ")
  cat(x$p-1)
  pp<-2^x$p-1
  cat(" covariates and ")
  cat(x$n)
  cat(" observations\n")
  cat("The potential covariates are:\n")
  cat(x$variables[-1])
  if(!is.null(x$time)){
    cat("\nIt has take ")
    cat(x$time)
    cat(" seconds to compute.\n")
  }
  if(x$method=="gibbs"){
    cat("\nAmong the visited models, the Highest posterior probability one is: \n")
    print(x$modelsprob)
  }else{
    if(pp<10){
      cat(paste("\nThe",pp,"most probable models and their probabilities are:\n",sep=" "))
      print(x$modelsprob)
    }else{
      cat("\nThe 10 most probable models and their probabilities are:\n")
      print(x$modelsprob[1:10,])
    }
    
  }
    
}
