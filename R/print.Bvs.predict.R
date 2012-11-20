print.Bvs.predict <-
function(x,...){
  if (!inherits(x, "Bvs.predict")) 
    warning("calling print.Bvs.predict(<fake-Bvs.predict-object>) ...")
  cat("Call:\n")
  print(x$call)
  cat(paste("\nPrediction for",x$nn,"observations:\n",sep=" "))
  print(x$prediction)
  cat("\n With desing matrix:\n")
  print(x$X)
  
  
}
