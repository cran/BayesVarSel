print.summary.Bvs <-
function(x,...){
  cat("Estimates and Inclusion Probabilities:\n")
  cat("Note: If a covariate has a * in HPM(MPM) indicates that it is\n included in the Higest Probability Model (Median Probability Model)\n ")
  print(x$summary)
  if(x$method=="gibbs"){
    cat("Note 2: For this problem you used Gibbs sampling so the HPM reported\n is the most probable among the visited models. Also\n the given inclusion probabilities are\n estimations of the real ones.\n")
  }
}
