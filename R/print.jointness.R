print.jointness <- function(x, ...){
  if(length(x)==2){
    
    cat("---------\n")
    cat(paste("The joint inclusion probability for", x[[2]][1], "and", x[[2]][2], "is: ",round(x[[1]][1,1],2),"\n", sep =" "))
    cat("---------\n")
    cat(paste("The ratio between the probability of including both covariates and the probability of including at least one of then is: ",round(x[[1]][2,1],2),"\n", sep=""))
    cat("---------\n")
    cat(paste("The probability of including both covariates together is",round(x[[1]][3,1],2), "times the probability of including one of them alone \n",sep=" "))
  }
  if(length(x)==3){
    cat("---------\n")
    cat(paste("The joint inclusion probability for All covariates are \n", sep =" "))
    print(x[[1]])
    cat("---------\n")
    cat(paste("The ratio between the probability of including two covariates together and the probability of including at least one of them is: \n", sep=""))
    print(x[[2]])
    cat("---------\n")
    cat(paste("The ratio between the probability of including two covariates together and the probability of including one of them alone is: \n",sep=" "))
    print(x[[3]])
    
  }
}

