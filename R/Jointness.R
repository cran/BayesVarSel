Jointness <- function(x,covariates="All"){
  if (!inherits(x, "Bvs")) 
    warning("calling Jointness(<fake-Bvs-object>) ...")


if(length(covariates)==1 && covariates=="All"){
joint<-x$jointinclprob
joint_LS1 <- x$jointinclprob
joint_LS2 <- x$jointinclprob
for(i in 1:dim(x$jointinclprob)[1]){
  for(j in i:dim(x$jointinclprob)[1]){
joint_LS1[i,j] <- x$jointinclprob[i,j]/(x$jointinclprob[i,i]+x$jointinclprob[j,j]-x$jointinclprob[i,j])

joint_LS1[j,i] <- x$jointinclprob[i,j]/(x$jointinclprob[i,i]+x$jointinclprob[j,j]-x$jointinclprob[i,j])

joint_LS2[i,j] <- x$jointinclprob[i,j]/(x$jointinclprob[i,i]+x$jointinclprob[j,j]-2*x$jointinclprob[i,j])

joint_LS2[j,i] <- x$jointinclprob[i,j]/(x$jointinclprob[i,i]+x$jointinclprob[j,j]-2*x$jointinclprob[i,j])
  }
  joint_LS2[i,i] <- NA
  }


jointness <- list(joint, joint_LS1, joint_LS2)
class(jointness) <- "jointness"

return(jointness)
}
  if(length(covariates)>1){
  
    
if(length(covariates)>2){stop("The number of covariates to obtain jointness measurements should be 2\n")}
if(length(covariates)<2){stop("The number of covariates to obtain jointness measurements should be 2\n")}
i<-0
j<-0

i <- which(names(x$jointinclprob)==covariates[1])
j <- which(names(x$jointinclprob)==covariates[2])

if(i==0 || j==0){
  stop("At least one of the covariates is not part of the analysis")
}

prob_joint <- x$jointinclprob[i,j]

joint_LS1 <- x$jointinclprob[i,j]/(x$jointinclprob[i,i]+x$jointinclprob[j,j]-x$jointinclprob[i,j])

joint_LS2 <- x$jointinclprob[i,j]/(x$jointinclprob[i,i]+x$jointinclprob[j,j]-2*x$jointinclprob[i,j])

jointness<-list(as.data.frame(c(prob_joint,joint_LS1,joint_LS2),row.names =c("prob_joint","joint_LS1","joint_LS2")), covariates)


names(jointness) <- "value"
class(jointness) <- "jointness"
jointness
}
}
