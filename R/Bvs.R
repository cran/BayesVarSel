Bvs <-
function(formula, data, prior.betas="Robust", prior.models="Constant", n.keep, time.test=TRUE){
#Let's define the result 
result<- list()

#Get a tempdir as working directory
wd<- tempdir()
#remove all the previous documents in the working directory
unlink(paste(wd,"*",sep="/"))

#Create the files with the design matrix and the dependent values
lm.obj = lm(formula, data, y=TRUE, x=TRUE)
Y<- lm.obj$y
X<- lm.obj$x
namesx<- dimnames(X)[[2]]
namesx[1]<- "Intercept" #namesx contains the name of variables including the intercept

p<- dim(X)[2]#Dimension of the full model including the intercept

n<- dim(X)[1]#Number of observations

#check if the number of models to save is correct
if(n.keep>2^(p-1))
  stop(paste("The number of models to keep (",n.keep, ") is larger than the total number of models (",2^(p-1),")",sep=""))

#write the data files in the working directory
write(Y, ncolumns=1, file=paste(wd,"/Dependent.txt",sep=""))
write(t(X), ncolumns=p, file=paste(wd,"/Design.txt",sep=""))


#prior for betas:
pfb<- substr(tolower(prior.betas),1,1)
  #check if the selected option exists
if (pfb!="g" && pfb!="r" && pfb!="z" && pfb!="l") stop("I am very sorry: prior for betas no valid\n")
#prior for model space:
pfms<- substr(tolower(prior.models),1,1)
  #check if the selected option exists
if (pfms!="c" && pfms!="s") stop("I am very sorry: prior for model space no valid\n")

method<- paste(pfb,pfms,sep="")

#Info:
cat("Info. . . .\n")
cat("Most complex model has ",p-1,"covariates plus the intercept\n")
cat("The problem has a total of", 2^(p-1), "competing models\n")
cat("Of these, the ", n.keep, "most probable (a posteriori) are kept\n")

#check if the number of covariates is too big. 
if (p>30){stop("Number of covariates too big. . . consider using GibbsBvs\n")}

#The previous test (for time)
estim.time<- 0

if (time.test && p>=18){
	cat("Time test. . . .\n")
	result<- switch(method,
	"gc"=.C("gConst", as.character(""), as.integer(n), as.integer(p), as.integer(1), 
	        as.integer(2^(p-2)-1999), as.integer(2^(p-2)+2000), as.character(wd), as.double(estim.time)),
	"gs"=.C("gSB", as.character(""), as.integer(n), as.integer(p), as.integer(1), 
	        as.integer(2^(p-2)-1999), as.integer(2^(p-2)+2000), as.character(wd), as.double(estim.time)),
	"rc"=.C("RobustConst", as.character(""), as.integer(n), as.integer(p), as.integer(1), 
	        as.integer(2^(p-2)-1999), as.integer(2^(p-2)+2000), as.character(wd), as.double(estim.time)),
	"rs"=.C("RobustSB", as.character(""), as.integer(n), as.integer(p), as.integer(1), 
	        as.integer(2^(p-2)-1999), as.integer(2^(p-2)+2000), as.character(wd), as.double(estim.time)),
	"lc"=.C("LiangConst", as.character(""), as.integer(n), as.integer(p), as.integer(1), 
	        as.integer(2^(p-2)-1999), as.integer(2^(p-2)+2000), as.character(wd), as.double(estim.time)),
	"ls"=.C("LiangSB", as.character(""), as.integer(n), as.integer(p), as.integer(1), 
	        as.integer(2^(p-2)-1999), as.integer(2^(p-2)+2000), as.character(wd), as.double(estim.time)),
	"zc"=.C("ZSConst", as.character(""), as.integer(n), as.integer(p), as.integer(1), 
	        as.integer(2^(p-2)-1999), as.integer(2^(p-2)+2000), as.character(wd), as.double(estim.time)),
	"zs"=.C("ZSSB", as.character(""), as.integer(n), as.integer(p), as.integer(1), 
	        as.integer(2^(p-2)-1999), as.integer(2^(p-2)+2000), as.character(wd), as.double(estim.time)))
	 
	estim.time<- result[[8]]*2^(p-1)/(60*4000) 
	cat("The problem would take ", estim.time, "minutes (approx.) to run\n")
	ANSWER <- readline("Do you want to continue?(y/n) then press enter.\n")
	while (substr(ANSWER, 1, 1) != "n" & substr(ANSWER, 1, 1) !="y"){
	  ANSWER <- readline("")
	}
	
	if (substr(ANSWER, 1, 1) == "n")
	{
	  return(NULL)
	}
	
}
#if the answer is yes work on the problem
cat("Working on the problem...please wait.\n")


#if the answer is yes work on the problem
result<- switch(method,
	"gc"=.C("gConst", as.character(""), as.integer(n), as.integer(p), as.integer(n.keep), as.integer(1), as.integer(2^(p-1)-1), as.character(wd),  as.double(estim.time)),
	"gs"=.C("gSB", as.character(""), as.integer(n), as.integer(p), as.integer(n.keep), 
	   as.integer(1), as.integer(2^(p-1)-1), as.character(wd),  as.double(estim.time)),
	"rc"=.C("RobustConst", as.character(""), as.integer(n), as.integer(p), as.integer(n.keep), 
	   as.integer(1), as.integer(2^(p-1)-1), as.character(wd),  as.double(estim.time)),
	"rs"=.C("RobustSB", as.character(""), as.integer(n), as.integer(p), as.integer(n.keep), 
	   as.integer(1), as.integer(2^(p-1)-1), as.character(wd),  as.double(estim.time)),
	"lc"=.C("LiangConst", as.character(""), as.integer(n), as.integer(p), as.integer(n.keep), 
	   as.integer(1), as.integer(2^(p-1)-1), as.character(wd),  as.double(estim.time)),
	"ls"=.C("LiangSB", as.character(""), as.integer(n), as.integer(p), as.integer(n.keep), 
	   as.integer(1), as.integer(2^(p-1)-1), as.character(wd),  as.double(estim.time)),
	"zc"=.C("ZSConst", as.character(""), as.integer(n), as.integer(p), as.integer(n.keep), 
	   as.integer(1), as.integer(2^(p-1)-1), as.character(wd),  as.double(estim.time)),
	"zs"=.C("ZSSB", as.character(""), as.integer(n), as.integer(p), as.integer(n.keep), 
	   as.integer(1), as.integer(2^(p-1)-1), as.character(wd),  as.double(estim.time)))
	   
  time <- result[[8]]

#a function to transform the number of the model into a binary number.
integer.base.b_C<-function(x, k){
  #x is the number we want to express in binary
  #k is the number positions we need
  if(x==0)
    return(rep(0,k))
  else{
    ndigits <- (floor(logb(x, base=2))+1)
    res<- rep(0, ndigits)
    for(i in 1:ndigits){#i <- 1
      res[i] <- (x %% 2)
      x <- (x %/% 2)
    }
    return(c(res,rep(0,k-ndigits)))}
}

#read the files given by C
models <- as.vector(read.table(paste(wd,"/MostProbModels",sep=""),colClasses="numeric"))
prob <- as.vector(read.table(paste(wd,"/PostProb",sep=""),colClasses="numeric"))
incl <- as.vector(read.table(paste(wd,"/InclusionProb",sep=""),colClasses="numeric"))
joint <- as.matrix(read.table(paste(wd,"/JointInclusionProb",sep=""),colClasses="numeric"))
dimen <- as.vector(read.table(paste(wd,"/ProbDimension",sep=""),colClasses="numeric"))
betahat<- as.vector(read.table(paste(wd,"/betahat",sep=""),colClasses="numeric"))


# data.frame with Most probable models
mod.mat <- as.data.frame(cbind(t(rep(0,(p+1)))))

names(mod.mat)<-c(namesx,"prob")

N<-n.keep

for(i in 1:N){
  mod.mat[i,2:p]<-integer.base.b_C(models[i,1],(p-1))
  mod.mat[i,1]<-1
  varnames.aux<-rep("",p)
  varnames.aux[mod.mat[i,1:p]==1]<-"*"
  mod.mat[i,1:p]<-varnames.aux
}

mod.mat[,(p+1)]<-prob[]


inclusion <- data.frame(Inclusion_Prob= incl[-1,])#inclusion probabilities except for the intercept
row.names(inclusion) <- namesx[-1] 
                 
#the final result
result <- list()

result$time <- time
result$lm <- lm.obj
result$variables <- namesx
result$n <- n
result$p <- p
result$HPMbin <- integer.base.b_C(models[1,1],(p-1))
result$modelsprob <- mod.mat
result$inclprob <- inclusion

result$jointinclprob <- data.frame(joint[2:p,2:p],row.names=namesx[-1])
names(result$jointinclprob) <- namesx[-1]

result$postprobdim <- data.frame(dimen, row.names=1:p)
names(result$postprobdim) <- "Prob"

result$betahat <- betahat
result$call <- match.call()
result$method <- "full"
class(result)<- "Bvs"
result 
}
