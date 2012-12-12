GibbsBvs <-
function(formula, data, prior.betas="Robust", prior.models="Constant",  n.iter, init.model="Null", n.burnin=50, time.test=TRUE){
result<- list()
#Get a tempdir as working directory
wd<- tempdir()
#remove possibly existing files:
unlink(paste(wd,"*",sep="/"))

#Create the files with the design matrix and the dependent values
lm.obj = lm(formula, data, y=TRUE, x=TRUE)
Y<- lm.obj$y
X<- lm.obj$x
namesx<- dimnames(X)[[2]]
namesx[1]<- "Intercept"
#Dimension of the full model
p<- dim(X)[2]
#Number of observations
n<- dim(X)[1]

write(Y, ncolumns=1, file=paste(wd,"/Dependent.txt",sep=""))
write(t(X), ncolumns=p, file=paste(wd,"/Design.txt",sep=""))

#The initial model:
if (is.character(init.model)==TRUE){
	im<- substr(tolower(init.model),1,1)
	if (im=="n"){init.model<- c(1,rep(0,p-1))}
	if (im=="f"){init.model<- rep(1,p)}
	if (im=="r"){init.model<- c(1,rbinom(n=p-1,size=1,prob=.5))}
}
init.model<- as.numeric(init.model>0)
if (init.model[1]!=1){stop("Initial model must contain the intercept\n")}	

write(init.model, ncolumns=1, file=paste(wd,"/initialmodel.txt",sep=""))

#Info:
cat("Info. . . .\n")
cat("Most complex model has ",p-1,"covariates plus the intercept\n")
cat("The problem has a total of", 2^(p-1), "competing models\n")
cat("Of these,", n.burnin+n.iter, "are sampled with replacement\n")
cat("Then,", n.iter, "are used to construct the summaries\n")


#prior for betas:
pfb<- substr(tolower(prior.betas),1,1)
if (pfb!="g" && pfb!="r" && pfb!="z" && pfb!="l") stop("I am very sorry: prior for betas no valid\n")
#prior for model space:
pfms<- substr(tolower(prior.models),1,1)
if (pfms!="c" && pfms!="s") stop("I am very sorry: prior for model space no valid\n")

method<- paste(pfb,pfms,sep="")

#The previous test (for time)
estim.time<- 0
if(p<=20){
  warning("The number of variables is small enough to visit every model. Consider Bvs (or pBvs for its parallel version).\n")
}
if (time.test&&p>20){
	cat("Time test. . . .\n")
result<- switch(method,
	"gc"=.C("GibbsgConst", as.character(""), as.integer(n), as.integer(p), 
		as.integer(49),as.character(wd),as.integer(1), as.double(estim.time)),
	"gs"=.C("GibbsgSB", as.character(""), as.integer(n), as.integer(p), 
		as.integer(49),as.character(wd),as.integer(1), as.double(estim.time)),
	"rc"=.C("GibbsRobustConst", as.character(""), as.integer(n), as.integer(p), 
		as.integer(49),as.character(wd),as.integer(1), as.double(estim.time)),
	"rs"=.C("GibbsRobustSB", as.character(""), as.integer(n), as.integer(p), 
		as.integer(49),as.character(wd),as.integer(1), as.double(estim.time)),
	"lc"=.C("GibbsLiangConst", as.character(""), as.integer(n), as.integer(p), 
		as.integer(49),as.character(wd),as.integer(1), as.double(estim.time)),
	"ls"=.C("GibbsLiangSB", as.character(""), as.integer(n), as.integer(p), 
		as.integer(49),as.character(wd),as.integer(1), as.double(estim.time)),
	"zc"=.C("GibbsZSConst", as.character(""), as.integer(n), as.integer(p), 
		as.integer(49),as.character(wd),as.integer(1), as.double(estim.time)),
	"zs"=.C("GibbsZSSB", as.character(""), as.integer(n), as.integer(p), 
		as.integer(49),as.character(wd),as.integer(1), as.double(estim.time)))
	 
	estim.time<- result[[7]]*(n.burnin+n.iter)/(60*50) 
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

#Call the corresponding function:
result<- switch(method,
	"gc"=.C("GibbsgConst", as.character(""), as.integer(n), as.integer(p), 
		as.integer(n.iter),as.character(wd),as.integer(n.burnin), as.double(estim.time)),
	"gs"=.C("GibbsgSB", as.character(""), as.integer(n), as.integer(p), 
		as.integer(n.iter),as.character(wd),as.integer(n.burnin), as.double(estim.time)),
	"rc"=.C("GibbsRobustConst", as.character(""), as.integer(n), as.integer(p), 
		as.integer(n.iter),as.character(wd),as.integer(n.burnin), as.double(estim.time)),
	"rs"=.C("GibbsRobustSB", as.character(""), as.integer(n), as.integer(p), 
		as.integer(n.iter),as.character(wd),as.integer(n.burnin), as.double(estim.time)),
	"lc"=.C("GibbsLiangConst", as.character(""), as.integer(n), as.integer(p), 
		as.integer(n.iter),as.character(wd),as.integer(n.burnin), as.double(estim.time)),
	"ls"=.C("GibbsLiangSB", as.character(""), as.integer(n), as.integer(p), 
		as.integer(n.iter),as.character(wd),as.integer(n.burnin), as.double(estim.time)),
	"zc"=.C("GibbsZSConst", as.character(""), as.integer(n), as.integer(p), 
		as.integer(n.iter),as.character(wd),as.integer(n.burnin), as.double(estim.time)),
	"zs"=.C("GibbsZSSB", as.character(""), as.integer(n), as.integer(p), 
		as.integer(n.iter),as.character(wd),as.integer(n.burnin), as.double(estim.time)))
	   
  time <- result[[7]]
  
  


models <- as.vector(read.table(paste(wd,"/MostProbModels",sep="")))
incl <- as.vector(read.table(paste(wd,"/InclusionProb",sep="")))
joint <- as.matrix(read.table(paste(wd,"/JointInclusionProb",sep="")))
dimen <- as.vector(read.table(paste(wd,"/ProbDimension",sep="")))
betahat<- as.vector(read.table(paste(wd,"/betahat",sep="")))
#Highest probability model
mod.mat <- as.data.frame(t(models))

names(mod.mat)<-namesx
row.names(mod.mat) <- ""
varnames.aux<-rep("",p)
varnames.aux[mod.mat[1,1:p]==1]<-"*"
mod.mat[1,1:p] <- varnames.aux

  
inclusion <- data.frame(Inclusion_Prob= incl[-1,])
row.names(inclusion) <- namesx[-1] 

result <- list()
result$time<-time
result$lm<-lm.obj
result$variables <- namesx
result$n <- n
result$p <- p
result$HPMbin <- as.data.frame(models[-1,])
names(result$HPMbin)<-"bin.mod"
result$modelsprob<- mod.mat
result$inclprob<- inclusion
result$jointinclprob<- data.frame(joint[2:p,2:p],row.names=namesx[-1])
names(result$jointinclprob)<- namesx[-1]
result$postprobdim<- data.frame(dimen, row.names=1:p)
names(result$postprobdim) <- "Prob"
result$betahat <- betahat
result$call <- match.call()
result$method <- "gibbs"
class(result)<- "Bvs"
result
}