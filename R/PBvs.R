PBvs <-
function(formula, data, prior.betas="Robust", prior.models="Constant", n.keep, n.nodes=2){

require(snow)#package for parallel computation

cl <- makeCluster(n.nodes, type = "SOCK") 

#Get the tempdir as working directory
wd<- tempdir()
#remove possibly existing files:
unlink(paste(wd,"*",sep="/"))
#Create the files with the design matrix andf the dependent values
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


#prior for betas:
pfb<- substr(tolower(prior.betas),1,1)
if (pfb!="g" && pfb!="r" && pfb!="z" && pfb!="l") stop("I am very sorry: prior for betas no valid\n")
#prior for model space:
pfms<- substr(tolower(prior.models),1,1)
if (pfms!="c" && pfms!="s") stop("I am very sorry: prior for model space no valid\n")

#Check if the number of models is large enough.
if(n.keep>2^(p-1)/n.nodes)
  stop("The number of models to keep (n.keep) should be smaller than the total number of models divided by the number of nodes (n.nodes)")

#Check if the number of covariates is not too large
if (p>30){stop("Number of covariates too big. . . consider using GibbsBvs\n")}


method<- paste(pfb,pfms,sep="")

#Info:
cat("Info. . . .\n")
cat("Most complex model has ",p-1,"covariates plus the intercept\n")
cat("The problem has a total of", 2^(p-1), "competing models\n")
cat("Of these, the ", n.keep, "most probable (a posteriori) are kept\n")
cat("Working on the problem...please wait\n")

#Call the corresponding function:

estim.time<- 0
myfun<- function(name.start.end, method){
	switch(method,
	"gc"=.C("gConst", as.character(name.start.end[1]), as.integer(n), as.integer(p), as.integer(n.keep), 
	   as.integer(name.start.end[2]), as.integer(name.start.end[3]), as.character(wd), as.double(estim.time)),
	"gs"=.C("gSB", as.character(name.start.end[1]), as.integer(n), as.integer(p), as.integer(n.keep), 
	   as.integer(name.start.end[2]), as.integer(name.start.end[3]), as.character(wd), as.double(estim.time)),
	"rc"=.C("RobustConst", as.character(name.start.end[1]), as.integer(n), as.integer(p), as.integer(n.keep), 
	   as.integer(name.start.end[2]), as.integer(name.start.end[3]), as.character(wd), as.double(estim.time)),
	"rs"=.C("RobustSB", as.character(name.start.end[1]), as.integer(n), as.integer(p), as.integer(n.keep), 
	   as.integer(name.start.end[2]), as.integer(name.start.end[3]), as.character(wd), as.double(estim.time)),
	"lc"=.C("LiangConst", as.character(name.start.end[1]), as.integer(n), as.integer(p), as.integer(n.keep), 
	   as.integer(name.start.end[2]), as.integer(name.start.end[3]), as.character(wd), as.double(estim.time)),
	"ls"=.C("LiangSB", as.character(name.start.end[1]), as.integer(n), as.integer(p), as.integer(n.keep), 
	   as.integer(name.start.end[2]), as.integer(name.start.end[3]), as.character(wd), as.double(estim.time)),
	"zc"=.C("ZSConst", as.character(name.start.end[1]), as.integer(n), as.integer(p), as.integer(n.keep), 
	   as.integer(name.start.end[2]), as.integer(name.start.end[3]), as.character(wd), as.double(estim.time)),
	"zs"=.C("ZSSB", as.character(name.start.end[1]), as.integer(n), as.integer(p), as.integer(n.keep), 
	   as.integer(name.start.end[2]), as.integer(name.start.end[3]), as.character(wd), as.double(estim.time)))
}

#Load the library in the different nodes
clusterEvalQ(cl, library(BayesVarSel))

#Calculate how to distribute the model space through the nodes:
if (n.nodes<2) stop("At least 2 nodes are needed\n")

iterperproc<- round((2^(p-1)-1)/n.nodes)
if (n.keep>iterperproc) stop("Number of kept models should be smaller than the number of models per node\n")
distrib<- list()
for (i in 1:(n.nodes-1)){distrib[[i]]<- c(i,(i-1)*iterperproc+1, i*iterperproc)}
distrib[[n.nodes]]<- c(n.nodes,(n.nodes-1)*iterperproc+1, 2^(p-1)-1)

clusterApply(cl, distrib, myfun, method=method)  
#myfun(method=method, startend=c(1, 2^(p-1)-1))  
  

stopCluster(cl)  

##############Put together the results

#next is the prior probability for the null model Pr(M_0)=p_0/sum(p_j)
if (pfms=="c"){
	PrM0<- 1/2^(p-1)
	#the unnormalized prior prob for M0:
	p0<- 1
}

if (pfms=="s"){
	PrM0<- 1/(p+1)
	#the unnormalized prior prob for M0:
	p0<- 1
}

fPostProb<- paste(wd,"PostProb", sep="/")
fInclusionProb<- paste(wd,"InclusionProb", sep="/")
fMostProbModels<- paste(wd,"MostProbModels", sep="/")
fNormConstant<- paste(wd,"NormConstant", sep="/")
fNormConstantPrior<- paste(wd,"NormConstantPrior", sep="/")
fProbDimension<- paste(wd,"ProbDimension", sep="/")
fJointInclusionProb<- paste(wd,"JointInclusionProb", sep="/")
fBetahat<- paste(wd,"betahat", sep="/")

#Obtain the normalizing constant (say E) for the prior probabilities:
#Pr(Ml)=p_l/E
E<- 0
for (i in 1:n.nodes){
	E<- E+scan(file=paste(fNormConstantPrior,i,sep=""), n=1, quiet=T)
}
E<- E-(n.nodes-1)*p0

#Obtain the normalizing constant (say D) for the posterior probabilities:
#Pr(Ml|data)=B_{l0}*Pr(M_l)/D, where B_{l0}=m_l(data)/m_0(data)
D<- 0
for (i in 1:n.nodes){
	D<- D+scan(file=paste(fNormConstant,i,sep=""), n=1, quiet=T)*scan(file=paste(fNormConstantPrior,i,sep=""), n=1, quiet=T)
}
D<- (D-(n.nodes-1)*PrM0)/E


#Now obtain the n.keep most probable models
i<- 1
thisNormConstant<- scan(file=paste(fNormConstant,i,sep=""), n=1, quiet=T)
thisNormConstantPrior<- scan(file=paste(fNormConstantPrior,i,sep=""), quiet=T)
#next is the Bayes factor times the (unnormalized) prior for this model
#(see the main.c code to see how is the unnormalized prior). So, if
#the unnormalized prior is=1, then next is the Bayes factor*1 and so on
thisUnnorPostProb<- read.table(file=paste(fPostProb,i,sep=""))[[1]]*thisNormConstant*thisNormConstantPrior

thisMostProbModels<- read.table(file=paste(fMostProbModels,1,sep=""), colClasses="character")[[1]]

for (i in 2:n.nodes){
	readNormConstant<- scan(file=paste(fNormConstant,i,sep=""), n=1, quiet=T)
        readNormConstantPrior<- scan(file=paste(fNormConstantPrior,i,sep=""), quiet=T)
	readUnnorPostProb<- read.table(file=paste(fPostProb,i,sep=""))[[1]]*readNormConstant*readNormConstantPrior
	readMostProbModels<- read.table(file=paste(fMostProbModels,i,sep=""), colClasses="character")[[1]]

	jointUnnorPostProb<- c(readUnnorPostProb, thisUnnorPostProb) 
        jointModels<- c(readMostProbModels, thisMostProbModels)
	reorder<- order(jointUnnorPostProb, decreasing=T)
	thisUnnorPostProb<- jointUnnorPostProb[reorder[1:n.keep]]
	thisMostProbModels<- jointModels[reorder[1:n.keep]]
}


#The inclusion probabilities
accum.InclusionProb<- read.table(file=paste(fInclusionProb,i,sep=""))[[1]]*0
for (i in 1:n.nodes){
	readNormConstant<- scan(file=paste(fNormConstant,i,sep=""), n=1, quiet=T)
    readNormConstantPrior<- scan(file=paste(fNormConstantPrior,i,sep=""), quiet=T)
	accum.InclusionProb<- accum.InclusionProb+read.table(file=paste(fInclusionProb,i,sep=""))[[1]]*readNormConstant*readNormConstantPrior
}

accum.InclusionProb<- accum.InclusionProb/(D*E)

#The joint inclusion probs:
accum.JointInclusionProb<- as.matrix(read.table(file=paste(fJointInclusionProb,i,sep="")))*0
for (i in 1:n.nodes){
	readNormConstant<- scan(file=paste(fNormConstant,i,sep=""), n=1, quiet=T)
    readNormConstantPrior<- scan(file=paste(fNormConstantPrior,i,sep=""), quiet=T)
	accum.JointInclusionProb<- accum.JointInclusionProb+
	  as.matrix(read.table(file=paste(fJointInclusionProb,i,sep="")))*readNormConstant*readNormConstantPrior
}

accum.JointInclusionProb<- accum.JointInclusionProb/(D*E)
#-----

#The dimension probabilities
accum.ProbDimension<- read.table(file=paste(fProbDimension,i,sep=""))[[1]]*0
for (i in 1:n.nodes){
	readNormConstant<- scan(file=paste(fNormConstant,i,sep=""), n=1, quiet=T)
        readNormConstantPrior<- scan(file=paste(fNormConstantPrior,i,sep=""), quiet=T)
	accum.ProbDimension<- accum.ProbDimension+read.table(file=paste(fProbDimension,i,sep=""))[[1]]*readNormConstant*readNormConstantPrior
}

accum.ProbDimension<- accum.ProbDimension/(D*E)

betahat<- read.table(file=paste(fBetahat,i,sep=""))[[1]]*0
ac<- 0
for (i in 1:n.nodes){
    readNormConstant<- scan(file=paste(fNormConstant,i,sep=""), n=1, quiet=T)
    readNormConstantPrior<- scan(file=paste(fNormConstantPrior,i,sep=""), quiet=T)
    betahat<- betahat+read.table(file=paste(fBetahat,i,sep=""))[[1]]*readNormConstant*readNormConstantPrior
	ac<- ac+readNormConstant
	}
betahat<- betahat/(D*E)

write.table(file=fMostProbModels, thisMostProbModels, row.names=F, col.names=F)
write.table(file=fPostProb, thisUnnorPostProb/(D*E), row.names=F, col.names=F)
write.table(file=fInclusionProb, accum.InclusionProb, row.names=F, col.names=F)
write.table(file=fProbDimension, accum.ProbDimension, row.names=F, col.names=F)
write.table(file=fNormConstant, D, row.names=F, col.names=F)
write.table(file=fNormConstant, D, row.names=F, col.names=F)
write.table(file=fNormConstantPrior, E, row.names=F, col.names=F)
write.table(file=fBetahat, betahat, row.names=F, col.names=F)
write.table(file=fJointInclusionProb, accum.JointInclusionProb, row.names=F, col.names=F)

##############End of put together the results

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


tempdir <- wd
models <- as.vector(read.table(paste(tempdir,"/MostProbModels",sep="")))
prob <- as.vector(read.table(paste(tempdir,"/PostProb",sep="")))
incl <- as.vector(read.table(paste(tempdir,"/InclusionProb",sep="")))
joint <- as.matrix(read.table(paste(tempdir,"/JointInclusionProb",sep="")))
dimen <- as.vector(read.table(paste(tempdir,"/ProbDimension",sep="")))
betahat<- as.vector(read.table(paste(wd,"/betahat",sep="")))

#Most probable models
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

###############Para depurarse:
inclusion <- data.frame(Inclusion_Prob= incl[-1,])
row.names(inclusion) <- namesx[-1] 


result<-list()
result$time <- NULL 
result$lm<-lm.obj
result$variables <- namesx
result$n <- n
result$p <- p
result$HPMbin <- integer.base.b_C(models[1,1],(p-1))
result$modelsprob<- mod.mat
result$inclprob<- inclusion
result$jointinclprob<- data.frame(joint[2:p,2:p],row.names=namesx[-1])
names(result$jointinclprob)<- namesx[-1]
result$postprobdim<- data.frame(dimen, row.names=1:p)
names(result$postprobdim) <- "Prob"
result$betahat <- betahat
result$call <- match.call()
result$method <- "parallel"
class(result)<- "Bvs"
result
}
