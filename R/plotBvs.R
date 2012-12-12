plotBvs <-
function(x,option="dimension"){
  
  if (!inherits(x, "Bvs")) 
    warning("calling predict.Bvs(<fake-Bvs-object>) ...")
  p <- x$p
  auxtp<- substr(tolower(option),1,1)
  
  if (auxtp!="d" && auxtp!="j" && auxtp!="c" && auxtp!="n") 
    stop("I am very sorry: type of plot not specified\n")
  
  #dimension probabilities
  if (auxtp=="d"){
    par(mar=c(5, 4, 4, 2) + 0.1,mfrow=c(1,1))
   if(x$method=="gibbs"){
     barplot(x$postprobdim$Prob,main="Estimated Posterior Dimension Probabilities",xlab="Number of covariates in the model",ylab="Probability",names.arg=0:(p-1))
   }else{
    barplot(x$postprobdim$Prob,main="Posterior Dimension Probabilities",xlab="Number of covariates in the model",ylab="Probability",names.arg=0:(p-1))
    }
  }
  #Special function for printing the joint and conditional probabilities. 
  
  
  myImagePlot <- function(x, scale,...){
     x<- as.matrix(x)
    
    
    #What do we do with the diagonal? 
    diag(x)<- 0
    #What do we do with the NA's? 
    x[is.na(x)]<- 0
    
    if (sum(is.na(x))>0)
      warning("The matrix contain NA's. The corresponding values have been set to 0")
   
    min <- min(x)
    max <- max(x)
    yLabels <- rownames(x)
    xLabels <- colnames(x)
    title <-c()
    # check for additional function arguments
    if( length(list(...)) ){
      Lst <- list(...)
      if( !is.null(Lst$zlim) ){
        min <- Lst$zlim[1]
        max <- Lst$zlim[2]
      }
      if( !is.null(Lst$yLabels) ){
        yLabels <- c(Lst$yLabels)
      }
      if( !is.null(Lst$xLabels) ){
        xLabels <- c(Lst$xLabels)
      }
      if( !is.null(Lst$title) ){
        title <- Lst$title
      }
    }
    # check for null values
    if( is.null(xLabels) ){
      xLabels <- c(1:ncol(x))
    }
    if( is.null(yLabels) ){
      yLabels <- c(1:nrow(x))
    }
    
    #layout(matrix(data=c(1,2), nrow=2, ncol=1), widths=c(1,1), heights=c(1,4.5))
    layout(matrix(data=c(1,2,3,4), nrow=2, ncol=2), widths=c(5,1), heights=c(1,4.5))
    
    ColorRamp<- gray.colors(40)[40:0]
    ColorLevels <- seq(0, 1, length=length(ColorRamp))
    
    # Reverse Y axis
    reverse <- nrow(x) : 1
    yLabels <- yLabels[reverse]
    x <- x[reverse,]
    
    # Inclusion probs:
    par(mar = c(3,5,2.5,2))
    image(1:length(scale),1,
          matrix(data=scale, ncol=1, nrow=length(scale)),
          xlab="",ylab="",
          col=ColorRamp,
          yaxt="n",xaxt="n")
        
     if( !is.null(title) ){
          title(main=title)
        }# Title on the top
    
    # Data Map
    par(mar = c(3,5,2.5,2))
    image(1:length(xLabels), 1:length(yLabels), t(x), col=ColorRamp, xlab="",
          ylab="", axes=FALSE, zlim=c(min,max))

    axis(side=3, at=1:length(xLabels), labels=xLabels, cex.axis=0.7, las=2)
    axis(side=2, at=1:length(yLabels), labels=yLabels, las=1,
         cex.axis=0.7)
    
    #Nothing
    par(mar=c(3,5,2.5,2))
    plot(1,1,type="n", axes=0, xaxt="n", yaxt="n",xlab="",ylab="")
    
    # color scale:
    par(mar = c(3,5,2.5,2))
    par(mar = c(3,2.5,2.5,2))
    image(1, ColorLevels,
          matrix(data=ColorLevels, ncol=length(ColorLevels),nrow=1),
          col=ColorRamp,
          xlab="",ylab="",
          xaxt="n")
    }
 
 #Joint posterior probabilities
  if(auxtp=="j"){
    if(x$method=="gibbs"){
      myImagePlot(x$jointinclprob, scale=as.vector(x$inclprob[,1]), title="Estimated Joint Inclusion Probabilities")
    }else{
      myImagePlot(x$jointinclprob, scale=as.vector(x$inclprob[,1]), title="Joint Inclusion Probabilities")
    }
    
    
    }

 #conditional posterior probabilities
  if(auxtp=="c"){
    prob_joint <- as.matrix(x$jointinclprob)
    prob_cond <- as.matrix(x$jointinclprob)
    for(i in 1:(p-1))
      prob_cond[i,]<-prob_joint[i,]/prob_joint[i,i]
    
    if(x$method=="gibbs"){
      myImagePlot(prob_cond, scale=as.vector(x$inclprob[,1]), title="Estimated Inclusion prob. of column var. given the row var. is included")
    }else{
      myImagePlot(prob_cond, scale=as.vector(x$inclprob[,1]), title="Inclusion prob. of column var. given the row var. is included")
    }
    
    }
  #conditional posterior probabilities given Not a variable  
  if(auxtp=="n"){
      #auxiliar function to compute de A Given not B matrix
      pAgivenNotB<- function(miobject){
        result<- 0*miobject$jointinclprob
        for (i in 1:dim(result)[1]){
          for (j in 1:dim(result)[1]){
            result[i,j]<-
              (miobject$inclprob[j,1]-miobject$jointinclprob[i,j])/(1-miobject$inclprob[i,1])#REVISAR  P(j|Not.i)=(1-P(i|j))P(j)/(1-P(i))=(P(j)-P(j,i))/(1-P(i))  
          }}
        colnames(result)<- colnames(miobject$jointinclprob)
        rownames(result)<- paste("Not.",colnames(miobject$jointinclprob),sep="")
        result
      }
      AgivenNotB <- pAgivenNotB(x)
     
      if(x$method=="gibbs"){
        myImagePlot(AgivenNotB, scale=as.vector(x$inclprob[,1]), title="Est. Incl. probability of column var. given the row var. is NOT included")
      }else{
        myImagePlot(AgivenNotB, scale=as.vector(x$inclprob[,1]), title="Incl. probability of column var. given the row var. is NOT included")
      }
     
  }
}
