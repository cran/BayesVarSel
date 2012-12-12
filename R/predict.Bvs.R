predict.Bvs <-
function (object, newdata, ...){
  pp<-object$p
  predictor <- as.matrix(object$betahat)
  prediction=list()
  prediction$predictor <- predictor
 
  #por si se fuerza a predict.Bvs con un tipo de objeto que no toca.  
  if (!inherits(object, "Bvs")) 
    warning("calling predict.Bvs(<fake-Bvs-object>) ...")
  
  #Si el valor missing no está
  if (missing(newdata) || is.null(newdata)) {
   X <- model.matrix(object$lm)
   prediction$X<-as.data.frame(X)
  }else{
    if (!is.data.frame(newdata)&&!is.vector(newdata)&&!is.matrix(newdata))
    stop("Not valid type of object for newdata. newdata should be a data.frame, matrix or vector")
    
    if(is.vector(newdata)){
      if(length(newdata)!=(pp-1)){
        stop(paste("Invalid number of variables,newdata should contain observations for ",pp-1," variables",sep=""))
        }else{
          X <- as.numeric(matrix(c(1,newdata),ncol=pp,nrow=1))
          prediction$X<-as.data.frame(t(X))
          }
      }
        
    if(is.data.frame(newdata)||is.matrix(newdata)){
      if(dim(newdata)[2]!=(pp-1)){
        
        stop(paste("Invalid number of variables,newdata should contain observations for ",pp-1," variables",sep=""))
            
        }else{
          if(dim(newdata)[1]==1){
            X <- as.numeric(matrix(c(1,newdata),ncol=pp,nrow=1))
            prediction$X<-as.data.frame(t(X))
          }else{
            nn<-dim(newdata)[1]
            X<-matrix(1,ncol=pp,nrow=nn)
            for(i in 1:nn)
              X[i,2:pp]<-as.numeric(newdata[i,])
            prediction$X<-as.data.frame(X)
          }
          
          }
      }
    }
   
 
  names(prediction$X)<- object$variables
  prediction$prediction<-as.data.frame(X%*%predictor)
  names(prediction$prediction) <- "Predicted Value"
  prediction$nn<-dim(prediction$prediction)[1]
  prediction$call <- match.call()
  class(prediction)<-"Bvs.predict"
  prediction
}
