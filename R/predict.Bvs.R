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
   prediction$X <- model.matrix(object$lm)
  }
  else{
    if (!is.data.frame(newdata)) stop("newdata should be a data.frame")
	
	tt<- terms(object$lm)
	Terms <- delete.response(tt)
	prediction$X<- model.matrix(Terms, model.frame(Terms, newdata))
  }  
    
  prediction$prediction<-as.data.frame(prediction$X%*%predictor)
  names(prediction$prediction) <- "Predicted Value"
  prediction$nn<-dim(prediction$prediction)[1]
  prediction$call <- match.call()
  prediction$X<- as.data.frame(prediction$X)
  class(prediction)<-"Bvs.predict"
  prediction
}
