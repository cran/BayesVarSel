\name{predict.Bvs}
\alias{predict.Bvs}
\alias{print.Bvs.predict}
\alias{predict}
\title{
Predict method por object of the class \code{Bvs}
}
\description{
For a given new set of variables returns the predicted value as the sum of betahat*xnew where
betahat is the expected posterior vector of the regression parameters. 
}
\usage{
\method{predict}{Bvs}(object, newdata, ...)
}
\arguments{
  \item{object}{An object of class \code{Bvs}
 }
  \item{newdata}{Either a \code{data.frame}, a numeric vector or a \code{matrix} with the new values of the
  	variables for which the prediction is performed. This \code{newdata} should not contain the 1 corresponding to
  	the intercept.}
  \item{\dots}{Not used}
}
\value{An object of the class \code{Bvs.predict} which is a list with the following components:
	\describe{
	\item{X}{The matrix given for prediction.}
	\item{prediction }{The prediction values for each of the rows in \code{X}}
	\item{nn }{The number of design points being predicted}
	\item{call }{The call to the function}
	}
}
\author{
Gonzalo Garcia-Donato and Anabel Forte

Maintainer: <anabel.forte@uji.es>
}
\seealso{
\code{\link[BayesVarSel]{Bvs}}, 
  	\code{\link[BayesVarSel]{PBvs}} or \code{\link[BayesVarSel]{GibbsBvs}}
}
\examples{
#read Hald data
data(Hald)
#run the main function:
hald.Bvs<- Bvs(formula="y~x1+x2+x3+x4", data=Hald, n.keep=16)

#Prediction at the mean value:
pred.mean<- predict(object=hald.Bvs, newdata=colMeans(Hald[,c("x1","x2","x3","x4")]))

#prediction over the same design matrix
predict(object=hald.Bvs, newdata=Hald[,c("x1","x2","x3","x4")])
}
