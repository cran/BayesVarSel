\name{predict.Bvs}
\alias{predict.Bvs}
\alias{print.Bvs.predict}
\alias{predict}
\title{
Predict method for object of the class \code{Bvs}
}
\description{
Optimal predicted values based on \code{Bvs} object.
}
\usage{
\method{predict}{Bvs}(object, newdata, ...)
}
\arguments{
  \item{object}{An object of class \code{Bvs}
 }
  \item{newdata}{An optional data frame in which to look for variables with
          which to predict. If omitted, the fitted values are used.}
  \item{\dots}{Not used}
}
\details{
The prediction is the scalar product of each new data point times the model averaged estimate of the regression parameters.
}
\value{An object of the class \code{Bvs.predict} which is a list with the following components:
	\describe{
	\item{X}{A data frame with the design matrix used in the prediction}
	\item{prediction }{The predictions for each of the rows in \code{X}}
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
pred.mean<- predict(object=hald.Bvs, newdata=as.data.frame(t(colMeans(Hald))))

#prediction at the original design matrix
predict(object=hald.Bvs)
}
