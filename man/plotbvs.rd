\name{plotBvs}
\alias{plotBvs}
\title{
A function for plotting summaries of an object of class \code{Bvs}
}
\description{
Four different plots to summarize graphically the results in an object of class \code{Bvs}.
}
\usage{
plotBvs(x, option = "dimension")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{An object of class \code{Bvs}
  	
}
  \item{option}{One of "dimension", "joint", "conditional" or "not"
}
}
\details{
If \code{option}="dimension"  this function returns a barplot of the posterior distribution of the dimension of the true model. If \code{option}="joint" an image plot of the joint inclusion probabilities is returned. If \code{option}="conditional" an image plot of the conditional inclusion probabilities. These should be read as the probabilty that the variable in the column is part of the true model if the corresponding variables on the row is. Finally, if \code{option}="not" the image plot that is returned is that of the the probabilty that the variable in the column is part of the true model if the corresponding variables on the row is not.
}
\author{
Gonzalo Garcia-Donato and Anabel Forte

Maintainer: <anabel.forte@uv.es>
}
\seealso{
See \code{\link[BayesVarSel]{Bvs}}, \code{\link[BayesVarSel]{PBvs}} and \code{\link[BayesVarSel]{GibbsBvs}} for creating objects of the class \code{Bvs}.
}

\examples{

#Analysis of Crime Data
#load data
data(UScrime)

#Default arguments are Robust prior for the regression parameters
#and constant prior over the model space
#Here we keep the 1000 most probable models a posteriori:
crime.Bvs<- Bvs(formula="y~.", data=UScrime, n.keep=1000)

#A look at the results:
crime.Bvs

summary(crime.Bvs)

#A plot with the posterior probabilities of the dimension of the
#true model:
plotBvs(crime.Bvs, option="dimension")

#An image plot of the joint inclusion probabilities:
plotBvs(crime.Bvs, option="joint")
 
#Two image plots of the conditional inclusion probabilities:
plotBvs(crime.Bvs, option="conditional")
plotBvs(crime.Bvs, option="not")

}

