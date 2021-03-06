% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Btest.R
\name{print.Btest}
\alias{print.Btest}
\title{Print an object of class \code{Btest}}
\usage{
\method{print}{Btest}(x, ...)
}
\arguments{
\item{x}{Object of class Btest}

\item{...}{Additional parameters to be passed}
}
\description{
Print an object of class \code{Btest}
}
\examples{
 \dontrun{
#Analysis of Crime Data
#load data
data(UScrime)
#Model selection among the following models: (note model1 is nested in all the others)
model1<- y ~ 1 + Prob
model2<- y ~ 1 + Prob + Time
model3<- y ~ 1 + Prob + Po1 + Po2
model4<- y ~ 1 + Prob + So
model5<- y ~ .

#Equal prior probabilities for models:
crime.BF<- Btest(models=list(basemodel=model1,
	ProbTimemodel=model2, ProbPolmodel=model3,
	ProbSomodel=model4, fullmodel=model5), data=UScrime)
	crime.BF
	}
}
\seealso{
See \code{\link[BayesVarSel]{Btest}} for creating objects of the class \code{Btest}.
}
