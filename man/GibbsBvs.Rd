\name{GibbsBvs}
\alias{GibbsBvs}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Bayesian Variable Selection for linear regression models using Gibbs sampling.
}
\description{
Approximate computation of summaries of the posterior distribution using a Gibbs sampling algorithm 
to explore the model space and frequency of "visits" to construct the estimates.
}
\usage{
GibbsBvs(formula, data, prior.betas = "Robust",
         prior.models = "Constant", n.iter, init.model = "Null",
         n.burnin = 50, time.test = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{formula }{Formula defining the most complex regression model in the analysis (package forces the intercept always to be included). See details.}
  \item{data }{data frame containing the data.}
  \item{prior.betas }{Prior distribution for regression parameters within each model. Possible choices include "Robust", "Liangetal", "gZellner" and "ZellnerSiow" (see details).}
  \item{prior.models }{Prior distribution over the model space. Possible choices are "Constant" and "ScottBerger"  (see details).}
  \item{n.iter }{The number of iterations (after the burn in process). All these models are kept and are used in the estimates.}
  \item{init.model }{The model at which the simulation process starts. Options include "null" (the model with just intercept), "full" (the model defined by \code{formula}) and a vector with p (the number of covariates in the model defined by \code{formula}) zeros and ones defining a model.}
  \item{n.burnin }{The number of iterations at the beginning of the MCMC that
          are thrown away.}
  \item{time.test }{If TRUE and the number of variables is large (>=21) a preliminary test to estimate computational time is performed.}
}
\details{
This is a heuristic approximation to the function \code{\link[BayesVarSel]{Bvs}} so the details there apply also here.

The algorithm implemented is a Gibbs sampling-based searching algorithm originally proposed by George and McCulloch (1997). Garcia-Donato and Martinez-Beneito (2013) have shown that this simple sampling strategy in combination with estimates based on frequency of visits (the one here implemented) provides very reliable results.
}

\value{\code{GibbsBvs} returns an object of class \code{Bvs} with the following elements:
\item{time }{The internal time consumed in solving the problem}
\item{lm }{The \code{lm} class object that results when the model defined by \code{formula} is fitted by \code{lm}}
\item{variables }{The name of all the potential explanatory variables}
\item{n }{Number of observations}
\item{p }{Total number of explanatory variables (including the intercept) in the most complex model.}
\item{HPMbin }{The binary expression of the most probable model found.}
\item{modelsprob }{\code{NULL} (posterior probabilities of single models are not well estimated using this method).}
\item{inclprob }{A \code{data.frame} with the estimates of the inclusion probabilities of all the variables.}
\item{jointinclprob }{A \code{data.frame} with the estimates of the joint inclusion probabilities of all the variables.}
\item{postprobdim }{Estimates of posterior probabilities of the dimension of the true model.}
\item{betahat }{Estimates of the model-averaged estimator of the regression parameters.}
\item{call }{The \code{call} to the function.}
\item{method }{\code{gibbs}} 

}

\references{
Garcia-Donato, G. and Martinez-Beneito, M.A. (2013) On sampling strategies in Bayesian variable selection problems with large model spaces. Journal of the American Statistical Association, 108: 340-352. 
  
George E. and McCulloch R. (1997) Approaches for Bayesian variable selection. Statistica Sinica, 7, 339:372.
}

\author{
Gonzalo Garcia-Donato and Anabel Forte
}

\keyword{ package }
\seealso{
\code{\link[BayesVarSel]{plotBvs}} for a plot of the object created.
}



%%\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
%%}
\examples{
\dontrun{
#Analysis of Ozone35 data

data(Ozone35)

#We use here the (Zellner) g-prior for
#regression parameters and constant prior
#over the model space
#In this Gibbs sampling scheme, we perform 10100 iterations,
#of which the first 100 are discharged (burnin)
#as initial model we use the null model (only intercept)
Oz35.GibbsBvs<- GibbsBvs(formula="y~.", data=Ozone35, prior.betas="gZellner",
prior.models="Constant", n.iter=10000, init.model="null", n.burnin=100, time.test = FALSE)

#Note: this is a heuristic approach and results are estimates
#of the exact answer.

#with the print we can see which is the most probable model
#among the visited
Oz35.GibbsBvs

#The estimation of inclusion probabilities and
#the model-averaged estimation of parameters:
summary(Oz35.GibbsBvs)

#Plots:
plotBvs(Oz35.GibbsBvs, option="conditional")
}
}
