\name{Bvs}
\alias{Bvs}
\alias{print.Bvs}
\alias{summary.Bvs}
\alias{print}
\alias{summary}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Bayesian Variable Selection for linear regression models
}
\description{
Exact computation of summaries of the posterior distribution using sequential computation.
}
\usage{
Bvs(formula, data, prior.betas = "Robust", 
    prior.models = "Constant", n.keep, time.test = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{formula }{Formula defining the most complex regression model in the analysis (package forces the intercept always to be included). See details.}
  \item{data }{data frame containing the data.}
  \item{prior.betas }{Prior distribution for regression parameters within each model. Possible choices include "Robust", "Liangetal", "gZellner" and "ZellnerSiow" (see details).}
  \item{prior.models }{Prior distribution over the model space. Possible choices are "Constant" and "ScottBerger" (see details).}
  \item{n.keep }{How many of the most probable models are to be kept?}
  \item{time.test }{If TRUE and the number of variables is moderately large (>=18) a preliminary test to estimate computational time is performed.}
}
\details{The model space is the set of all models, Mi, that contain the intercept and are nested in that specified by \code{formula}. The simplest of such models, M0, contains only the intercept. Then \code{Bvs} provides exact summaries of the posterior distribution over this model space, that is, summaries of the 
discrete distribution which assigns to each model Mi its probability given the data: 

Pr(Mi | \code{data})=Pr(Mi)*Bi/C, 

where Bi is the Bayes factor of Mi to M0, Pr(Mi) is the prior probability of Mi and C is the normalizing constant. 

The Bayes factor B_i depends on the prior assigned for the regression parameters in Mi and \code{Bvs} implements a number of popular choices plus the "Robust" prior recently proposed by Bayarri et al (2012). The "Robust" prior is the default choice for both theoretical (see the reference for details) and computational reasons since it produces Bayes factors with closed-form expressions. The "gZellner" prior implemented corresponds to the prior in Zellner (1986) with g=n while the "Liangetal" prior is the hyper-g/n with a=3 (see the original paper Liang et al 2008, for details). Finally, "ZellnerSiow" is the multivariate Cauchy prior proposed by Zellner and Siow (1980, 1984), further studied by Bayarri and Garcia-Donato (2007).

With respect to the prior over the model space Pr(Mi) two possibilities are implemented: "Constant", under which every model has the same prior probability and "ScottBerger" under which Pr(Mi)=1/((p-1)*choose(p-1, ki-1)) (assuming p and ki are the number of explanatory variables in the most complex model and Mi respectively). The "ScottBerger" prior was proposed by Scott and Berger (2010) and 
controls for the potential pernicious effect on the posterior probabilities of a too large complex model.

This function works for problems of up to p=31. 

Warnings and limitations: the current version of the package does not allow for factors to be used as explanatory variables and results obtained are essentially unpredictable. The Bayes factors can be extremely big numbers if the sample size is large, potentially leading to NA's.
}


\value{\code{Bvs} returns an object of class \code{Bvs} with the following elements:
\item{time }{The internal time consumed in solving the problem}
\item{lm }{The \code{lm} class object that results when the model defined by \code{formula} is fitted by \code{lm}}
\item{variables }{The name of all the potential explanatory variables.}
\item{n }{Number of observations}
\item{p }{Total number of explanatory variables (including the intercept) in the most complex model}
\item{HPMbin }{The binary expression of the Highest Posterior Probability model}
\item{modelsprob }{A \code{data.frame} which summaries the \code{n.keep} most probable, a posteriori models, and their associated probability.}
\item{inclprob }{A \code{data.frame} with the inclusion probabilities of all the variables.}
\item{jointinclprob }{A \code{data.frame} with the joint inclusion probabilities of all the variables.}
\item{postprobdim }{Posterior probabilities of the dimension of the true model}
\item{betahat }{The model-averaged estimator of the regression parameters.}
\item{call }{The \code{call} to the function}
\item{method }{\code{full}} 
}

\author{
Gonzalo Garcia-Donato and Anabel Forte

Maintainer: <anabel.forte@uji.es>
}
\references{
  Bayarri, M.J., Berger, J.O., Forte, A. and Garcia-Donato, G. (2012) Criteria for Bayesian Model choice with Application to Variable Selection. The Annals of Statistics. 40: 1550-1557.

  Bayarri, M.J. and Garcia-Donato, G. (2007). Extending conventional priors for testing general hypotheses in linear models. Biometrika, 94:135-152.
  
  Barbieri, M and Berger, J (2004). Optimal Predictive Model Selection. The Annals of Statistics, 32, 870-897.

  Liang, F., Paulo, R., Molina, G., Clyde, M. and  Berger,
  J.O. (2008) Mixtures of  g-priors for Bayesian Variable
  Selection. Journal of the American Statistical Association. 103:410-423.%  \cr \url{http://dx.doi.org/10.1198/016214507000001337}
  
  Zellner, A. and Siow, A. (1980). Posterior Odds Ratio for Selected Regression Hypotheses. In Bayesian Statistics 1 (J.M. Bernardo, M. H. DeGroot, D. V. Lindley and A. F. M. Smith, eds.) 585-603. Valencia: University Press. 
  
  Zellner, A. and Siow, A. (1984). Basic Issues in Econometrics. Chicago: University of
Chicago Press.
 
  Zellner, A. (1986). On Assessing Prior Distributions and Bayesian Regression Analysis with g-prior Distributions. In Bayesian Inference and Decision techniques: Essays in Honor of Bruno de Finetti (A. Zellner, ed.) 389-399. Edward Elgar Publishing
Limited. 

}

\keyword{ package }
\seealso{
\code{\link[BayesVarSel]{plotBvs}} for several plots of the result.

\code{\link[BayesVarSel]{PBvs}} for a parallelized version of \code{Bvs}.

\code{\link[BayesVarSel]{GibbsBvs}} for a heuristic approximation based on Gibbs sampling (recommended when p>20, no other possibilities when p>31).
}

\examples{
\dontrun{
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

#Two image plots of the conditional inclusion probabilities:
plotBvs(crime.Bvs, option="conditional")
plotBvs(crime.Bvs, option="not")
}
}

