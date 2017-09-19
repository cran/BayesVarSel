\name{histBMA}
\alias{histBMA}
\title{
A function for histograms-like representations of objects of class \code{bma.coeffs}
}
\description{
The columns in \code{bma.coeffs} are simulations of the model averaged posterior distribution. This normally is a mixture of a discrete (at zero) and several continuous distributions. This plot provides a convenient graphical summary of such distributions.
}
\usage{
histBMA(x, covariate, n.breaks=100, text = TRUE, gray.0 = 0.6, gray.no0=0.8)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{An object of class \code{bma.coeffs}}
  
  \item{covariate}{The name of an explanatory variable whose accompanying coefficient is to be represented. This must be the name of one of the columns in \code{x}}
  
  \item{n.breaks}{The number of equally lentgh bars for the histogram}
  
  \item{text}{If set to TRUE the probability of the coefficient being zero is added in top of the bar at zero. Note: this probability is based on the models used in \code{bma.coeffs} (see details in that function)}
  
  \item{gray.0}{A numeric value between 0 and 1 that specifies the darkness, in a gray scale (0 is white and 1 is black) of the bar at zero}
  
   \item{gray.no0}{A numeric value between 0 and 1 that specifies the darkness, in a gray scale (0 is white and 1 is black) of the bars different from zero}

  
}

\details{This function produces a histogram but with the peculiarity that the zero values in the simulation are represented as bar centered at zero. The area of all the bars is one and of these, the area of the bar at zero (colored with \code{gray.0}) is, conditionally on the retained models (see details in \code{\link[BayesVarSel]{BMAcoeff}}), the probability of that coefficient be exactly zero. This number is included in the top of the zero bar if \code{text} is set to TRUE.
}
\author{
Gonzalo Garcia-Donato and Anabel Forte

Maintainer: <anabel.forte@uv.es>
}
\seealso{
See \code{\link[BayesVarSel]{BMAcoeff}}. Also see \code{\link[BayesVarSel]{Bvs}}, \code{\link[BayesVarSel]{PBvs}} and \code{\link[BayesVarSel]{GibbsBvs}} for creating objects of the class \code{BMAcoeff}.
}

\examples{

\dontrun{
	
#Analysis of Crime Data
#load data
data(UScrime)

crime.Bvs<- Bvs(formula="y~.", data=UScrime, n.keep=1000)
crime.Bvs.BMA<- BMAcoeff(crime.Bvs, n.sim=10000)
#the best 1000 models are used in the mixture

#Observe the bimodality of the coefficient associated with regressor M
histBMA(crime.Bvs.BMA, "M")

#Note 1:
#The value in top of the bar at zero (0.251 in this case) is the probability of beta_M is
#zero conditional on a model space containing the 1000 models used in the mixture. This value
#should be closed to the exact value
#1-crime.Bvs$inclprob["M"]
#which in this case is 0.2954968
#if n.keep above is close to 2^15

#Note 2:
#The BMA posterior distribution of beta_M has two modes approximately located at 0 and 10
#If we summarize this distribution using the mean
mean(crime.Bvs.BMA[ ,"M"])
#or median
median(crime.Bvs.BMA[ ,"M"])
#we obtain values around 7 (or 7.6) which do not represent this distribution.

#With the Gibbs algorithms:
data(Ozone35)

Oz35.GibbsBvs<- GibbsBvs(formula="y~.", data=Ozone35, prior.betas="gZellner",
prior.models="Constant", n.iter=10000, init.model="Full", n.burnin=100, 
time.test = FALSE)
Oz35.GibbsBvs.BMA<- BMAcoeff(Oz35.GibbsBvs, n.sim=10000)

histBMA(Oz35.GibbsBvs.BMA, "x6.x7")
#In this case (Gibbs sampling), the value in top of the bar at zero (0.366) 
#basically should coincide (if n.sim is large enough)
#with the estimated complement of the inclusion probability
#1-Oz35.GibbsBvs$inclprob["x6.x7"]
#which in this case is 0.3638
}
}
