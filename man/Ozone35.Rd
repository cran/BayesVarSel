\name{Ozone35}
\alias{Ozone35}
\docType{data}
\title{
Ozone35 dataset
}
\description{
Polution data
}
\usage{data(Ozone35)}
\format{
  A data frame with 178 observations on the following 36 variables.
  \describe{
    \item{\code{y}}{Response = Daily maximum 1-hour-average ozone reading (ppm) at Upland, CA}
    \item{\code{x4}}{500-millibar pressure height (m) measured at Vandenberg AFB}
    \item{\code{x5}}{Wind speed (mph) at Los Angeles International Airport (LAX)}
    \item{\code{x6}}{Humidity (percentage) at LAX}
    \item{\code{x7}}{Temperature (Fahrenheit degrees) measured at Sandburg, CA}
    \item{\code{x8}}{Inversion base height (feet) at LAX}
    \item{\code{x9}}{Pressure gradient (mm Hg) from LAX to Daggett, CA}
    \item{\code{x10}}{Visibility (miles) measured at LAX}
        \item{\code{x4.x4}}{=x4*x4}
    \item{\code{x4.x5}}{=x4*x5}
    \item{\code{x4.x6}}{=x4*x6}
    \item{\code{x4.x7}}{=x4*x7}
    \item{\code{x4.x8}}{=x4*x8}
    \item{\code{x4.x9}}{=x4*x9}
    \item{\code{x4.x10}}{=x4*x10}
    \item{\code{x5.x5}}{=x5*x5}
    \item{\code{x5.x6}}{=x5*x6}
    \item{\code{x5.x7}}{=x5*x7}
    \item{\code{x5.x8}}{=x5*x8}
    \item{\code{x5.x9}}{=x5*x9}
    \item{\code{x5.x10}}{=x5*x10}
    \item{\code{x6.x6}}{=x6*x6}
    \item{\code{x6.x7}}{=x6*x7}
    \item{\code{x6.x8}}{=x6*x8}
    \item{\code{x6.x9}}{=x6*x9}
    \item{\code{x6.x10}}{=x6*x10}
    \item{\code{x7.x7}}{=x7*x7}
    \item{\code{x7.x8}}{=x7*x8}
    \item{\code{x7.x9}}{=x7*x9}
    \item{\code{x7.x10}}{=x7*x10}
    \item{\code{x8.x8}}{=x8*x8}
    \item{\code{x8.x9}}{=x8*x9}
    \item{\code{x8.x10}}{=x8*x10}
    \item{\code{x9.x9}}{=x9*x9}
    \item{\code{x9.x10}}{=x9*x10}
    \item{\code{x10.x10}}{=x10*x10}
  }
}
\details{
This dataset has been used by Garcia-Donato and Martinez-Beneito (2013) to illustrate the potential of the 
Gibbs sampling method (in \code{BayesVarSel} implemented in \code{\link[BayesVarSel]{GibbsBvs}}).
	
This data were previously used by Casella and Moreno (2006) and Berger and Molina (2005) and concern N = 178 measures of ozone concentration in the atmosphere. Of the 10 main effects originally considered, we only make use of those with an atmospheric meaning x4 to x10, as was done by Liang et al. (2008). We then have 7 main effects which, jointly with the quadratic terms and second order interactions, produce the above-mentioned p = 35 possible regressors.
}
\references{
	Berger, J. and Molina, G. (2005)<DOI:j.1467-9574.2005.00275.x> Posterior model probabilities via path-based pairwise priors. Statistica Neerlandica, 59:3-15.
	
	Casella, G. and Moreno, E. (2006)<DOI:10.1198/016214505000000646> Objective Bayesian variable selection. Journal of the American Statistical Association, 101(473).
	
	Garcia-Donato, G. and Martinez-Beneito, M.A. (2013)<DOI:10.1080/01621459.2012.742443> On sampling strategies in Bayesian variable selection problems with large model spaces. Journal of the American Statistical Association, 108: 340-352. 
	
	 Liang, F., Paulo, R., Molina, G., Clyde, M. and  Berger, J.O. (2008)<DOI:10.1198/016214507000001337> Mixtures of  g-priors for Bayesian Variable Selection. Journal of the American Statistical Association. 103:410-423.  
	 }
\examples{
data(Ozone35)
}
\keyword{datasets}
