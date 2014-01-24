\name{BayesVarSel-package}
\alias{BayesVarSel-package}
\alias{BayesVarSel}
\docType{package}
\title{
Bayesian Variable selection in Linear Models
}
\description{
This package provides specific tools for the analysis of the variable selection problem in linear regression models from a Bayesian perspective. It provides simple and intuitive methods to explore and synthesize the results and allows the calculations to be performed either exactly (sequential or parallel computation) or heuristically, using a Gibbs sampling algorithm studied in Garcia-Donato and Martinez-Beneito (2013). 

The default implementation takes advantage of a closed-form expression for the posterior probabilities that the "Robust" prior in Bayarri et al (2012) produces. Also, other  priors like Zellner (1986) g-prior, Zellner-Siow (1980,1984) or Liang et al (2008) prior can be used. See \code{\link[BayesVarSel]{Bvs}} for a more precise definition of the priors implemented.
}
\details{
\tabular{ll}{
Package: \tab BayesVarSel\cr
Type: \tab Package\cr
Version: \tab 1.4\cr
Date: \tab 2012-11-12\cr
License: \tab GPL-2\cr
}
}
\author{
Gonzalo Garcia-Donato and Anabel Forte

Maintainer: Anabel Forte \email{forte@uji.es}
}
\references{

  Bayarri, M.J., Berger, J.O., Forte, A. and Garcia-Donato, G. (2012) Criteria for Bayesian Model choice with Application to Variable Selection. The Annals of Statistics. 40: 1550-1577

   Garcia-Donato, G. and Martinez-Beneito, M.A. (2013) On sampling strategies in Bayesian variable selection problems with large model spaces. Journal of the American Statistical Association. 108: 340-352.
  
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
\code{\link[BayesVarSel]{Bvs}},
\code{\link[BayesVarSel]{PBvs}},
\code{\link[BayesVarSel]{GibbsBvs}}
}
\examples{
demo(BayesVarSel.Hald)
}
