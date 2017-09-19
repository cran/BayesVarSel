\name{BayesVarSel-package}
\alias{BayesVarSel-package}
\alias{BayesVarSel}
\docType{package}
\title{
Bayes Factors, Model Choice And Variable Selection In Linear Models
}
\description{
Hypothesis testing, model selection and model averaging are important statistical problems that have in common the explicit consideration of the uncertainty about which is the true model. The formal Bayesian tool to solve such problems is the Bayes factor (Kass and Raftery, 1995) that reports the evidence in the data favoring each of the entertained hypotheses/models and can be easily translated to posterior probabilities. 

This package has been specifically conceived to calculate Bayes factors in linear models and then to provide a formal Bayesian answer to testing and variable selection problems. From a theoretical side, the emphasis in the package is placed on the prior distributions (a very delicate issue in this context) and BayesVarSel allows using a wide range of them: Jeffreys-Zellner-Siow (Jeffreys, 1961; Zellner and Siow, 1980,1984) Zellner (1986); Fernandez et al. (2001), Liang et al. (2008) and Bayarri et al. (2012).

The interaction with the package is through a friendly interface that syntactically mimics the well-known lm command of R. The resulting objects can be easily explored providing the user very valuable information (like marginal, joint and conditional inclusion probabilities of potential variables; the highest posterior probability model, HPM; the median probability model, MPM) about the structure of the true -data generating- model. Additionally, BayesVarSel incorporates abilities to handle problems with a large number of potential explanatory variables through parallel and heuristic versions (Garcia-Donato and Martinez-Beneito 2013) of the main commands.
}
\details{
\tabular{ll}{
Package: \tab BayesVarSel\cr
Type: \tab Package\cr
Version: \tab 1.7.1\cr
Date: \tab 2017-09-19\cr
License: \tab GPL-2\cr
}
}
\author{
Gonzalo Garcia-Donato and Anabel Forte

Maintainer: Anabel Forte \email{anabel.forte@uv.es}
}
\references{

  Bayarri, M.J., Berger, J.O., Forte, A. and Garcia-Donato, G. (2012)<DOI:10.1214/12-aos1013> Criteria for Bayesian Model choice with Application to Variable Selection. The Annals of Statistics. 40: 1550-1577

  Fernandez, C., Ley, E. and Steel, M.F.J. (2001)<DOI:10.1016/s0304-4076(00)00076-2> Benchmark priors for Bayesian model averaging. Journal of Econometrics, 100, 381-427.  

   Garcia-Donato, G. and Martinez-Beneito, M.A. (2013)<DOI:10.1080/01621459.2012.742443> On sampling strategies in Bayesian variable selection problems with large model spaces. Journal of the American Statistical Association. 108: 340-352.
  
  Liang, F., Paulo, R., Molina, G., Clyde, M. and  Berger, J.O. (2008)<DOI:10.1198/016214507000001337> Mixtures of  g-priors for Bayesian Variable Selection. Journal of the American Statistical Association. 103:410-423.
  
  Zellner, A. and Siow, A. (1980)<DOI:10.1007/bf02888369>. Posterior Odds Ratio for Selected Regression Hypotheses. In Bayesian Statistics 1 (J.M. Bernardo, M. H. DeGroot, D. V. Lindley and A. F. M. Smith, eds.) 585-603. Valencia: University Press. 
  
  Zellner, A. and Siow, A. (1984) Basic Issues in Econometrics. Chicago: University of Chicago Press.
 
  Zellner, A. (1986)<DOI:10.2307/2233941> On Assessing Prior Distributions and Bayesian Regression Analysis with g-prior Distributions. In Bayesian Inference and Decision techniques: Essays in Honor of Bruno de Finetti (A. Zellner, ed.) 389-399. Edward Elgar Publishing
Limited. 

}
\keyword{ package }
\seealso{
\code{\link[BayesVarSel]{Btest}},
\code{\link[BayesVarSel]{Bvs}},
\code{\link[BayesVarSel]{PBvs}},
\code{\link[BayesVarSel]{GibbsBvs}},
\code{\link[BayesVarSel]{BMAcoeff}},
\code{\link[BayesVarSel]{predictBvs}}
}
\examples{
demo(BayesVarSel.Hald)
}
