#' Bayes factors and posterior probabilities for linear regression models
#'
#' It Computes the Bayes factors and posterior probabilities of a list of linear
#' regression models proposed to explain a common response variable over the
#' same dataset
#'
#' The Bayes factors, BFi0, are expressed in relation with the simplest model
#' (the one nested in all the others). Then, the posterior probabilities of the
#' entertained models are obtained as
#'
#' Pr(Mi | \code{data})=Pr(Mi)*BFi0/C,
#'
#' where Pr(Mi) is the prior probability of model Mi and C is the normalizing
#' constant.
#'
#' The Bayes factor B_i depends on the prior assigned for the parameters in the regression
#' models Mi and \code{Bvs} implements a number of popular choices. The "Robust"
#' prior by Bayarri, Berger, Forte and Garcia-Donato (2012) is the recommended (and default) choice.
#' This prior prior can be implemented in a more stable way using the derivations in Greenaway (2019) 
#' and that are available in BayesVarSel since version 2.2.x setting the argument to "Robust.G". 
#'
#' Additional options are "gZellner" a prior which 
#' corresponds to the
#' prior in Zellner (1986) with g=n. Also "Liangetal" prior is the
#' hyper-g/n with a=3 (see the original paper Liang et al 2008, for details).
#' "ZellnerSiow" is the multivariate Cauchy prior proposed by Zellner and Siow
#' (1980, 1984), further studied by Bayarri and Garcia-Donato (2007).
#' "FLS" is the (benchmark) prior recommended by Fernandez, Ley and Steel (2001) which is
#' the prior in Zellner (1986) with g=max(n, p*p) p being the number of
#' covariates to choose from (the most complex model has p+number of fixed
#' covariates).
#' "intrinsic.MGC" is the intrinsic prior derived by Moreno, Giron, Casella (2015)
#' and "IHG" corresponds to the intrinsic hyper-g prior derived in Berger, Garcia-Donato, 
#' Moreno and Pericchi (2022).
#'
#' With respect to the prior over the model space Pr(Mi) three possibilities
#' are implemented: "Constant", under which every model has the same prior
#' probability and "User". With this last option, the prior probabilities are
#' defined through the named list \code{priorprobs}. These probabilities can be
#' given unnormalized.
#'
#' Limitations: the error "A Bayes factor is infinite.". Bayes factors can be
#' extremely big numbers if i) the sample size is even moderately large or if
#' ii) a model is much better (in terms of fit) than the model taken as the
#' null model. We are currently working on more robust implementations of the
#' functions to handle these problems. In the meanwhile you could try using the
#' g-Zellner prior (which is the most simple one and results, in these cases,
#' should not vary much with the prior) and/or using more accurate definitions
#' of the simplest model.
#'
#' @aliases Btest
#' @export
#' @param models A named list with the entertained models defined with their
#' corresponding formulas. If the list is unnamed, default names are given by
#' the routine. One model must be nested in all the others.
#' @param data data frame containing the data.
#' @param prior.betas Prior distribution for regression parameters within each
#' model (to be literally specified). Possible choices include "Robust", "Robust.G", "Liangetal", "gZellner",
#' "ZellnerSiow", "FLS", "intrinsic.MGC" and "IHG" (see details).
#' @param prior.models Type of prior probabilities of the models (to be literally specified). Possible choices are
#' "Constant" and "User" (see details).
#' @param priorprobs A named vector ir list (same length and names as in argument
#' \code{models}) with the prior probabilities of the models (used in combination
#' of \code{prior.models="User"}). If the provided object
#' is not named, then the order in the list of \code{models} is used to assign the prior
#' probabilities
#' @param null.model The name of the null model (eg. the one nested in all the others).
#' By default, the names of covariates in the different
#' models are used to identify the null model. An error is produced if such identification fails.
#' This identification
#' is not performed if the definition of the null model is provided, with this argument,
#' by the user.
#' Note that the (the \code{null.model} must coincide with that model with the largest
#' sum of squared errors and should be smaller in dimension to any other model).
#' @return \code{Btest} returns an object of type \code{Btest} which is a
#' \code{list} with the following elements: \item{BFio }{A vector with the
#' Bayes factor of each model to the simplest model.} \item{PostProbi }{A
#' vector with the posterior probabilities of each model.} \item{models }{A
#' list with the entertained models. } \item{nullmodel}{The position of the
#' simplest model. }
#' \item{prior.betas}{prior.betas}
#' \item{prior.models}{prior.models}
#' \item{priorprobs}{priorprobs}

#' @author Gonzalo Garcia-Donato and Anabel Forte
#'
#' Maintainer: <anabel.forte@@uv.es>
#' @seealso \code{\link[BayesVarSel]{Bvs}} for variable selection within linear
#' regression models
#' @references Bayarri, M.J., Berger, J.O., Forte, A. and Garcia-Donato, G.
#' (2012)<DOI:10.1214/12-aos1013> Criteria for Bayesian Model choice with
#' Application to Variable Selection. The Annals of Statistics. 40: 1550-1557.
#'
#' Bayarri, M.J. and Garcia-Donato, G. (2007)<DOI:10.1093/biomet/asm014>
#' Extending conventional priors for testing general hypotheses in linear
#' models. Biometrika, 94:135-152.
#'
#' Barbieri, M and Berger, J (2004)<DOI:10.1214/009053604000000238> Optimal
#' Predictive Model Selection. The Annals of Statistics, 32, 870-897.
#'
#' Berger, J., Garcıa-Donato, G., Moreno, E., and Pericchi, L. (2022). 
#' The intrinsic hyper-g prior for normal linear models. in preparation.
#'
#' Fernandez, C., Ley, E. and Steel, M.F.J.
#' (2001)<DOI:10.1016/s0304-4076(00)00076-2> Benchmark priors for Bayesian
#' model averaging. Journal of Econometrics, 100, 381-427.
#'
#' Greenaway, M. (2019) Numerically stable approximate Bayesian methods for
#' generalized linear mixed models and linear model selection. Thesis (Department of Statistics,
#' University of Sydney).
#'
#' Liang, F., Paulo, R., Molina, G., Clyde, M. and Berger,J.O.
#' (2008)<DOI:10.1198/016214507000001337> Mixtures of g-priors for Bayesian
#' Variable Selection. Journal of the American Statistical Association.
#' 103:410-423
#'
#' Zellner, A. and Siow, A. (1980)<DOI:10.1007/bf02888369> Posterior Odds Ratio
#' for Selected Regression Hypotheses. In Bayesian Statistics 1 (J.M. Bernardo,
#' M. H. DeGroot, D. V. Lindley and A. F. M. Smith, eds.) 585-603. Valencia:
#' University Press.
#'
#' Zellner, A. and Siow, A. (1984) Basic Issues in Econometrics. Chicago:
#' University of Chicago Press.
#'
#' Zellner, A. (1986)<DOI:10.2307/2233941> On Assessing Prior Distributions and
#' Bayesian Regression Analysis with g-prior Distributions. In Bayesian
#' Inference and Decision techniques: Essays in Honor of Bruno de Finetti (A.
#' Zellner, ed.) 389-399. Edward Elgar Publishing Limited.
#' @keywords package
#' @examples
#'
#' \dontrun{
#' #Analysis of Crime Data
#' #load data
#' data(UScrime)
#' #Model selection among the following models: (note model1 is nested in all the others)
#' model1<- y ~ 1 + Prob
#' model2<- y ~ 1 + Prob + Time
#' model3<- y ~ 1 + Prob + Po1 + Po2
#' model4<- y ~ 1 + Prob + So
#' model5<- y ~ .
#'
#' #Equal prior probabilities for models:
#' crime.BF<- Btest(models=list(basemodel=model1,
#' 	ProbTimemodel=model2, ProbPolmodel=model3,
#' 	ProbSomodel=model4, fullmodel=model5), data=UScrime)
#'
#' #Another configuration of prior probabilities of models:
#' crime.BF2<- Btest(models=list(basemodel=model1, ProbTimemodel=model2,
#' 	ProbPolmodel=model3, ProbSomodel=model4, fullmodel=model5),
#' 	data=UScrime, prior.models = "User", priorprobs=list(basemodel=1/8,
#' 	ProbTimemodel=1/8, ProbPolmodel=1/2, ProbSomodel=1/8, fullmodel=1/8))
#' #same as:
#' #crime.BF2<- Btest(models=list(basemodel=model1, ProbTimemodel=model2,
#' 	#ProbPolmodel=model3,ProbSomodel=model4, #fullmodel=model5), data=UScrime,
#' 	#prior.models = "User", priorprobs=list(basemodel=1, ProbTimemodel=1,
#' 	#ProbPolmodel=4, #ProbSomodel=1, fullmodel=1))
#' }
Btest <-
  function(models,
           data,
           prior.betas = "Robust",
           prior.models = "Constant",
           priorprobs = NULL,
		       null.model = NULL) {

	#N is the number of models:
  N <- length(models)

	if (!is.list(models)) stop("Argument models should be a list\n")

	#If competing models come wihtout a name, give one by default:
	if (is.null(names(models))){
		if (!is.null(null.model)) stop("Please provide a name for the competing models. The null model must be in that list\n")
		names(models) <- paste("model", 1:N, sep="")
	}

	#Check if the given null model is one of the competing models:
	if (!is.null(null.model)){
	  relax.nest = TRUE
		pos.user.null.model <- which(null.model == names(models))
		if (length(pos.user.null.model) == 0) stop("The null model provided is not in the list of competing models.\n")
	}
  else relax.nest = FALSE


    #n is the sample size
    n <- dim(data)[1]
    #SSE is a vector with SSE's for each model; Dim with the dimension (number of regressors in each)
    SSE <- rep(0, N)
    Dim <- rep(0, N)
    BFi0 <- rep(0, N)
    PostProbi <- rep(0, N)

    #prior for betas: (#pfb will contain the internal representation of the prior.betas argument)
	pfb <- check.prior.betas(prior.betas) 
	
    #prior for model space:
    pfms <- substr(tolower(prior.models), 1, 1)
    if (pfms != "c" &&
        pfms != "u")
      stop("I am very sorry: prior for models not supported\n")
    if (pfms == "u" &&
        is.null(priorprobs)) {
      stop("A valid vector of prior probabilities must be provided\n")
    }
    if (pfms == "u" &&
        length(priorprobs) != N) {
      stop("Vector of prior probabilities with incorrect length\n")
    }
    if (pfms == "u" &&
        sum(priorprobs < 0) > 0) {
      stop("Prior probabilities must be positive\n")
    }
    #If priorprobs are unnamed, give this vector the name of the models
    if (pfms == "u" &&
        is.null(names(priorprobs))) {
          names(priorprobs)<- names(models)
    }

    #Prior probabilities of models:
    PriorModels <- rep(0, N)
    if (prior.models == "Constant") {
      PriorModels <- rep(1, N)
    }
    if (prior.models == "User") {
      #should coincide with the length of prior.models
      for (i in 1:N) {
        PriorModels[i] <- priorprobs[[names(models)[i]]]
      }
    }

    #list that contains the names of the covariates in each model
    covar.list <- list()
    for (i in 1:N) {
      temp <-
        lm(
          formula = as.formula(models[[i]]),
          data = data,
          y = TRUE,
          x = TRUE
        )
      SSE[i] <- sum(temp$residuals ^ 2)
      Dim[i] <- length(temp$coefficients)
      covar.list[[i]] <- dimnames(temp$x)[[2]]
    }
    ordered.SSE <- sort(SSE, index.return = TRUE, decreasing = TRUE)
    #Which acts as null model:
    nullmodel <- ordered.SSE$ix[1]

	#Check that the given null model and the model with largest SSE
	#(which in turns defines the null model)
	#coincides
	if (!is.null(null.model)){
		if (nullmodel != pos.user.null.model){
		stop("The given null model does not coincide with the one with the largest sum of squared error (and it should).\n")
		}
	}

    if (pfb != "f") {
      for (i in (1:N)[-nullmodel]) {
        #check if the "null" model is nested in all the others
        if (!relax.nest &
            sum(covar.list[[nullmodel]] %in% covar.list[[i]]) < Dim[nullmodel]) {
            stop("I suspect that perhaps the simplest (null) model is not nested in all the others. Define explicitly the simplest model if you are sure it is the case\n.\n")
        }
        Qi0 <- SSE[i] / SSE[nullmodel]
        #The .C to be used:
        if (pfb == "g")
          BFi0[i] <-
          .C(
            "gBF",
            as.integer(n),
            as.integer(Dim[i]),
            as.integer(Dim[nullmodel]),
            as.double(Qi0),
            as.double(0.0)
          )[5][[1]]
        if (pfb == "r")
          BFi0[i] <-
          .C(
            "RobustBF",
            as.integer(n),
            as.integer(Dim[i]),
            as.integer(Dim[nullmodel]),
            as.double(Qi0),
            as.double(0.0)
          )[5][[1]]
        if (pfb == "z")
          BFi0[i] <-
          .C(
            "ZSBF",
            as.integer(n),
            as.integer(Dim[i]),
            as.integer(Dim[nullmodel]),
            as.double(Qi0),
            as.double(0.0)
          )[5][[1]]
        if (pfb == "l")
          BFi0[i] <-
          .C(
            "LiangBF",
            as.integer(n),
            as.integer(Dim[i]),
            as.integer(Dim[nullmodel]),
            as.double(Qi0),
            as.double(0.0)
          )[5][[1]]
	      if (pfb == "i")
	        BFi0[i] <-
	        .C(
	          "intrinsicBF",
	          as.integer(n),
	          as.integer(Dim[i]),
	          as.integer(Dim[nullmodel]),
	          as.double(Qi0),
	          as.double(0.0)
	        )[5][[1]]
  	      if (pfb == "x")
  	        BFi0[i] <-
  	        .C(
  	          "geointrinsicBF",
  	          as.integer(n),
  	          as.integer(Dim[i]),
  	          as.integer(Dim[nullmodel]),
  	          as.double(Qi0),
  	          as.double(0.0)
  	        )[5][[1]]
    	  if (pfb == "y")
    	    BFi0[i] <-
    	    .C(
    	      "geointrinsic2BF",
    	      as.integer(n),
    	      as.integer(Dim[i]),
    	      as.integer(Dim[nullmodel]),
    	      as.double(Qi0),
    	      as.double(0.0)
    	    )[5][[1]]
        }
    }

    if (pfb == "f") {
      p <- max(Dim) - min(Dim)
      for (i in (1:N)[-nullmodel]) {
        #check if the "null" model is nested in all the others
        if (!relax.nest &
            sum(covar.list[[nullmodel]] %in% covar.list[[i]]) < Dim[nullmodel]) {
          stop("There is no a model nested in all the others\n")
        }
        Qi0 <- SSE[i] / SSE[nullmodel]
        BFi0[i] <-
          .C(
            "flsBF",
            as.integer(p),
            as.integer(n),
            as.integer(Dim[i]),
            as.integer(Dim[nullmodel]),
            as.double(Qi0),
            as.double(0.0)
          )[6][[1]]
      }
    }


    BFi0[nullmodel] <- 1
    names(BFi0) <-
      paste(names(models), ".to.", names(models)[nullmodel], sep = "")
    PostProbi <- BFi0 * PriorModels / sum(BFi0 * PriorModels)
    names(PostProbi) <- names(models)
    result <- list()
    result$BFi0 <- BFi0
    result$PostProbi <- PostProbi
    result$models <- models
    result$nullmodel <- nullmodel
	
	result$prior.betas<- prior.betas	
	result$prior.models<- prior.models
	result$priorprobs<- priorprobs
	
    class(result) <- "Btest"
    result
  }



#' Print an object of class \code{Btest}
#'
#'Print an object of class \code{Btest}
#' @export
#' @param x Object of class Btest
#' @param ... Additional parameters to be passed
#'
#' @seealso See \code{\link[BayesVarSel]{Btest}} for creating objects of the class \code{Btest}.
#'
#'
#' @examples
#'  \dontrun{
#' #Analysis of Crime Data
#' #load data
#' data(UScrime)
#' #Model selection among the following models: (note model1 is nested in all the others)
#' model1<- y ~ 1 + Prob
#' model2<- y ~ 1 + Prob + Time
#' model3<- y ~ 1 + Prob + Po1 + Po2
#' model4<- y ~ 1 + Prob + So
#' model5<- y ~ .
#'
#' #Equal prior probabilities for models:
#' crime.BF<- Btest(models=list(basemodel=model1,
#' 	ProbTimemodel=model2, ProbPolmodel=model3,
#' 	ProbSomodel=model4, fullmodel=model5), data=UScrime)
#' 	crime.BF
#' 	}
print.Btest <- function(x, ...){
  cat("---------\n")
  cat("Models:\n")
  print(x$models)
  cat("---------\n")
  cat(paste("Bayes factors (expressed in relation to ",names(x$models)[x$nullmodel],")\n", sep=""))
  print(x$BFi0)
  cat("---------\n")
  cat("Posterior probabilities:\n")
  print(round(x$PostProbi,3))
}
