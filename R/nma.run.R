#' Run NMA model
#' @description Takes bugs code from an object produced by \code{nma.model} and runs model using \code{jags}.
#' 
#' @param model A \code{BUGSnetModel} object produced by running \code{nma.model}.
#' @param monitor A vector of all variables that you would like to monitor. Default is "DEFAULT" which will monitor the relative treatment effects \code{d} 
#' as well as \code{sigma} when a random effects model is used and the regression coefficients \code{beta} when meta-regression is used.
#' @param DIC Default is TRUE and nodes required to calculate the DIC and other fit statistics are monitored. Otherwise you may set it to FALSE. 
#' @param n.adapt Number of adaptations for the mcmc chains.
#' @param n.burnin Number of burnin iterations for the mcmc chains.
#' @param n.iter Number of iterations for the mcmc chains.
#' @param thin Thinning factor for the mcmc chains. Default is 1.
#' @param n.chains Number of mcmc chains. Default is 3.
#' @param inits Specifies initial values and random number generator options for each chain. The "DEFAULT" option uses the R random seed to set the JAGS random
#' seeds. Non-default options are passed directly to \code{\link{jags.model}}. In order to use the JAGS default initialization, set inits to NULL. 
#' See \code{\link{jags.model}} for more info.
#' 
#' @return \code{nma.run} returns an object of class \code{BUGSnetRun} which is a list containing the following components:
#' @return \code{samples} - The MCMC samples produced by running the BUGS model.
#' @return \code{model} - The \code{BUGSnetModel} object obtained from \code{nma.model} which was used to run \code{jags}.
#' @return \code{scale} - The scale of the outcome, based on the chosen family and link function.
#' @return \code{trt.key} - Treatments mapped to numbers, used to run BUGS code.
#' @return \code{family} - Family that was used for the NMA model (e.g normal, binomial, poisson)
#' @return \code{link} - Link function that was used for the NMA model (e.g normal, binomial, poisson)
#' @export
#' @seealso \code{\link{nma.model}}, \code{\link{nma.fit}}, \code{\link{nma.league}}, \code{\link{nma.rank}}, \code{\link{nma.forest}}, \code{\link{nma.regplot}}, \code{\link{nma.trace}}, \code{\link{jags.model}}
#' @examples
#' data(diabetes.sim)
#' 
#' diabetes.slr <- data.prep(arm.data = diabetes.sim, 
#' varname.t = "Treatment", 
#' varname.s = "Study")
#' 
#' #Random effects, consistency model.
#' #Binomial family, cloglog link. This implies that the scale will be the Hazard Ratio.
#'diabetes.re.c <- nma.model(data = diabetes.slr,
#'        outcome = "diabetes", 
#'        N = "n",
#'        reference = "Placebo",
#'        family = "binomial",
#'        link = "cloglog",
#'        effects = "random",
#'        type="consistency",
#'        time="followup"
#'        )
#'  
#'diabetes.re.c.res <- nma.run(diabetes.re.c,
#'n.adapt=100,
#'n.burnin=0,
#'n.iter=100)
#'        


nma.run <- function(model,
                    monitor="DEFAULT",
                    DIC=TRUE,
                    n.adapt = 1000, 
                    n.burnin = floor(n.iter / 2), 
                    n.iter, 
                    thin=1,
                    n.chains=3,
                    inits = "DEFAULT"){
  
  if(class(model) != "BUGSnetModel")
    stop("\'model\' must be a valid BUGSnetModel object created using the nma.model function.")
  
  if (!is.null(inits) && inits == "DEFAULT")
  {
    seeds <- sample(.Machine$integer.max, n.chains, replace = FALSE)
    inits <- list()
    for (i in 1:n.chains)
      inits[[i]] <- list(.RNG.seed = seeds[i], .RNG.name = "base::Mersenne-Twister")
  }
  
  jagsmodel <- jags.model(textConnection(model$bugs),        #Create a connection so JAGS can access the variables
                          model$data,
                          n.chains=n.chains,
                          n.adapt=n.adapt,
                          inits = inits)
  
  if(n.burnin!=0) jagsburnin <- update(jagsmodel, n.iter=n.burnin)
  
  if(length(monitor)==1){
    if(monitor=="DEFAULT"){
      make.monitor <- "d"
      if(model$effects=="random") make.monitor <- c(make.monitor, "sigma")
      if(!is.null(model$covariate)) make.monitor <- c(make.monitor, "beta")
    }else make.monitor <- unique(c(monitor, "d"))
  }else make.monitor <- unique(c(monitor, "d"))
  
  if(DIC==TRUE){
    if (model$family == "binomial"){
      DIC.monitor <- c("dev", "r", "n","totresdev","rhat")
    } else if(model$family == "poisson"){
      DIC.monitor <- c("dev", "r","totresdev","theta")
    } else if(model$family == "normal"){
      DIC.monitor <- c("theta", "prec", "y", "totresdev", "dev")
    }
    new.monitor <- unique(c(make.monitor, DIC.monitor))
  } else if(DIC==FALSE){
    new.monitor <- make.monitor
  }
  
  jagssamples <- coda.samples(jagsmodel, 
                              variable.names=new.monitor,
                              n.iter=n.iter,
                              thin=thin)
  
  
  # print("The baseline treatment was ...")
  
  brun <- structure(list(samples=jagssamples,
                         model=model, 
                         scale=model$scale,
                         family=model$family,
                         link =model$link,
                         "trt.key"=as.character(t(model$trt.key[1]))),
                    class = "BUGSnetRun")
  return(brun)
  
}