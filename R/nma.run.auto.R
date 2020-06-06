#' Run NMA model
#' @description Takes bugs code from an object produced by \code{nma.model} and runs model using \code{jags}.
#' 
#' @param model A \code{BUGSnetModel} object produced by running \code{nma.model}.
#' @param monitor A vector of all variables that you would like to monitor. Default is "DEFAULT" which will monitor the relative treatment effects \code{d} 
#' as well as \code{sigma} when a random effects model is used and the regression coefficients \code{beta} when meta-regression is used.
#' @param DIC Default is TRUE and nodes required to calculate the DIC and other fit statistics are monitored. Otherwise you may set it to FALSE. 
#' @param inits Specifies initial values and random number generator options for each chain. The "DEFAULT" option uses the R random seed to set the JAGS random
#' seeds. Non-default options are passed directly to \code{\link{jags.model}}. In order to use the JAGS default initialization, set inits to NULL. 
#' See \code{\link{jags.model}} for more info.
#' @param mode String specifying what mode to use for the analysis. Choices are "quick", "report", or "paper".
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


nma.run.auto <- function(model,
                    monitor="DEFAULT",
                    DIC=TRUE,
                    inits = "DEFAULT",
                    mode = "report"){
  
  n.chains <- sma_set_mode(mode)$n.chains
  
  if(class(model) != "BUGSnetModel")
    stop("\'model\' must be a valid BUGSnetModel object created using the nma.model function.")
  
  if (!is.null(inits) && inits == "DEFAULT")
  {
    seeds <- sample(.Machine$integer.max, n.chains, replace = FALSE)
    inits <- list()
    for (i in 1:n.chains)
      inits[[i]] <- list(.RNG.seed = seeds[i], .RNG.name = "base::Mersenne-Twister")
  }
  
  # Set up monitor variables
  
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
  
  
  results <- sma_analyse_bayesian(sims = as_nlists(as_nlist(model$data)), # get data from BUGSnetModel object and put into compatible format
                                  code = model$bugs,
                                  monitor = new.monitor,
                                  inits=inits,
                                  mode=sma_set_mode(mode = mode),
                                  deviance = FALSE, # I think deviance is false bc the names of the things are different or something
                                  # path = ".",
                                 #  analysis = "analysis0000001",
                                  progress = FALSE,
                                  options = furrr::future_options(seed = TRUE))
  
  
  jagssamples <- as.mcmc.list(results$mcmcr1) # converting to mcmc.list for compatibility
  
  brun <- structure(list(samples=jagssamples,
                         model=model, 
                         scale=model$scale,
                         family=model$family,
                         link =model$link,
                         "trt.key"=as.character(t(model$trt.key[1]))),
                    class = "BUGSnetRun")
  return(brun)
  
}
