#' Run NMA model
#' @description Takes bugs code from an object produced by \code{nma.model} and runs model using \code{jags}.
#' 
#' @param model Object produced by running \code{nma.model}.
#' @param monitor A list of all variables that you would like to monitor. The default is simply the treatment
#' effect samples ("d"). But you may want to monitor the deviance ("dev") as well.
#' @param n.adapt Number of adaptations for the mcmc chains.
#' @param n.burnin Number of burnin iterations for the mcmc chains.
#' @param n.iter Number of iterations for the mcmc chains.
#' @param thin Thinning factor for the mcmc chains. Default is 1.
#' @param n.chains Number of mcmc chains. Default is 3.
#' 
#' @return \code{model} - The object obtained from \code{nma.model} that was used to run \code{jags}.
#' @return \code{data} - The data used with the BUGS code.
#' @return \code{scale} - The scale of the outcome, based on the chosen family and link function.
#' @return \code{trt.key} - Treatments mapped to numbers, used to run BUGS code.
#' @return \code{family} - Family that was used for the NMA model (e.g normal, binomial, poisson)
#' @return \code{link} - Link function that was used for the NMA model (e.g normal, binomial, poisson)
#' 

nma.run <- function(model,
                         monitor=c("d"),
                         n.adapt, 
                         n.burnin=0, 
                         n.iter, 
                         thin=1,
                         n.chains=3){
  
  
  jagsmodel <- jags.model(textConnection(model$bugs),
                          model$data,
                          n.chains=n.chains,
                          n.adapt=n.adapt)
  
  if(n.burnin!=0) jagsburnin <- update(jagsmodel, n.iter=n.burnin)
  
  jagssamples <- coda.samples(jagsmodel, 
                              variable.names=monitor,
                              n.iter=n.iter,
                              thin=thin)
  
  
  # print("The baseline treatment was ...")
  
  return(list(samples=jagssamples,
              model=model, 
              scale=model$scale,
              family=model$family,
              link =model$link,
              "trt.key"=as.character(t(model$trt.key[1]))))
  
}