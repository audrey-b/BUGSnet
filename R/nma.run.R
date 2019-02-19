#' Run NMA model
#' @description Takes bugs code from an object produced by \code{nma.model} and runs model using \code{rjags}.
#' 
#' @param bugs Object produced by running \code{nma.model}.
#' @param monitor A list of all variables that you would like to monitor. The default is simply the treatment
#' effect samples ("d"). But you may want to monitor the deviance ("dev") as well.
#' 
#' @return \code{model} - A long character string containing the model that was run in \code{rjags}.
#' @return \code{data} - The data used in the BUGS code.
#' @return \code{scale} - The scale of the outcome, based on the chosen family and link function.
#' @return \code{trt.key} - Treatments mapped to numbers, used to run BUGS code.
#' @return \code{family} - Family that was used for the NMA model (e.g normal, binomial, poisson)
#' @return \code{link} - Link function that was used for the NMA model (e.g normal, binomial, poisson)
#' 

nma.run <- function(bugs,
                         monitor=c("d"),
                         n.adapt, 
                         n.burnin, 
                         n.iter, 
                         thin=1,
                         n.chains=3){
  
  
  jagsmodel <- jags.model(textConnection(bugs$model),
                          bugs$data,
                          n.chains=n.chains,
                          n.adapt=n.adapt)
  
  jagsburnin <- update(jagsmodel,
                       n.iter=n.burnin)
  
  jagssamples <- coda.samples(jagsmodel, 
                              variable.names=monitor,
                              n.iter=n.iter,
                              thin=thin)
  
  
  # print("The baseline treatment was ...")
  
  return(list(model=jagsmodel, 
              samples=jagssamples, 
              scale=bugs$scale,
              family=bugs$family,
              link =bugs$link,
              "trt.key"=as.character(t(bugs$trt.key[1]))))
  
}