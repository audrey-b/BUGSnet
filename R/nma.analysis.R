#' Run NMA model
#' @description Takes bugs code from an object produced by \code{nma.bugs} and runs model using \code{rjags}.
#' 
#' @param bugs Object produced by running \code{nma.bugs}.
#' @param monitor A list of all variables that you would like to monitor. The default is simply the treatment
#' effect samples ("d"). But you may want to monitor the deviance ("dev") as well.
#' 
#' @return \code{model} - A long character string containing the model that was run in \code{rjags}.
#' @return \code{bugsdata2}
#' @return \code{scale}
#' @return \code{trt.map.table} - Treatments mapped to integer numbers, used to run BUGS code.
#' @return \code{family} - Family that was used for the NMA model (e.g normal, binomial, poisson)
#' @return \code{link} - Link function that was used for the NMA model (e.g normal, binomial, poisson)
#' 

nma.analysis <- function(bugs,
                         monitor=c("d"),
                         n.adapt, 
                         n.burnin, 
                         n.iter, 
                         thin=1,
                         n.chains=3){
  
  
  jagsmodel <- jags.model(textConnection(bugs$model.str),
                          bugs$bugsdata2,
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
              trt.key=as.character(t(bugs$trt.map.table[1]))))
  
}