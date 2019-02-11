#' Traceplot of MCMC chains
#' @description Produces traceplots of the MCMC chains obtained from \code{nma.run()}
#' @param jagsoutput An output produced by \code{nma.run()}
#' @param n 

nma.trace <- function(jagsoutput, n="all"){
  samples <- do.call(rbind, jagsoutput$samples) %>% data.frame() %>% select(starts_with("d."))
  
  if (n=="all"){
    num_plots = ncol(samples)
  } else {
      num_plots = n
  }

  par(mfrow = c(num_plots, 2))
  for (i in 1:num_plots){
    plot(samples[,i], type="l", main = paste(names(samples[i])))
    plot(density(samples[,i]), main = paste(names(samples[i])))
  }
   
}