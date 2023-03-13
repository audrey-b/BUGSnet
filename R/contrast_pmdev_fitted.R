# Computing pmdev_fitted for contrast data
# y is the data.frame of y from MCMC sampling
# ytilde is a numeric vector of length sum(na)
# Omega data.frame of Omega from MCMC sampling
# nma BUGSnetRun object (has info on data)

#' @noRd
contrast_pmdev_fitted <- function(y, ytilde, Omega, nma) {
  
  ns <- nma$model$data$ns_c
  na <- nma$model$data$na_c
  pmdev_fitted <- numeric(ns)
  
  for (i in 1:ns) {
    
    # construct omega
    Om <- matrix(nrow = na[i]-1, ncol = na[i]-1)
    
    for(j in 1:(na[i]-1)) {
      for(k in 1:(na[i]-1)) {
        Om[j,k] <- Omega[1, names(Omega) == paste0("Omega.", i, ".", j, ".", k, ".")]
      }
    } 
    
    # construct y_i
    yi <- matrix(nrow = na[i]-1, ncol = 1)
    for(j in 2:na[i]) {
      yi[j-1] <- y[1, names(y) == paste0("y_c.", i, ".", j, ".")]
    }
    
    # construct ytilde_i
    ytildei <- yi
    for(j in 2:na[i]) {
      ytildei[j-1] <- ytilde[names(ytilde) == paste0("theta_c.", i, ".", j, ".")]
    }
    
    #calculate pmdev_fittedi
    pmdev_fitted[i] <- t(yi-ytildei) %*% Om %*% (yi - ytildei)
    
  }
  names(pmdev_fitted) <- paste0("pmdev_fitted.", seq(1, ns, by = 1))
  
  return(pmdev_fitted)
  
}