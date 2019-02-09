#' Assess Model Fit
#' @description Computes the Deviance Information Criteria (DIC) and produces a leverage plot (as in the NICE Technical Support Document 2)
#' for a given model. These can be used to assess and compare the fit of different models (i.e fixed vs random effects, consistency vs
#' inconsistency). \code{nma.fit} also produces a plot comparing the leverage of each data point against their contribution to
#' the total posterior deviance. Any point lying outside the purple dotted line is considered to be an "outlier";
#' contributing to the model's poor fit.
#' 
#' @param jagsoutput Resulting dataset from running \code{nma.analysis()}. 
#' @param ... Graphical arguments such as main=, ylab=, and xlab= may be passed as in \code{plot()}. These arguments will only effect the
#' leverage plot.
#' 
#' @return \code{DIC} - A number indicating the Deviance Information Criteria. A larger DIC is indicative of a worse model fit.
#' @return \code{leverage} - A vector indicating the leverage of each data point (study arm). Leverage is defined as the posterior mean of the
#' residual deviance minus the deviance at the posterior mean of the fitted values.
#' @return \code{w} - A vector with 1 value per study arm. The magnitude of \code{w} represents the data point's contribution to the
#' posterior mean deviance of the model. The sign indicates whether the data is being over or under estimated by the model.
#' @return \code{pmdev} - 
#' 
#' @examples
#' #Compare fixed vs random effects via leverage plots and DIC
#' #fixed_effects_results are the outputs of nma.analysis() where the parameter effects="fixed"
#' #random_effects_results are the outputs of nma.analysis() where the parameter effects="random"
#' 
#' par(mfrow=c(1,2))
#' nma.fit(fixed_effects_results, main = "Fixed Effects Model" )
#' nma.fit(random_effects_results, main= "Random Effects Model")

nma.fit  <- function(jagsoutput, plot.pD=TRUE, plot.DIC=TRUE, plot.Dres=TRUE, ...){
  jagssamples <- jagsoutput$samples
  
  if (class(jagssamples) != "mcmc.list"){stop('Object jagssamples must be of class mcmc.list')}
  
  #stack all chains on top of each other
  samples = do.call(rbind, jagssamples) %>% data.frame()
  
  r <- samples %>% select(., starts_with("r.")) 
  n <- samples %>% select(., starts_with("n"))
  dev <- samples %>% select(., starts_with("dev"))
  totresdev <- samples$totresdev %>% mean()
  pmdev <- colMeans(dev)

  if (jagsoutput$family == "binomial") {
    rhat <- samples %>% select(., starts_with("rhat"))
    rtilde <- rhat %>%
      colMeans() %>%
      matrix(nrow=1, ncol=ncol(rhat)) %>%
      data.frame() %>%
      slice(rep(1:n(), each = nrow(rhat)))
    
    pmdev_fitted <- 2*(r*log(r/rtilde)+(n-r)*log((n-r)/(n-rtilde)))[1,]
    
  } else if (jagsoutput$family == "poisson"){
    rhat <- samples %>% select(., starts_with("theta"))
    thetatilde <- rhat %>%
      colMeans() %>%
      matrix(nrow=1, ncol=ncol(rhat)) %>%
      data.frame() %>%
      slice(rep(1:n(), each = nrow(rhat)))
    pmdev_fitted <- 2*((thetatilde-r) + r*log(r/thetatilde))[1,]
  }
  
  
  leverage = pmdev-pmdev_fitted
  
  DIC = sum(leverage) + totresdev
  
  sign = sign(colMeans(r)-colMeans(rhat))
  
  w = sign*sqrt(as.numeric(pmdev))
  
  pD = sum(leverage)
  eq = function(x,c){c-x^2}
  x=seq(-3, 3, 0.001)
  c1=eq(x, c=rep(1, 6001))
  
  plot(w, leverage, xlab=expression('w'[ik]), ylab=expression('leverage'[ik]),
       ylim=c(0, 4.5), xlim=c(-3,3), ...)
  lines(x, c1,   lty=1, col="firebrick3",    lwd =2)
  lines(x, c1+1, lty=2, col="chartreuse4",   lwd =2)
  lines(x, c1+2, lty=3, col="mediumpurple3", lwd =2)
  lines(x, c1+3, lty=4, col="deepskyblue3",  lwd =2)
  if (plot.pD ==TRUE){text(2, 4.3, paste("pD=", round(pD, 2)), cex = 0.8)}
  if (plot.Dres == TRUE){text(2, 3.9, paste("Dres=", round(totresdev,2)), cex = 0.8)}
  if (plot.DIC==TRUE){text(2, 3.5, paste("DIC=", round(DIC,2)), cex = 0.8)}
  
  return(list(DIC=DIC,
              leverage=leverage,
              w=w,
              pmdev=pmdev,
              Dres = totresdev,
              pD = pD)
  )
}

#' Consistency vs Inconsistency plot
#' @description Plots the posterior mean deviance of a consistency model vs an inconsistency model. This plot can help identify loops
#' where inconsistency is present. Ideally, both models will contribute approximately 1 to the posterior mean deviance.
#' 
#' @param consistency.model.fit Results of \code{nma.fit()} of an consistency model.
#' @param inconsistency.model.fit Results of \code{nma.fit()} of an inconsistency model.
#' @param ... Graphical arguments such as main=, ylab=, and xlab= may be passed as in \code{plot()}.
#' 
#' @examples
#' # Assess model fit for a both an inconsistency model and consistency model using nma.fit()
#' assess.consistency <- nma.fit(consistency_results)
#' assess.inconsistency <- nma.fit(inconsistency_results)
#' 
#' #Plot the results against each other to assess inconsistency
#' inconsistency.plot(consistency_results, inconsistency_results)


nma.compare <- function(consistency.model.fit, inconsistency.model.fit, ...){
  x=as.numeric(consistency.model.fit$pmdev)
  y=as.numeric(inconsistency.model.fit$pmdev)
  plot(x, y,
       xlim=c(0, max(0, max(x,y) + 0.5)),
       ylim=c(0, max(0, max(x,y) + 0.5)),
       ylab="Inconsistency model", xlab="Consistency model", pch=16, ...)
  lines(-1:100,-1:100, lty = 2 )
}