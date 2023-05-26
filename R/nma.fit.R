#' @title
#' Assess Model Fit
#' 
#' @description
#' Computes the Deviance Information Criteria (DIC) and produces a leverage plot (as in the NICE Technical Support Document 2)
#' for a given model. These can be used to assess and compare the fit of different models (i.e fixed vs random effects, consistency vs
#' inconsistency). \code{nma.fit} also produces a plot comparing the leverage of each data point against their contribution to
#' the total posterior deviance. Points lying outside the purple dotted line are generally identified as contributing to the model's poor fit.
#' Points with high leverage are influencial i.e. they have a stong influence on the estimates.
#' 
#' @param nma A \code{BUGSnetRun} object produced by running \code{nma.run()}. 
#' @param plot.pD Whether to include pD on the plot. Default is TRUE.
#' @param plot.DIC Whether to include DIC on the plot. Default is TRUE.
#' @param plot.Dres Whether to include Dres on the plot. Default is TRUE.
#' @param ... Graphical arguments such as main=, ylab=, and xlab= may be passed as in \code{plot()}. These arguments will only effect the
#' leverage plot.
#' 
#' @return \code{DIC} - A number indicating the Deviance Information Criteria. The DIC is calculated as the sum of \code{Dres} and \code{pD}. 
#' A larger DIC is indicative of a worse model fit.
#' @return \code{leverage} - A vector with one value per study arm indicating the leverage of each data point (study arm). Leverage is defined as \code{pmdev} minus the 
#' deviance at the posterior mean of the fitted values.
#' @return \code{w} - A vector with one value per study arm. The magnitude of \code{w} represents the data point's contribution to the posterior mean deviance of the 
#' model and is simply the square root of \code{pmdev}. The sign indicates whether the data point is being over (negative sign) or under (positive sign) estimated by 
#' the model and is calculated as the sign of the difference of the observed outcome minus the predicted outcome.
#' @return \code{pmdev} - A vector with one value per study arm representing the posterior mean residual deviance for each data point (study arm).
#' @return \code{Dres} - The posterior mean of the residual deviance.
#' @return \code{pD} - The effective number of parameters, calculated as the sum of the leverages.
#' 
#' @seealso
#' \code{\link{nma.run}}
#' 
#' @importFrom dplyr ends_with select slice starts_with
#' @importFrom graphics points text
#' @importFrom magrittr %>%
#' 
#' @examples
#' data(diabetes.sim)
#' 
#' diabetes.slr <- data.prep(
#'   arm.data = diabetes.sim, 
#'   varname.t = "Treatment", 
#'   varname.s = "Study"
#' )
#' 
#' #Random effects, consistency model.
#' #Binomial family, cloglog link. This implies that the scale will be the Hazard Ratio.
#' diabetes.re.c <- nma.model(
#'   data = diabetes.slr,
#'   outcome = "diabetes", 
#'   N = "n",
#'   reference = "Placebo",
#'   family = "binomial",
#'   link = "cloglog",
#'   effects = "random",
#'   type = "consistency",
#'   time = "followup"
#' )
#'  
#' diabetes.re.c.res <- nma.run(
#'   model = diabetes.re.c,
#'   n.adapt = 100,
#'   n.burnin = 0,
#'   n.iter = 100
#' )
#'
#' #Fixed effects, consistency model.
#' #Binomial family, cloglog link. This implies that the scale will be the Hazard Ratio.
#' diabetes.fe.c <- nma.model(
#'   data = diabetes.slr,
#'   outcome = "diabetes", 
#'   N = "n",
#'   reference = "Placebo",
#'   family = "binomial",
#'   link = "cloglog",
#'   effects = "fixed",
#'   type="consistency",
#'   time="followup"
#' )
#'  
#' diabetes.fe.c.res <- nma.run(
#'   model = diabetes.fe.c,
#'   n.adapt = 100,
#'   n.burnin = 0,
#'   n.iter = 100
#' )  
#' 
#' #Compare fixed vs random effects via leverage plots and DIC 
#' par(mfrow=c(1,2))
#' nma.fit(diabetes.fe.c.res, main = "Fixed Effects Model")
#' nma.fit(diabetes.re.c.res, main= "Random Effects Model")

#' @export
nma.fit  <- function(
  nma,
  plot.pD=TRUE,
  plot.DIC=TRUE,
  plot.Dres=TRUE,
  ...
){
  
  if (!inherits(nma, 'BUGSnetRun'))
    stop("\'nma\' must be a valid BUGSnetRun object created using the nma.run function.")
  
  jagssamples <- nma$samples
  
  # check what type of data included
  contrast <- ("y_c[1,1]" %in%colnames(jagssamples[[1]]))
  arm <- ("dev_a[1,1]" %in%colnames(jagssamples[[1]]))
  
  if (!inherits(jagssamples, 'mcmc.list'))
    stop('Object jagssamples must be of class mcmc.list')
  
  #stack all chains on top of each other
  samples = do.call(rbind, jagssamples) %>% data.frame()
  
  
  dev_a <- samples %>% select(., starts_with("dev_a"))
  dev_c <- samples %>% select(., starts_with("dev_c"))
  totresdev <- samples$totresdev %>% mean()
  pmdev_a <- colMeans(dev_a)
  pmdev_c <- colMeans(dev_c)
  
  if(arm) {
    
    if (nma$family == "binomial") {
      rhat <- samples %>% select(., starts_with("rhat"))
      r <- samples %>% select(., starts_with("r."))
      n <- samples %>% select(., starts_with("n"))
      rtilde <- rhat %>%
        colMeans() %>%
        matrix(nrow=1, ncol=ncol(rhat)) %>%
        data.frame() %>%
        slice(rep(1:n(), each = nrow(rhat)))
      
      pmdev_fitted <- 2*(r*log(r/rtilde)+(n-r)*log((n-r)/(n-rtilde)))[1,]
      
      if(TRUE %in% c(nma$model$data$r==0)) warning("Leverage cannot be calculated for zero cells.")
      
    } else if (nma$family == "poisson"){
      rhat <- samples %>% select(., starts_with("theta"))
      r <- samples %>% select(., starts_with("r.")) 
      n <- samples %>% select(., starts_with("n"))
      thetatilde <- rhat %>%
        colMeans() %>%
        matrix(nrow=1, ncol=ncol(rhat)) %>%
        data.frame() %>%
        slice(rep(1:n(), each = nrow(rhat)))
      pmdev_fitted <- 2*((thetatilde-r) + r*log(r/thetatilde))[1,]
      
      if(TRUE %in% c(nma$model$data$r==0)) warning("Leverage cannot be calculated for zero cells.")
      
    } else if (nma$family == "normal"){
      rhat <- samples %>% select(., starts_with("theta_a"))
      r <- samples %>% select(., starts_with("y."))
      prec <- samples %>% select(., starts_with("prec"))
      ytilde <- rhat %>%
        colMeans() %>%
        matrix(nrow=1, ncol=ncol(rhat)) %>%
        data.frame() %>%
        slice(rep(1:n(), each = nrow(rhat)))
      pmdev_fitted <- ((r-ytilde)*(r-ytilde)*prec)[1,]
    } 
    
    # calculate leverage for arms
    leverage_a = pmdev_a-pmdev_fitted
    colnames(leverage_a) <- colnames(pmdev_a)
    sign = sign(colMeans(r)-colMeans(rhat))
    w_a = sign*sqrt(as.numeric(pmdev_a))
  } else {
    leverage_a <- 0
    w_a <- NULL
  }

  if (contrast) {
    theta <- samples %>% select(., starts_with("theta_c"))
    y <- samples %>% select(., starts_with("y_c"))
    Omega <- samples %>% select(., starts_with("Omega"))
    ytilde <- colMeans(theta) # dim should be sum(na)
    pmdev_fitted_c <- contrast_pmdev_fitted(y, ytilde, Omega, nma)
    r <- y %>% select(., !ends_with(".1.")) # get rid of zeros
    rhat <- theta %>% select(., !ends_with(".1.")) # get rid of zeros
    
    # calculate leverage
    leverage_c <- pmdev_c-pmdev_fitted_c
    
    w_c <- sqrt(as.numeric(pmdev_c))
    
  } else {
    leverage_c <- 0
    w_c <- NULL
  }
  
  # Calculate effective parameters and DIC
  pD = sum(leverage_a, leverage_c)

  DIC = pD + totresdev
  
  
  
  eq <- function(x,c){c-x^2}
  x=seq(-3, 3, 0.001)
  c1=eq(x, c=rep(1, 6001))
  
  if(arm) {
    
    plot(w_a, leverage_a, xlab=expression('w'[ik]), ylab=expression('leverage'[ik]),
         ylim=c(0, max(c1+3, na.rm=TRUE)*1.15), xlim=c(-3,3), ...)
    
    if(contrast) {
      
      points(w_c, leverage_c) # add contrast points after plotting arm points
      # return both arm and contrast leverages, w, pmdev
      leverage <- c(as.numeric(leverage_a[1,]), leverage_c) 
      w <- c(w_a, w_c)
      pmdev <- c(pmdev_a, pmdev_c)
    } else {
      # only return arm leverage, w, pmdev
      leverage <- leverage_a[1,]
      w <- w_a
      pmdev <- pmdev_a
    }
    
  } else {
    
    plot(w_c, leverage_c, xlab=expression('w'[ik]), ylab=expression('leverage'[ik]),
         ylim=c(0, max(c1+3, na.rm=TRUE)*1.15), xlim=c(-3,3), ...)
    # only return contrast leverage, w, pmdev
    leverage <- leverage_c
    w <- w_c
    pmdev <- pmdev_c
  } 

  points(x, ifelse(c1 < 0, NA, c1),   lty=1, col="firebrick3",    type="l")
  points(x, ifelse(c1 < -1, NA, c1)+1, lty=2, col="chartreuse4",   type="l")
  points(x, ifelse(c1 < -2, NA, c1)+2, lty=3, col="mediumpurple3", type="l")
  points(x, ifelse(c1 < -3, NA, c1)+3, lty=4, col="deepskyblue3",  type="l")
  if (plot.pD ==TRUE){text(2, 4.3, paste("pD=", round(pD, 2)), cex = 0.8)}
  if (plot.Dres == TRUE){text(2, 3.9, paste("Dres=", round(totresdev,2)), cex = 0.8)}
  if (plot.DIC==TRUE){text(2, 3.5, paste("DIC=", round(DIC,2)), cex = 0.8)}
  
  # leverage <- c(as.numeric(leverage_a[1,]), leverage_c)
  # w <- c(w_a, w_c)
  # pmdev <- c(pmdev_a, pmdev_c)
  
  return(list(DIC=DIC,
              Dres = totresdev,
              pD = pD,
              leverage=leverage,
              w=w,
              pmdev=pmdev)
  )
}


#' @title
#' Consistency vs Inconsistency plot
#' 
#' @description
#' Plots the posterior mean deviance of a consistency model vs an inconsistency model. This plot can help identify loops
#' where inconsistency is present. Ideally, both models will contribute approximately 1 to the posterior mean deviance.
#' 
#' @param consistency.model.fit Results of \code{nma.fit()} of an consistency model.
#' @param inconsistency.model.fit Results of \code{nma.fit()} of an inconsistency model.
#' @param ... Graphical arguments such as main=, ylab=, and xlab= may be passed in \code{plot()}.
#' 
#' @seealso
#' \code{\link{nma.run}}
#' 
#' @importFrom graphics points
#' 
#' @examples
#' data(diabetes.sim)
#' diabetes.slr <- data.prep(
#'   arm.data = diabetes.sim, 
#'   varname.t = "Treatment", 
#'   varname.s = "Study"
#' )
#' 
#' #Random effects, consistency model.
#' #Binomial family, cloglog link. This implies that the scale will be the Hazard Ratio.
#' diabetes.re.c <- nma.model(
#'   data = diabetes.slr,
#'   outcome = "diabetes",
#'   N = "n",
#'   reference = "Placebo",
#'   family = "binomial",
#'   link = "cloglog",
#'   effects = "random",
#'   type = "consistency",
#'   time = "followup"
#' )
#'  
#' diabetes.re.c.res <- nma.run(
#'   model = diabetes.re.c,
#'   n.adapt = 100,
#'   n.burnin = 0,
#'   n.iter = 100
#' )
#'
#' #Random effects, inconsistency model.
#' #Binomial family, cloglog link. This implies that the scale will be the Hazard Ratio.
#' diabetes.re.i <- nma.model(
#'   data = diabetes.slr,
#'   outcome = "diabetes", 
#'   N = "n",
#'   reference = "Placebo",
#'   family = "binomial",
#'   link = "cloglog",
#'   effects = "random",
#'   type = "inconsistency",
#'   time = "followup"
#' )
#'  
#' diabetes.re.i.res <- nma.run(
#'   model = diabetes.re.i,
#'   n.adapt = 100,
#'   n.burnin = 0,
#'   n.iter = 100
#' )  
#' 
#' # Assess model fit for a both an inconsistency model and consistency model using nma.fit()
#' assess.consistency <- nma.fit(diabetes.re.c.res)
#' assess.inconsistency <- nma.fit(diabetes.re.i.res)
#' 
#' #Plot the results against each other to assess inconsistency
#' nma.compare(assess.consistency, assess.inconsistency)

#' @export
nma.compare <- function(
  consistency.model.fit,
  inconsistency.model.fit,
  ...
){
  x=as.numeric(consistency.model.fit$pmdev)
  y=as.numeric(inconsistency.model.fit$pmdev)
  upp <- max(0, max(x,y)*1.04)
  ptsline <- seq(0,upp, length.out=2000)
  plot(x, y,
       xlim=c(0, upp),
       ylim=c(0, upp),
       ylab="Inconsistency model", xlab="Consistency model", pch=16, ...)
  points(ptsline, ptsline, type = "l", lty=2)
}
