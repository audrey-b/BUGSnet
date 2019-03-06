#' Traceplot of MCMC chains
#' @description Produces traceplots of the MCMC chains obtained from \code{nma.run()}
#' @param nma An output produced by \code{nma.run()}
#' @param n Integer which limits the number of printed variables to the first \code{n}. Default is "all" which plots every variable.
#' @export
#' @seealso \code{\link{nma.run}}
#' 
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
#'n.adapt=1000,
#'n.burnin=1000,
#'n.iter=10000)
#'
#'nma.trace(diabetes.re.c.res, n=4)
#'  

nma.trace <- function(nma, n="all"){
  samples <- do.call(rbind, nma$samples) %>% data.frame()
  
  # samples.vars <- colnames(samples)
  # scalars <- vars[which(vars %in% samples.vars)]
  # samples.matr <- setdiff(samples.vars, scalars)
  # 
  # samples %<>% select(starts_with("d."), scalars)
  
  
  if("sigma" %in% colnames(samples)) {samples %<>% select(starts_with("d."), "sigma")
    } else {samples %<>% select(starts_with("d."))}
  
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