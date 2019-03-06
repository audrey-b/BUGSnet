#' Traceplot of MCMC chains
#' @description Produces traceplots of the MCMC chains obtained from \code{nma.run()}
#' @param nma An output produced by \code{nma.run()}
#' @param n Limits the number of printed variables to the first \code{n}
#' @export

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