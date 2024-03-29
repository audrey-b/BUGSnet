#' @title
#' Deprecated. Please use \code{nma.diag} instead.
#' 
#' @description
#' Deprecated. Produces traceplots of the MCMC chains obtained from \code{nma.run()}
#' 
#' @param nma A \code{BUGSnetRun} object produced by running \code{nma.run()}.
#' @param n Integer which limits the number of printed variables to the first \code{n}. Default is "all" which plots every variable.
#' @param thin Thinning factor for the mcmc chains. Default is 1.
#' @param colours An optional vector of colors, one for each chain.
#' 
#' @seealso
#' \code{\link{nma.diag}}
#' 
#' @importFrom dplyr select starts_with
#' @importFrom graphics lines par
#' @importFrom magrittr %>% %<>%
#' @importFrom stats density

#' @export
nma.trace <- function(
  nma,
  n = "all",
  thin = 1,
  colours = "DEFAULT"
){
  .Deprecated("nma.diag")
  
  if (!inherits(nma, 'BUGSnetRun'))
    stop("\'nma\' must be a valid BUGSnetRun object created using the nma.run function.")
  
  samples <- do.call(rbind, nma$samples) %>% data.frame()
  
  n.iter <- nrow(nma$samples[[1]]) 
  n.chains <- nrow(samples)/n.iter
  
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

  if (colours=="DEFAULT"){
    colors=c("black", "red", "blue", "green", "purple", "yellow", "pink", "grey")#colors to be used in trace.plots
  } else {
    colors = colours
  }
  
  if(length(colors)<n.chains) stop("length(colours) must be no smaller than n.chains.")
  
  par(mfrow = c(num_plots, 2), mar=c(1,1,1,1))
  for (i in 1:num_plots){
    plot(seq(thin, n.iter, thin), #ensures x-axis labels correspond to correct iteration
         samples[seq(thin, n.iter, thin),i], 
         type="l", 
         main = paste(names(samples[i])), 
         ylab="",
         xlab="iteration",
         col = colors[1],
         ylim=c(min(samples[seq(thin, nrow(samples), thin),i]), 
                max(samples[seq(thin, nrow(samples), thin),i]))) #ensure that each chain fits in plot
    
    if (n.chains>1){ 
      for (j in 2:n.chains){
          lines(seq(thin, n.iter, thin), samples[seq((j-1)*n.iter+thin, j*n.iter, thin), i ], col=colors[j])
      }
    }  
    plot(density(samples[,i]), main = paste(names(samples[i])))
  }
}