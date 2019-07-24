
#' Trace plots and convergence diagnostics for MCMC chains
#' @description Produces trace plots and Gelman-Rubin and Geweke convergence diagnostics for the MCMC chains obtained from 
#' \code{nma.run()}. The Gelman-Rubin and Geweke diagnostics are implemented using functions from the \code{coda} package. 
#' See \code{\link{gelman.diag}} and \code{\link{geweke.diag}} for more details.
#' @param nma An output produced by \code{nma.run()}
#' @param trace If TRUE, outputs trace plots. Default is TRUE.
#' @param gelman.rubin If TRUE, runs Gelman-Rubin diagnostic. Default is TRUE.
#' @param geweke If TRUE, runs Geweked diagnostic. Default is TRUE.
#' @param params Integer vector which specifies which parameters to produce trace plots for. Default is "all" which plots every parameter.
#' @param thin Thinning factor for the mcmc chains. Default is 1.
#' @param colours An optional vector of colors, one for each chain. 
#' @param plot_prompt If true, prompts the user to hit enter before plotting each additional batch of trace plots. Default is TRUE.
#' @param geweke_frac1 Fraction to use from beginning of chain. Default is 0.1.
#' @param geweke_frac2 Fraction to use from end of chain. Default is 0.5.
#' 
#' @return \code{gelman.rubin} An object of class \code{gelman.rubin.results} containing the Gelman-Rubin diagnostic results. 
#' A formatted table with custom PSRF threshold can be printed using \code{print(x, gelman.rubin.threshold = 1.2)}.
#' @return \code{geweke} An object of class \code{geweke.results} containing the Geweke diagnostic results. A formatted table 
#' with custom significance level can be printed using \code{print(x, alpha = 0.05)}.
#' 
#' @examples 
#' data(thrombolytic)
#' dich.slr <- data.prep(arm.data = cbind(thrombolytic, age), varname.t = "treatment", varname.s = "study")
#' random_effects_model <- nma.model(data=dich.slr, outcome="events", N="sampleSize", reference="SK",
#'                                   family="binomial", link="log", effects="random")
#' random_effects_results <- nma.run(random_effects_model, n.adapt=1000, n.burnin=1000, n.iter=10000)
#' nma.diag(random_effects_results)
#' @export
#' @seealso \code{\link{nma.run}}, \code{\link{gelman.diag}}, \code{\link{geweke.diag}}

nma.diag <- function(nma, 
                     trace = TRUE,
                     gelman.rubin = TRUE, 
                     geweke = TRUE,
                     params = "all", 
                     thin = 1, 
                     colours = "DEFAULT",
                     plot_prompt = TRUE,
                     geweke_frac1 = 0.1,
                     geweke_frac2 = 0.5)
{
  if (trace == FALSE && gelman.rubin == FALSE && geweke == FALSE)
    stop("At least one of the \'trace\', \'gelman_rubin\' or \'geweke\' parameters must be set to TRUE")
  
  #pull out column indices for the 'd' and 'sigma' parameters
  clist <- grep("(^d\\[[0-9]+\\])|(^sigma$)", colnames(nma$samples[[1]]))
  clist <- clist[-which(clist == grep("d\\[1\\]", colnames(nma$samples[[1]])))]
  mcmc.obj <- nma$samples[, clist, drop = FALSE]
  pnames <- colnames(mcmc.obj[[1]])
  
  #run gelman.diag from the coda package and produced formatted output
  if (gelman.rubin == TRUE)
  {
    gr <- gelman.diag(mcmc.obj)
    gr.obj <- structure(list(psrf = gr$psrf,
                             mpsrf = gr$mpsrf), class = "gelman.rubin.results")
  }

  #run geweke.diag from the coda package and produced formatted output
  if (geweke == TRUE)
  {
    gw <- geweke.diag(mcmc.obj, frac1 = geweke_frac1, frac2 = geweke_frac2)
    
    gwtbl <- list()
    for (i in 1:length(gw))
      gwtbl[[i]] <- gw[[i]]$z
    gwtbl <- do.call(cbind, gwtbl)
    colnames(gwtbl) <- paste0("Chain ", 1:length(gw))
    
    gw.obj <- structure(list(stats = gwtbl,
                             frac1 = geweke_frac1,
                             frac2 = geweke_frac2), class = "geweke.results")
  }
  
  #produce trace plots
  if (trace == TRUE)
  {
    samples <- do.call(rbind, mcmc.obj) %>% data.frame()
    n.iter <- nrow(mcmc.obj[[1]]) 
    n.chains <- nrow(samples)/n.iter
    
    if (params=="all"){
      params = 1:ncol(samples)
    } else if (!is.integer(params) || min(params) < 1 || max(params) > ncol(samples)) {
      stop(paste0("\'params\' argument must be \"all\" or an integer vector containing values from 1 to K ", 
                  "where K = [# treatment params] - 1 + [1 if random effects model]"))
    }
    
    if (colours=="DEFAULT"){
      colors=c("black", "red", "blue", "green", "purple", "yellow", "pink", "grey")#colors to be used in trace.plots
    } else {
      colors = colours
    }
    
    if(length(colors)<n.chains) stop("length(colours) must be no smaller than n.chains.")
    
    if (plot_prompt == TRUE)
      par(mfrow = c(min(4, length(params)), 2))
    else
      par(mfrow = c(length(params), 2))
    row_count <- 0
    for (i in params){
      if (plot_prompt == TRUE && row_count >= 4)
      {
        readline(prompt = "Press [ENTER] to continue plotting trace plots")
        dev.off()
        par(mfrow = c(min(4, length(params)), 2))
        row_count <- 0
      }
      
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
      row_count <- row_count + 1
    }
  }
  
  return(list(gelman.rubin = gr.obj, geweke = gw.obj))
}

print.gelman.rubin.results <- function(obj, gelman.rubin.threshold = 1.2)
{
  psrf <- obj$psrf
  psrf[,1] <- ifelse(psrf[,1] > gelman.rubin.threshold, sprintf("%.2f*", psrf[,1]), sprintf("%.2f", psrf[,1]))
  psrf[,2] <- sprintf("%.2f", as.numeric(psrf[,2]))
  
  mpsrf <- obj$mpsrf
  mpsrf <- ifelse(mpsrf > gelman.rubin.threshold, sprintf("%.2f*", mpsrf), sprintf("%.2f", mpsrf))
  
  cat("Gelman and Rubin's PSRF Convergence Diagnostic:\n--------------------------------\n")
  print.table(psrf)
  cat(paste0("--------------------------------\nMultivariate PSRF: ", mpsrf,
             "\n--------------------------------\n", paste0("* PSRF > ", gelman.rubin.threshold), "\n\n"))
}

print.geweke.results <- function(obj, alpha = 0.05)
{
  gwtbl <- list()
  for (i in 1:ncol(obj$stats))
    gwtbl[[i]] <- ifelse(abs(obj$stats[,i]) > qnorm(1 - alpha / 2, mean = 0, sd = 1), 
                       sprintf("%.2f*", obj$stats[,i]), sprintf("%.2f", obj$stats[,i]))
  gwtbl <- do.call(cbind, gwtbl)
  colnames(gwtbl) <- colnames(obj$stats)
  
  cat("Geweke\'s Convergence Diagnostic Results:\n--------------------------------\n")
  print.table(gwtbl)
  cat("--------------------------------",
      "\nFraction in 1st window = ", obj$frac1,
      "\nFraction in 2nd window = ", obj$frac2,
      paste0("\n* Statistically significant at the ", 100 * (1 - alpha), 
             "% confidence level"), sep = "")
}