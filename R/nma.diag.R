
#' Trace plots and convergence diagnostics for MCMC chains
#' @description Produces trace plots and Gelman-Rubin and Geweke convergence diagnostics for the MCMC chains obtained from 
#' \code{nma.run()}. The Gelman-Rubin and Geweke diagnostics are implemented using functions from the \code{coda} package. 
#' @param nma A \code{BUGSnetRun} object produced by \code{nma.run()}
#' @param trace If TRUE, outputs trace plots. Default is TRUE.
#' @param gelman.rubin If TRUE, runs Gelman-Rubin diagnostic. Default is TRUE.
#' @param geweke If TRUE, runs Geweke diagnostic. Default is TRUE.
#' @param params Integer or character vector which specifies which parameters to produce trace plots for when trace is set to TRUE. 
#' Default is "all" which plots every monitored parameter.
#' @param thin Thinning factor for the mcmc chains when producing trace plots. Default is 1.
#' @param nrow Number rows in each batch of trace plots
#' @param ncol Number of columns in each batch of trace plots
#' @param plot_prompt If TRUE, prompts the user to hit enter before plotting each additional batch of trace plots. Default is TRUE.
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
#' dich.slr <- data.prep(arm.data = thrombolytic, varname.t = "treatment", 
#'                       varname.s = "study")
#' random_effects_model <- nma.model(data=dich.slr, outcome="events", 
#'                                   N="sampleSize", reference="SK",
#'                                   family="binomial", link="log", 
#'                                   effects="random")
#' random_effects_results <- nma.run(random_effects_model, n.adapt=100, 
#'                                   n.burnin=0, n.iter=100)
#' nma.diag(random_effects_results)
#' @export
#' @seealso \code{\link{nma.run}}

nma.diag <- function(nma, 
                     trace = TRUE,
                     gelman.rubin = TRUE, 
                     geweke = TRUE,
                     params = "all",
                     thin = 1, 
                     ncol = 1,
                     nrow = 3,
                     plot_prompt = TRUE,
                     geweke_frac1 = 0.1,
                     geweke_frac2 = 0.5)
{
  
  # Bind variables to function
  str_sub <- NULL
  str_extract <- NULL
  gelman.diag <- NULL
  geweke.diag <- NULL
  
  if (class(nma) != "BUGSnetRun")
    stop("\'nma\' must be a valid BUGSnetRun object created using the nma.run function.")
  
  if (trace == FALSE && gelman.rubin == FALSE && geweke == FALSE)
    stop("At least one of the \'trace\', \'gelman_rubin\' or \'geweke\' parameters must be set to TRUE")
  
  #pull out column indices for the 'd' and 'sigma' parameters
  if ((length(params) == 1 && params == "all") || is.numeric(params) == TRUE) {
    clist <- grep("(^d\\[[0-9,]+\\])|(^sigma$)", colnames(nma$samples[[1]]))
    clist <- clist[colnames((nma$samples[[1]])[,clist]) != "d[1]"]
    #quick fix for inconsistency model
    if (nma$model$type == "inconsistency")
    {
      i1 <- str_sub(str_extract(colnames((nma$samples[[1]])[,clist]), "d\\[[0-9]+"), 3, -1)
      i2 <- str_sub(str_extract(colnames((nma$samples[[1]])[,clist]), ",[0-9]+\\]"), 2, -2)
      clist <- clist[is.na(i1) | (i1 != i2)]
    }
    if (is.numeric(params) == TRUE) 
      clist <- clist[params]
    mcmc.obj <- nma$samples[, clist, drop = FALSE]
  } else {
    mcmc.obj <- nma$samples[, params, drop = FALSE]
  }
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
    samples$chain <- as.factor(rep(1:n.chains, rep(n.iter, n.chains)))
    samples$iteration <- rep(1:n.iter, n.chains)
    
    thinned_index <- NULL
    for (i in 1:n.chains)
      thinned_index <- c(thinned_index, (i-1) * n.iter + seq(1, n.iter, thin))
    samples <- samples[thinned_index,]
    
    traceplots <- list()
    for (i in 1:length(pnames))
    {
      traceplots[[2*i-1]] <- ggplot(samples, aes_string(x = "iteration", y = make.names(pnames)[i], col = "chain")) +
        geom_line() + xlab("") + ylab("") + ggtitle(pnames[i]) + 
        theme(plot.title = element_text(hjust = 0.5), legend.position = "none")
      traceplots[[2*i]] <- ggplot(samples, aes_string(x = make.names(pnames)[i], col = "chain")) +
        geom_density() + xlab("") + ggtitle(pnames[i]) + 
        theme(plot.title = element_text(hjust = 0.5), legend.position = "none")
    }
    
    for (i in 1:(ceiling(length(pnames) / (nrow * ncol))))
    {
      if (plot_prompt == TRUE && i > 1)
      {
        pstr <- readline(prompt = "Press [ENTER] to continue plotting trace plots (or type \'stop\' to end plotting)> ")
        if (trimws(tolower(pstr), which = "both") == "stop")
          break
      }
      gridExtra::grid.arrange(grobs = traceplots[((i-1)*2*nrow*ncol+1):min((i*2*nrow*ncol),length(traceplots))], ncol = 2 * ncol)
    }
  }
  
  if(gelman.rubin == TRUE && geweke == TRUE){
    return(list(gelman.rubin = gr.obj, geweke = gw.obj))
  }else if(gelman.rubin == TRUE){
    return(list(gelman.rubin = gr.obj))
  }else if(geweke == TRUE){
    return(list(geweke = gw.obj))
  }

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