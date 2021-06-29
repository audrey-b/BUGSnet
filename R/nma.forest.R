#' Forest plot
#' @description Produces a forest plot of point estimates and 95% credible intervals obtained with the quantile method.
#' @param nma A \code{BUGSnetRun} object produced by running \code{nma.run()}.
#' @param comparator The treatment to use as a comparator
#' @param central.tdcy The posterior statistic used in order to measure relative effectiveness. The options are "mean" and "median". Default is median.
#' @param order Optional. A vector of strings representing the order in which to display the treatments.
#' @param lwd Line width relative to the default (default=1).
#' @param x.trans Optional. A string indicating a transformation to apply to the x-axis. Setting this parameter to "log" is useful when there are extreme values or to allow an easier interpretation of odds ratios and relative ratios (if e.g. treatment B is twice as far from the line y=1 then treatment A then it's OR/RR is twice that of treatment A.) 
#' @param log.scale If TRUE, odds ratios, relative risk or hazard ratios are reported on the log scale. Default is FALSE.
#' @param cov.value  Must be specified for meta-regression. This is the value of the covariate for which to report the results.

#' @return \code{forestplot} - A forest plot.
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
#'n.adapt=100,
#'n.burnin=0,
#'n.iter=100)
#'      
#'        
#' #make forest plot
#' nma.forest(nma = diabetes.re.c.res, comparator="Placebo")
#' @export
#' @seealso \code{\link{nma.run}}, \code{\link{nma.league}}, \code{\link{nma.rank}} 


nma.forest <- function(nma, 
                       comparator, 
                       central.tdcy = "median", 
                       order = NULL,
                       log.scale=FALSE,
                       lwd=1,
                       x.trans=NULL,
                       cov.value=NULL) {
  
  # Bind variables to function
  key <- NULL
  value <- NULL
  trt <- NULL
  lci <- NULL
  uci <- NULL
  
  if (class(nma) != "BUGSnetRun")
    stop("\'nma\' must be a valid BUGSnetRun object created using the nma.run function.")
  
  if(!is.null(nma$model$covariate) & is.null(cov.value)) stop("cov.value must be specified for meta-regression")
  
  if(nma$model$type=="inconsistency") {
    stop('This function has not been implemented yet for the inconsistency model.')
  } 
  
  x2 <- do.call(rbind, nma$samples) %>% data.frame() %>% select(starts_with("d."))
  trt.names <- nma$trt.key
  colnames(x2) <- trt.names
  
  if(!is.null(cov.value)){#meta-regression
    
    betamat <- do.call(rbind, nma$samples) %>% data.frame() %>% select(starts_with("beta."))
    trt.names <- nma$trt.key
    colnames(betamat) <- trt.names
    
    x2 <- x2+betamat*(cov.value-nma$model$mean.cov)
  }
  
  x3 <- x2
  new.vars <- paste0(colnames(x2), "-", comparator)
  for(i in 1:ncol(x2)) {
    x3[[new.vars[i]]] <- x2[, i] - x2[,comparator]
  }
  
  x3 %<>% select(new.vars)
  colnames(x3) <- trt.names
  x3 %<>% select(-comparator)
  
  if(log.scale==FALSE & nma$link!="identity"){
    
    if(central.tdcy=="mean"){
      tmp.mean <- x3 %>%  
        summarise_all(list(mean = exp.mean)) %>% gather() %>%
        rename(trt = key, mean = value) %>%
        mutate(trt = sub("_mean", "", trt))
    }else if(central.tdcy=="median"){
      tmp.mean <- x3 %>%  
        summarise_all(list(mean = exp.median)) %>% gather() %>%
        rename(trt = key, mean = value) %>%
        mutate(trt = sub("_mean", "", trt))
    }else stop("central.tdcy must be mean or median")
    
    tmp.lci <- x3 %>%  
      summarise_all(list(lci = exp.lci)) %>% gather() %>%
      rename(trt = key, lci = value) %>%
      mutate(trt = sub("_lci", "", trt))
    
    tmp.uci <- x3 %>%  
      summarise_all(list(uci = exp.uci)) %>% gather() %>%
      rename(trt = key, uci = value) %>%
      mutate(trt = sub("_uci", "", trt))
    
    
    null.value <- 1
    log.str<-""
    
  } else{
    
    if(central.tdcy=="mean"){
      tmp.mean <- x3 %>%  
        summarise_all(list(mean = id.mean)) %>% gather() %>%
        rename(trt = key, mean = value) %>%
        mutate(trt = sub("_mean", "", trt))
    }else if(central.tdcy=="median"){
      tmp.mean <- x3 %>%  
        summarise_all(list(mean = id.median)) %>% gather() %>%
        rename(trt = key, mean = value) %>%
        mutate(trt = sub("_mean", "", trt))
    }else stop("central.tdcy must be mean or median")
    
    
    
    tmp.lci <- x3 %>%  
      summarise_all(list(lci = id.lci)) %>% gather() %>%
      rename(trt = key, lci = value) %>%
      mutate(trt = sub("_lci", "", trt))
    
    tmp.uci <- x3 %>%  
      summarise_all(list(uci = id.uci)) %>% gather() %>%
      rename(trt = key, uci = value) %>%
      mutate(trt = sub("_uci", "", trt))
    
    null.value <- 0
    
    if(nma$link=="identity"){
      log.str <- ""
    } else{
      log.str <- "Log "  
    }
    
  }
  
  tmp1 <- left_join(tmp.mean, tmp.lci, by = "trt") %>%
    left_join(., tmp.uci, by = "trt") %>% data.frame()
  
  if(!is.null(cov.value)){#meta-regression
    cov.str <- paste0(" when ",nma$model$covariate,"=",cov.value)
  }else cov.str <- ""
  
  #Sys.setlocale("LC_COLLATE","C") #for sorting mix of capital and small cap
  if(is.null(order)){order <- sort(tmp1$trt, decreasing=TRUE)
  }else if(comparator %in% order){order <- order[-which(order==comparator)]
  }
  
  f.plot <- ggplot(tmp1, aes(x=trt, y=mean, ymin=lci, ymax=uci)) +
    geom_pointrange(size=lwd) +
    geom_hline(yintercept=null.value,lty=2) +
    scale_x_discrete(limits = order) +
    xlab("Treatment") +
    ylab(paste0(log.str, nma$scale, " relative to ",comparator, cov.str,
                "\n(showing posterior ", central.tdcy," with 95% CrI)")) +
    coord_flip() +
    theme_classic() #+
  #labs(caption = paste("note: each treatment compared to", comparator))
  
  if(is.null(x.trans)){
    f.plot <- f.plot +
      scale_y_continuous(breaks = scales::pretty_breaks(n=10 ))
  }
  else{
    
    f.plot <- f.plot +
      scale_y_continuous(trans=x.trans,
                         breaks = scales::pretty_breaks(n=10))
  }
  
  return("forestplot"=f.plot)
  
}

