#' Forest plot
#' @description Produces a forest plot of point estimates and 95% credible intervals obtained with the quantile method.
#' @param nma An output produced by \code{nma.run()}
#' @param comparator The treatment to use as a comparator
#' @param central.tdcy The posterior statistic used in order to measure relative effectiveness. The options are "mean" and "median". Default is median.
#' @param lwd Line width relative to the default (default=1).
#' @param x.trans Optional. A string indicating a transformation to apply to the x-axis. Setting this parameter to "log" is useful when there are extreme values or to allow an easier interpretation of odds ratios and relative ratios (if e.g. treatment B is twice as far from the line y=1 then treatment A then it's OR/RR is twice that of treatment A.) 
#' @param log.scale If TRUE, odds ratios, relative risk or hazard ratios are reported on the log scale. Default is FALSE.
#'
#' @return \code{forestplot} - A forest plot.
#'
#' @examples
#' 
#' #make forest plot
#' nma.forest(nma = nma_results, comparator="Placebo")
#' 


nma.forest <- function(nma, 
                           comparator, 
                           central.tdcy = "median", 
                       log.scale=FALSE,
                           line.size=1,
                           x.trans=NULL) {
  
  x2 <- do.call(rbind, nma$samples) %>% data.frame() %>% select(starts_with("d."))
  trt.names <- nma$trt.key
  colnames(x2) <- trt.names
  
  x3 <- x2
  new.vars <- paste0(colnames(x2), "-", comparator)
  for(i in 1:ncol(x2)) {
    x3[[new.vars[i]]] <- x2[, i] - x2[,comparator]
  }
  
  x3 %<>% select(new.vars)
  colnames(x3) <- trt.names
  x3 %<>% select(-comparator)
  
  if(log.scale==FALSE & nma$link!="identity"){

  tmp.mean <- x3 %>%  
    summarise_all(funs(mean = exp.mean)) %>% gather() %>%
    rename(trt = key, mean = value) %>%
    mutate(trt = sub("_mean", "", trt))
  
  tmp.lci <- x3 %>%  
    summarise_all(funs(lci = exp.lci)) %>% gather() %>%
    rename(trt = key, lci = value) %>%
    mutate(trt = sub("_lci", "", trt))
  
  tmp.uci <- x3 %>%  
    summarise_all(funs(uci = exp.uci)) %>% gather() %>%
    rename(trt = key, uci = value) %>%
    mutate(trt = sub("_uci", "", trt))
  
  null.value <- 1
  log.str<-""
  } else{
    
    tmp.mean <- x3 %>%  
      summarise_all(funs(mean = id.mean)) %>% gather() %>%
      rename(trt = key, mean = value) %>%
      mutate(trt = sub("_mean", "", trt))
    
    tmp.lci <- x3 %>%  
      summarise_all(funs(lci = id.lci)) %>% gather() %>%
      rename(trt = key, lci = value) %>%
      mutate(trt = sub("_lci", "", trt))
    
    tmp.uci <- x3 %>%  
      summarise_all(funs(uci = id.uci)) %>% gather() %>%
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
  


f.plot <- ggplot(tmp1, aes(x=trt, y=mean, ymin=lci, ymax=uci)) +
     geom_pointrange(size=line.size) +
     geom_hline(yintercept=null.value,lty=2) +
     scale_x_discrete(limits = sort(tmp1$trt, decreasing=TRUE)) +
     xlab("Treatment") +
     ylab(paste0(log.str, nma$scale, " relative to ",comparator,
                 "\n(showing posterior ", central.tdcy," with 95% CrI)")) +
     coord_flip() +
     theme_classic() #+
     #labs(caption = paste("note: each treatment compared to", comparator))

if(is.null(x.trans)){
f.plot <- f.plot +
  scale_y_continuous(breaks = scales::pretty_breaks(c(tmp1$lci,tmp1$uci), n = 10))
}else{
  f.plot <- f.plot +
    scale_y_continuous(trans=x.trans,
                       breaks = scales::pretty_breaks(c(min(tmp1$lci),max(tmp1$uci)), n = 10))
}

return("forestplot"=f.plot)

}

