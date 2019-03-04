#' Plot of relative treatment effects vs covariate values for meta-regression.
#' @description Produces a plot of relative treatment effects on the linear scale vs covariate values for meta-regression.
#' 
#' @param nma Resulting output from running \code{nma.run()}.
#' @param lwd Line width relative to the default (default=1).
#' @param palette A string indicating the colour set from RcolorBrewer for the plot. "set1" is great, but you may need a different colour set if 
#' there are many treatments in your network.
#' 
#' @return \code{regplot} - A plot of the relative treatment effects vs covariate values for meta-regression.


nma.regplot <- function(nma, x.range=NULL, lwd=1, palette="Set1"){
  
  if(is.null(x.range)){
    x.min=min(nma$model$data$x, na.rm=TRUE) + nma$model$mean.cov
    x.max=max(nma$model$data$x, na.rm=TRUE) + nma$model$mean.cov
    x.range=c(x.min,x.max)
  }
  
  
  dmat <- do.call(rbind, nma$samples) %>% data.frame() %>% select(starts_with("d."))
  trt.names <- nma$trt.key
  colnames(dmat) <- trt.names
  
  dmat %<>% select(-matches(nma$model$reference))
  
  d.vec <- dmat %>% 
    summarise_all(mean) %>%
    gather(key="Treatment", value="d")
  
  
  #d.vec <- colMeans(dmat)
  
  betamat <- do.call(rbind, nma$samples) %>% data.frame() %>% select(starts_with("beta."))
  trt.names <- nma$trt.key
  colnames(betamat) <- trt.names
  
  betamat %<>% select(-matches(nma$model$reference))
  
  beta.vec <- betamat %>% 
    summarise_all(mean) %>% 
    gather(key="Treatment", value="beta")
  
  #beta.vec <- colMeans(betamat) %>% as.data.frame()
  #colnames(beta.vec) <- c("Treatment", "beta")
  
  xvalues <- seq(x.range[1],x.range[2],length.out=1000) %>% as.data.frame
  colnames(xvalues) <- "x"
  
  summary_values <- full_join(d.vec, beta.vec, by="Treatment")
  
  plotting_values <- left_join(summary_values %>% mutate(one=1), 
                               xvalues %>% mutate(one=1),
                               by="one") %>%
    mutate(y= d+beta*(x-nma$model$mean.cov))
  
  ##Extend color palette if necessary
  n.trts <- dim(d.vec)[1]
  max.colors <- brewer.pal.info[palette,]$maxcolors
  
  if (max.colors < n.trts){
    tmp.colors <- brewer.pal(max.colors, palette)
    plot.colors <- colorRampPalette(tmp.colors)(n.trts)
  } else{
    plot.colors <- brewer.pal(n.trts, palette)
  }
  
  g <- ggplot(plotting_values, aes(x=x, y=y, group=Treatment)) +
    geom_line(aes(color=Treatment), size=lwd) +
    scale_color_manual(values = plot.colors)+
    theme_bw() +
    labs(x=paste0("Values of the covariate ", nma$model$covariate), 
         y=paste0("Effect relative to ",
                  nma$model$reference, 
                  " on the linear scale"),
         color="Treatment")
    
  return("regplot"=g)
  
}