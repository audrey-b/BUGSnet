# Forest plot
# @description Produces traceplots of the MCMC chains obtained from \code{nma.run()}
# @param jagsoutput An output produced by \code{nma.run()}
# @param comparator The treatment to use as a comparator
# @param central.tdcy The statistic that you want to use in order to measure relative effectiveness. The options are "mean" and "median".
# @param line.size 
# @param x.trans Optional. A string indicating a transformation to apply to the x-axis. Setting this parameter to "log" is useful when there are extreme values. It also allows an easier interpretation of odds ratios and relative ratios because if e.g. treatment B is twice as far from the line y=1 then treatment A then it's OR/RR is twice the one of treatment A. 

nma.forest <- function(jagsoutput, 
                           comparator, 
                           central.tdcy = "median", 
                           line.size=1,
                           x.trans=NULL) {
  
  x2 <- do.call(rbind, jagsoutput$samples) %>% data.frame() %>% select(starts_with("d."))
  trt.names <- jagsoutput$trt.key
  colnames(x2) <- trt.names
  
  x3 <- x2
  new.vars <- paste0(colnames(x2), "-", comparator)
  for(i in 1:ncol(x2)) {
    x3[[new.vars[i]]] <- x2[, i] - x2[,comparator]
  }
  
  x3 %<>% select(new.vars)
  colnames(x3) <- trt.names
  x3 %<>% select(-comparator)

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
  
  tmp1 <- left_join(tmp.mean, tmp.lci, by = "trt") %>%
    left_join(., tmp.uci, by = "trt") %>% data.frame()

f.plot <- ggplot(tmp1, aes(x=trt, y=mean, ymin=lci, ymax=uci)) +
     geom_pointrange(size=line.size) +
     geom_hline(yintercept=1,lty=2) +
     scale_x_discrete(limits = sort(tmp1$trt, decreasing=TRUE)) +
     xlab("Treatment") +
     ylab(paste0(jagsoutput$scale, " relative to ",comparator,
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

return(f.plot)

}

