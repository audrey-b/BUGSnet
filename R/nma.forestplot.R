nma.forest <- function(jagsoutput, 
                           base.trt, 
                           central.tdcy = "median", 
                           line.size=1,
                           x.trans=NULL) {
  
  x2 <- do.call(rbind, jagsoutput$samples) %>% data.frame() %>% select(starts_with("d."))
  trt.names <- jagsoutput$trt.key
  colnames(x2) <- trt.names
  
  x3 <- x2
  new.vars <- paste0(colnames(x2), "-", base.trt)
  for(i in 1:ncol(x2)) {
    x3[[new.vars[i]]] <- x2[, i] - x2[,base.trt]
  }
  
  x3 %<>% select(new.vars)
  colnames(x3) <- trt.names
  x3 %<>% select(base.trt)

  tmp.mean <- x3 %>%  
    summarise_all(funs(mean = e.mean)) %>% gather() %>%
    rename(trt = key, mean = value) %>%
    mutate(trt = sub("_mean", "", trt))
  
  tmp.lci <- x3 %>%  
    summarise_all(funs(lci = e.lci)) %>% gather() %>%
    rename(trt = key, lci = value) %>%
    mutate(trt = sub("_lci", "", trt))
  
  tmp.uci <- x3 %>%  
    summarise_all(funs(uci = e.uci)) %>% gather() %>%
    rename(trt = key, uci = value) %>%
    mutate(trt = sub("_uci", "", trt))
  
  tmp1 <- left_join(tmp.mean, tmp.lci, by = "trt") %>%
    left_join(., tmp.uci, by = "trt") %>% data.frame()

f.plot <- ggplot(tmp1, aes(x=trt, y=mean, ymin=lci, ymax=uci)) +
     geom_pointrange(size=line.size) +
     geom_hline(yintercept=1,lty=2) +
     scale_x_discrete(limits = sort(tmp1$trt, decreasing=TRUE)) +
     xlab("Treatment") +
     ylab(paste0(jagsoutput$scale, " relative to ",base.trt,
                 "\n(showing posterior ", central.tdcy," with 95% CrI)")) +
     coord_flip() +
     theme_classic() #+
     #labs(caption = paste("note: each treatment compared to", base.trt))

if(is.null(x.trans)){
f.plot <- f.plot +
  scale_y_continuous(breaks = scales::pretty_breaks(c(tmp1$lci,tmp1$uci), n = 10))
}
else{

  f.plot <- f.plot +
    scale_y_continuous(trans=x.trans,
                       breaks = scales::pretty_breaks(c(min(tmp1$lci),max(tmp1$uci)), n = 10))
}

return(f.plot)

}

