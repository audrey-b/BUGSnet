inconsistency.PMA <- function(slr, 
                          jagsoutput,
                          outcome,
                          N,
                          base.trt,
                          model="random", 
                          central.tdcy="median", 
                          line.size=1,
                          x.trans=NULL) {
  
  data <- slr$arm.data

  names(data)[names(data) == slr$varname.t] <- "trt"
  names(data)[names(data) == slr$varname.s] <- "trial"
  
  # subset data where there are pairwise comparisons vs. baseline trt
  tmp.study <- data %>% 
    filter(trt == base.trt) %>% 
    select(trial) %>%
    left_join(., data, by = "trial")
  
  # select trt names minus baseline trt
  trt.names <- tmp.study %>% select(trt) %>% distinct() %>% arrange(trt) %>% pull()
  trt.names <- trt.names[trt.names!=base.trt]
  
  if ("n" %in% colnames(data)) {
    count <- tmp.study %>%
      select(-n) %>%
      filter(trt != base.trt) %>%
      count(trt) %>% arrange(trt) %>% pull(n) 
  }else {  count <- tmp.study %>%
    filter(trt != base.trt) %>%
    count(trt) %>% arrange(trt) %>% pull(n)}

  # empty vectors
  est <- vector(mode="numeric", length=n_distinct(trt.names))
  lci <- vector(mode="numeric", length=n_distinct(trt.names))
  uci <- vector(mode="numeric", length=n_distinct(trt.names))

source("pairwise.R")
  
for(i in 1:n_distinct(trt.names)) {
  pairwise.output <- pairwise(slr, 
                              base.trt, 
                              trt.names[i], 
                              outcome,
                              N,
                              method = "MH",
                              method.tau="DL")
  
if(model=="random") {
est[i] <- pairwise.output$summary$re.estimate
lci[i] <- pairwise.output$summary$re.lci
uci[i] <- pairwise.output$summary$re.uci
}

if(model=="fixed") {
  est[i] <- pairwise.output$summary$fe.estimate
  lci[i] <- pairwise.output$summary$fe.lci
  uci[i] <- pairwise.output$summary$fe.uci
}
}

  # pairwise results
  tmp.pairwise <- data.frame(trt=trt.names,
                             estimate=round((est), digits=2),
                             lower=round((lci), digits=2), 
                             upper=round((uci), digits=2)) %>%
    mutate(trt = as.character(trt), 
           nma = 0)
 
    # extract results from NMA
   source("league.table.R")
  
  tmp.nma <- leaguetable(jagsoutput, central.tdcy="median", layout="default") 
  tmp.nma2 <- data.frame(trt=rownames(tmp.nma), value=tmp.nma[base.trt, ]) %>%
    filter(trt!=base.trt) 
  
  tmp.nma2$value <- gsub(" \\(", "@", tmp.nma2$value)
  tmp.nma2$value <- gsub("\\)", "", tmp.nma2$value)
  tmp.nma2$value <- gsub(" to ", "@", tmp.nma2$value)
  tmp.nma2$value <- gsub("  ", "@", tmp.nma2$value)

  tmp.nma2 <- tmp.nma2 %>%
    separate(., value, into=c("estimate", "lower", "upper"), sep="@") %>%
    mutate(trt = as.character(trt),
           estimate = as.numeric(estimate),
           lower = as.numeric(lower),
           upper = as.numeric(upper),
           nma = 1)
    
    
  # merge datasets  
 tmp.merge <- bind_rows(tmp.nma2, tmp.pairwise) 
 
 p.names <- tmp.merge %>% select(trt) %>% distinct() %>% pull()

 # forest plot
 f.plot <- ggplot(tmp.merge, aes(x=factor(trt), y=estimate, ymin=lower, ymax=upper, color=factor(nma))) +
   geom_pointrange(size=line.size, alpha=0.5, position = position_dodge(width=0.5)) + 
   scale_colour_discrete(labels=c(paste0("Pairwise Meta-Analysis,\n", model,"-effects model,\npooled estimate with 95% CI"),
                                  paste0("Network Meta-Analysis,\nposterior ", central.tdcy," with 95% CrI"))) +
   coord_flip() +
   geom_hline(yintercept=1,lty=2) +
   xlab("Treatment") +
   labs(color="") +
   theme_classic() +
   theme(legend.position = "bottom") +
   scale_x_discrete(limits = sort(p.names, decreasing=TRUE))
 #+
 #labs(caption = paste("note: each treatment compared to", base.trt))
 
 
 if(jagsoutput$scale == "RR"){
   f.plot <- f.plot + ylab(paste0("RR relative to ",base.trt))
 } else if(jagsoutput$scale == "OR"){
   f.plot <- f.plot + ylab(paste0("OR relative to ",base.trt))
 }
 
 if(is.null(x.trans)){
   f.plot <- f.plot + 
     scale_y_continuous(breaks = scales::pretty_breaks(c(tmp.merge$lower,tmp.merge$upper), n = 10))
 }else{
     f.plot <- f.plot + 
     scale_y_continuous(trans=x.trans,
                        breaks = scales::pretty_breaks(c(min(tmp.merge$lower),max(tmp.merge$upper)), n = 10))
 }
 
 tmp.nma2 %<>% transmute(Treatment=trt, 
                         Comparator=base.trt,
                         NMA=paste(estimate, 
                                              " (", 
                                              lower,
                                              " to ", 
                                            upper, 
                                              ")", sep=""))
 
 tmp.pairwise %<>% transmute(Treatment=trt, 
                             Comparator=base.trt,
                             Pairwise=paste(estimate, 
                                              " (", 
                                              lower,
                                              " to ", 
                                              upper, 
                                              ")", sep=""))

 comp.table <- full_join(tmp.nma2, tmp.pairwise, by=c("Treatment","Comparator"))
 
 return(list(table=comp.table, plot=f.plot))
  
}


  
  
