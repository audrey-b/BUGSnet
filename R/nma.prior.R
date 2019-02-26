nma.prior <- function(data.nma, outcome, scale, N, sd=NULL, time = NULL){
  if (scale =="OR" ){type.outcome = "binomial"}
  else if (scale =="RR"){type.outcome = "binomial"}
  else if (scale =="MD"){type.outcome = "continuous"}
  else if (scale =="HR"){type.outcome = "rate"}
  else if (scale =="Rate Ratio"){type.outcome = "rate2"}
  
  table <- by.comparison(data.nma, outcome, type.outcome = type.outcome, N, sd=sd, time = time)
  names(table)[names(table) == paste0(outcome,".e")] <- "outcome.e"
  names(table)[names(table) == paste0(outcome,".c")] <- "outcome.c"
  
  if (scale == "OR"){
    names(table)[names(table) == paste0(N,".e")] <- "N.e"
    names(table)[names(table) == paste0(N,".c")] <- "N.c"
    
    deltas <- table %>% 
      mutate(adj_r.e = (outcome.e+0.5)/(N.e + 1), # add 0.5 to ensure ratio is non-zero
             adj_r.c = (outcome.c+0.5)/(N.c + 1)) %>%
      mutate(theta.e = log(adj_r.e/(1-adj_r.e)),
             theta.c = log(adj_r.c/(1-adj_r.c))) %>%
      mutate(delta = theta.e-theta.c) %>%
      select(delta)
    
  } else if (scale =="RR"){
    names(table)[names(table) == paste0(N,".e")] <- "N.e"
    names(table)[names(table) == paste0(N,".c")] <- "N.c"
    
    deltas <- table %>% 
      mutate(adj_r.e = (outcome.e+0.5)/(N.e + 1),
             adj_r.c = (outcome.c+0.5)/(N.c + 1)) %>%
      mutate(theta.e = log(adj_r.e),
             theta.c = log(adj_r.c)) %>%
      mutate(delta = theta.e-theta.c) %>%
      select(delta)
    
  } else if (scale == "MD"){
    
    deltas <- table %>% 
      mutate(delta = as.numeric(outcome.e)-as.numeric(outcome.c)) %>%
      select(delta)
    
  } else if (scale == "Rate Ratio"){
    names(table)[names(table) == paste0(N,".e")] <- "N.e"
    names(table)[names(table) == paste0(N,".c")] <- "N.c"
    names(table)[names(table) == paste0(time,".e")] <- "time.e"
    names(table)[names(table) == paste0(time,".c")] <- "time.c"
    
    deltas <- table %>% 
      mutate(adj_r.e = (outcome.e+0.5)/((N.e+1)*time.e),
             adj_r.c = (outcome.c+0.5)/((N.c+1)*time.c)) %>%
      mutate(theta.e = log(-log(adj_r.e)),
             theta.c = log(-log(adj_r.c))) %>%
      mutate(delta = theta.e-theta.c) %>%
      select(delta)
    
  } else if (scale =="HR"){
    names(table)[names(table) == paste0(time,".e")] <- "time.e"
    names(table)[names(table) == paste0(time,".c")] <- "time.c"
    
    deltas <- table %>% 
      mutate(adj_r.e = (outcome.e+0.5)/(time.e+1),
             adj_r.c = (outcome.c+0.5)/(time.c+1)) %>%
      mutate(theta.e = log(adj_r.e),
             theta.c = log(adj_r.c)) %>%
      mutate(delta = theta.e-theta.c) %>%
      select(delta)
    
  }
  
  return(max(abs(deltas)))
}

#nma.prior(data.nma = dich.slr, scale = "RR", outcome = "responders", N = "sampleSize")



  