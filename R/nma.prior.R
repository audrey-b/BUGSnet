nma.prior <- function(data.nma, outcome, scale, N, sd=NULL, time = NULL){
  
  #Bind variables to function
  outcome.e <- NULL
  N.e <- NULL
  outcome.c <- NULL
  N.c <- NULL
  adj_r.e <- NULL
  adj_r.c <- NULL
  theta.e <- NULL
  theta.c <- NULL
  delta <- NULL
  time.e <- NULL
  time.c <- NULL
  
  
  if (scale =="Odds Ratio" ){
    type.outcome = "binomial"
  } else if (scale =="Risk Ratio"){
    type.outcome = "binomial"
  } else if (scale =="Mean Difference"){
    type.outcome = "continuous"
  } else if (scale =="Rate Ratio"){
    type.outcome = "rate"
  } else if (scale =="Hazard Ratio"){
    type.outcome = "rate2"
  }
  
  table <- by.comparison(data.nma, outcome, type.outcome = type.outcome, N, sd=sd, time = time)
  names(table)[names(table) == paste0(outcome,".e")] <- "outcome.e"
  names(table)[names(table) == paste0(outcome,".c")] <- "outcome.c"
  
  if (scale == "Odds Ratio"){
    names(table)[names(table) == paste0(N,".e")] <- "N.e"
    names(table)[names(table) == paste0(N,".c")] <- "N.c"
    
    deltas <- table %>% 
      mutate(adj_r.e = (outcome.e+0.5)/(N.e + 1), # add 0.5 to ensure ratio is non-zero
             adj_r.c = (outcome.c+0.5)/(N.c + 1)) %>%
      mutate(theta.e = log(adj_r.e/(1-adj_r.e)),
             theta.c = log(adj_r.c/(1-adj_r.c))) %>%
      mutate(delta = theta.e-theta.c) %>%
      select(delta)
    
  } else if (scale =="Risk Ratio"){
    names(table)[names(table) == paste0(N,".e")] <- "N.e"
    names(table)[names(table) == paste0(N,".c")] <- "N.c"
    
    deltas <- table %>% 
      mutate(adj_r.e = (outcome.e+0.5)/(N.e + 1),
             adj_r.c = (outcome.c+0.5)/(N.c + 1)) %>%
      mutate(theta.e = log(adj_r.e),
             theta.c = log(adj_r.c)) %>%
      mutate(delta = theta.e-theta.c) %>%
      select(delta)
    
  } else if (scale == "Mean Difference"){
    
    deltas <- table %>% 
      mutate(delta = as.numeric(outcome.e)-as.numeric(outcome.c)) %>%
      select(delta)
    
  } else if (scale == "Hazard Ratio"){
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
    
  } else if (scale =="Rate Ratio"){
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

#nma.prior(data.nma = dich.slr, scale = "Risk Ratio", outcome = "responders", N = "sampleSize")



  