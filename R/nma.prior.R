# Calculates estimate of max-delta to be used in prior calculation

nma.prior <- function(data_arm, data_contrast, outcome, differences, scale, N, sd=NULL, time = NULL){
  
  if(!is.null(data_arm)){
  
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

  # build table of treatment comparisons for which we have data (for arm-based)
  table <- by.comparison(data_arm, outcome, type.outcome = type.outcome, N, sd=sd, time = time)
  names(table)[names(table) == paste0(outcome,".e")] <- "outcome.e"
  names(table)[names(table) == paste0(outcome,".c")] <- "outcome.c"
  
  
  #calculate MLEs of deltas from arm-based data table
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
    # make numeric
    deltas <- pull(deltas, delta)
    
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
    # make numeric
    deltas <- pull(deltas, delta)
    
  } else if (scale == "Mean Difference"){
    
    deltas <- table %>% 
      mutate(delta = as.numeric(outcome.e)-as.numeric(outcome.c)) %>%
      select(delta)
    # make numeric
    deltas <- pull(deltas, delta)
    
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
    # make numeric
    deltas <- pull(deltas, delta)
    
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
    # make numeric
    deltas <- pull(deltas, delta)
    
  }
  } else {deltas <- 0}
  
  # calculate deltas from contrast-based data
  
  if(!is.null(data_contrast)) {
    
    deltas2 <- data_contrast$arm.data %>% select(differences) # the differences are already reported, we just need to select the outcome column
    deltas2 <- pull(deltas2, differences)
    deltas2 <- as.numeric(deltas2[!(deltas2%in%c(NA,"NA"))]) # drop NAs

  } else {deltas2 <- 0}

  # return maximum delta for priors
  return(max(abs(c(deltas,deltas2))))

}



  