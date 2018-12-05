by.comparison <- function(slr, outcome, type.outcome="binomial", N, sd=NULL, time = NULL){
  data <- slr$raw.data
  names(data)[names(data) == slr$varname.t] <- "trt"
  names(data)[names(data) == slr$varname.s] <- "trial"
  #names(data)[names(data) == outcome] <- "outcome"
  #names(data)[names(data) == N] <- "N"
  
  if (type.outcome=="continuous"){

    names(data)[names(data) == slr$sd] <- "sd"
    
    #List all possible comparisons per study
    
    data %<>% select(trial, trt, outcome, N, sd)
    data.st <- select(data, trial, trt)
    data.st %<>% nest(trt, .key="treatments")
    data.st %<>%
      mutate(comparisons = map(treatments, function(x) combn(x[[1]], 2) %>% t() %>% as_tibble)) %>%
      select(-treatments) %>%
      unnest() %>%
      rename(trt.e=V1, trt.c=V2) %>%
      left_join(data %>% select(trial, outcome, N, trt, sd), by = c("trial", "trt.e" = "trt")) %>%
      left_join(data %>% select(trial, outcome, N, trt, sd), by = c("trial", "trt.c" = "trt"), suffix=c(".e",".c"))%>%
      mutate(comparison = ifelse(trt.e < trt.c, 
                                 paste(trt.e, trt.c, sep = " vs. "),
                                 paste(trt.c, trt.e, sep = " vs. ")))
    
  } else if (type.outcome=="binomial"){
    
    data %<>% select(trial, trt, outcome, N)
    data.st <- select(data, trial, trt) 
    data.st %<>% nest(trt, .key="treatments")
    data.st %<>% 
      mutate(comparisons = map(treatments, function(x) combn(x[[1]], 2) %>% t() %>% as_tibble)) %>%
      select(-treatments) %>% 
      unnest() %>%
      rename(trt.e=V1, trt.c=V2) %>%
      left_join(data %>% select(trial, outcome, N, trt), by = c("trial", "trt.e" = "trt")) %>%
      left_join(data %>% select(trial, outcome, N, trt), by = c("trial", "trt.c" = "trt"), suffix=c(".e",".c"))%>%
      mutate(comparison = ifelse(trt.e < trt.c, 
                                 paste(trt.e, trt.c, sep = " vs. "),
                                 paste(trt.c, trt.e, sep = " vs. "))) 
    
  } else if (type.outcome %in% c("rate", "rate2")){
    
    names(data)[names(data) == slr$time] <- "time"
    
    #List all possible comparisons per study
    
    data %<>% select(trial, trt, outcome, N, time)
    data.st <- select(data, trial, trt)
    data.st %<>% nest(trt, .key="treatments")
    data.st %<>%
      mutate(comparisons = map(treatments, function(x) combn(x[[1]], 2) %>% t() %>% as_tibble)) %>%
      select(-treatments) %>%
      unnest() %>%
      rename(trt.e=V1, trt.c=V2) %>%
      left_join(data %>% select(trial, outcome, N, trt, time), by = c("trial", "trt.e" = "trt")) %>%
      left_join(data %>% select(trial, outcome, N, trt, time), by = c("trial", "trt.c" = "trt"), suffix=c(".e",".c"))%>%
      mutate(comparison = ifelse(trt.e < trt.c, 
                                 paste(trt.e, trt.c, sep = " vs. "),
                                 paste(trt.c, trt.e, sep = " vs. "))) 
  }
  return(data.st)
}
