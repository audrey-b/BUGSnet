pairwise.all <- function(slr, 
                          outcome,
                          N,
                          method.tau="DL",
                          sm) {
  
  data <- slr$arm.data

  names(data)[names(data) == slr$varname.t] <- "trt"
  names(data)[names(data) == slr$varname.s] <- "trial"
  names(data)[names(data) == slr$outcome] <- "outcome"
  names(data)[names(data) == slr$N] <- "N"
  
  #List all possible comparisons per study
  
  data %<>% select(trial, trt, all_of(outcome), all_of(N))
  data.st <- select(data, trial, trt) 
  data.st %<>% nest(treatments=c(trt))
  data.st %<>% 
    mutate(comparisons = map(treatments, function(x) combn(x[[1]], 2) %>% t() %>% as_tibble)) %>%
    select(-treatments) %>% 
    unnest() %>%
    rename(treatment=V1, comparator=V2) %>%
    left_join(data %>% select(trial, all_of(outcome), all_of(N), trt), by = c("trial", "treatment" = "trt")) %>%
    left_join(data %>% select(trial, all_of(outcome), all_of(N), trt), by = c("trial", "comparator" = "trt"), suffix=c(".t",".c"))

  meta1 <- metabin(data.st[,paste0(outcome,".t")] %>% t() %>% as.vector,
                   data.st[,paste0(N,".t")] %>% t() %>% as.vector,
                   data.st[,paste0(outcome,".c")] %>% t() %>% as.vector,
                   data.st[,paste0(N,".c")] %>% t() %>% as.vector,
                   method = "MH",
                   method.tau=method.tau,
                   sm=sm)
  
  #data.st %<>% bind_cols(effect=meta1$TE, se.effect=meta1$seTE, lower=meta1$lower, upper=meta1$upper, pval=meta1$pval)
  
  #max(data.st$se.effect)
  
 return(data.st)
}

pairwise.all(slr = slr, outcome="r1", N="N", sm="OR")

