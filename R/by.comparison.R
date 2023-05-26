#' @importFrom dplyr left_join mutate rename select
#' @importFrom magrittr %>%
#' @importFrom purrr map
#' @importFrom tibble as_tibble
#' @importFrom tidyr nest unnest
#' @importFrom utils combn
#' @importFrom magrittr %<>%
#' @importFrom tidyselect all_of


#' @noRd
by.comparison <- function(
  data.nma,
  outcome,
  type.outcome = "binomial",
  N,
  sd = NULL,
  time = NULL
){
  
  #Binding variables to function
  trial <- NULL
  trt <- NULL
  treatments <- NULL
  comparisons <- NULL
  V1 <- NULL
  V2 <- NULL
  trt.e <- NULL
  trt.c <- NULL
  
  data <- data.nma$arm.data %>% select(all_of(outcome), data.nma$varname.t, data.nma$varname.s, all_of(N), all_of(sd), all_of(time))
  names(data)[names(data) == data.nma$varname.t] <- "trt"
  names(data)[names(data) == data.nma$varname.s] <- "trial"
  
  if (type.outcome=="continuous"){

    names(data)[names(data) == data.nma$sd] <- "sd"
    
    data %<>% select(trial, trt, all_of(outcome), all_of(N), all_of(sd))
    data.st <- select(data, trial, trt)
    data.st %<>% nest(treatments=c(trt))
    data.st %<>%
      mutate(comparisons = map(treatments, function(x) combn(x[[1]], 2) %>% t() %>% `colnames<-`(c("V1", "V2")) %>% as_tibble)) %>%
      select(-treatments) %>%
      unnest(cols = c(comparisons)) %>%
      rename(trt.e=V1, trt.c=V2) %>%
      left_join(data %>% select(trial, all_of(outcome), all_of(N), trt, all_of(sd)), by = c("trial", "trt.e" = "trt")) %>%
      left_join(data %>% select(trial, all_of(outcome), all_of(N), trt, all_of(sd)), by = c("trial", "trt.c" = "trt"), suffix=c(".e",".c"))%>%
      mutate(comparison = ifelse(trt.e < trt.c, 
                                 paste(trt.e, trt.c, sep = " vs. "),
                                 paste(trt.c, trt.e, sep = " vs. ")))
    
  } else if (type.outcome=="binomial"){
    
    data %<>% select(trial, trt, all_of(outcome), all_of(N))
    data.st <- select(data, trial, trt) 
    data.st %<>% nest(treatments=c(trt))
    data.st %<>% 
      mutate(comparisons = map(treatments, function(x) combn(x[[1]], 2) %>% t() %>% `colnames<-`(c("V1", "V2")) %>% as_tibble)) %>%
      select(-treatments) %>% 
      unnest(cols = c(comparisons)) %>%
      rename(trt.e=V1, trt.c=V2) %>%
      left_join(data %>% select(trial, all_of(outcome), all_of(N), trt), by = c("trial", "trt.e" = "trt")) %>%
      left_join(data %>% select(trial, all_of(outcome), all_of(N), trt), by = c("trial", "trt.c" = "trt"), suffix=c(".e",".c"))%>%
      mutate(comparison = ifelse(trt.e < trt.c, 
                                 paste(trt.e, trt.c, sep = " vs. "),
                                 paste(trt.c, trt.e, sep = " vs. "))) 
    
  } else if (type.outcome %in% c("rate", "rate2")){
    
    names(data)[names(data) == data.nma$time] <- "time"
    
    data %<>% select(trial, trt, all_of(outcome), all_of(N), all_of(time))
    data.st <- select(data, trial, trt)
    data.st %<>% nest(treatments=c(trt))
    data.st %<>%
      mutate(comparisons = map(treatments, function(x) combn(x[[1]], 2) %>% t() %>% `colnames<-`(c("V1", "V2")) %>% as_tibble)) %>%
      select(-treatments) %>%
      unnest(cols = c(comparisons)) %>%
      rename(trt.e=V1, trt.c=V2) %>%
      left_join(data %>% select(trial, all_of(outcome), all_of(N), trt, all_of(time)), by = c("trial", "trt.e" = "trt")) %>%
      left_join(data %>% select(trial, all_of(outcome), all_of(N), trt, all_of(time)), by = c("trial", "trt.c" = "trt"), suffix=c(".e",".c"))%>%
      mutate(comparison = ifelse(trt.e < trt.c, 
                                 paste(trt.e, trt.c, sep = " vs. "),
                                 paste(trt.c, trt.e, sep = " vs. "))) 
  }
  return(data.st)
}
