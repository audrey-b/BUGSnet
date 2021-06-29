network.structure <- function(data.nma) {
  
  # Bind Variables to function
  node.weight <- NULL
  from <- NULL
  to <- NULL
  trt <- NULL
  flag <- NULL
  mtchvar <- NULL
  
  
  trial <- rlang::quo(!! as.name(data.nma$varname.s))
  varname.t.quo <- rlang::quo(!! as.name(data.nma$varname.t))
  
  if ("n" %in% colnames(data.nma$arm.data)) {
    nodes <- data.nma$arm.data %>% select(-n) %>% count(!! varname.t.quo) %>% rename(node.weight = n) %>%
      mutate(id = as.character(1:n())) %>% select(data.nma$varname.t, node.weight)
  }else {
    nodes <- data.nma$arm.data %>% count(!! varname.t.quo) %>% rename(node.weight = n) %>%
      mutate(id = as.character(1:n())) %>% select(data.nma$varname.t, node.weight)
  }
  
  studytrt <- data.nma$arm.data %>%
    select(data.nma$varname.s, data.nma$varname.t) %>%
    nest(data=c(data.nma$varname.t))
  
  cnt <- data.nma$arm.data %>%
    select(data.nma$varname.s, data.nma$varname.t) %>%
    count(across(data.nma$varname.s))
  tmp1 <- bind_cols(studytrt, cnt) %>%
    filter(n>1)
  
  pairs <- tmp1[1,"data"] %>% unlist %>% sort %>% combn(2)
  
  for(i in 2:nrow(tmp1)){
    pairs <- tmp1[i,"data"] %>% unlist %>% sort %>% combn(2) %>% cbind(pairs)
  }
  
  pairs2 <- data.frame(from = pairs[1,],
                       to = pairs[2,]) %>%
    group_by(from, to) %>%
    mutate(edge.weight = max(1:n())) %>%
    arrange(from, to) %>%
    distinct() %>%
    mutate(mtchvar = 1)
  studylabs <- studytrt %>%
    group_by(!! trial) %>%
    mutate(trt = paste( unlist(data), collapse=';')) %>%
    select(!! trial, trt) %>%
    mutate(mtchvar = 1)
  edges <- left_join(pairs2, studylabs, by="mtchvar") %>%
    ungroup() %>%
    mutate(trt = as.character(trt),
           from = as.character(from),
           to = as.character(to)) %>%
    mutate(flag = ifelse(stringr::str_detect(trt, stringr::coll(from)) &
                           stringr::str_detect(trt, stringr::coll(to)), 1, 0)) %>%
    filter(flag == 1) %>%
    select(-c(mtchvar, flag, trt)) %>%
    nest(data=c(!! trial)) %>%
    group_by(from, to) %>%
    mutate(study = paste(unlist(data), collapse=', \n')) %>%
    select(-data)
  
  return(list("edges"=edges, "nodes"=nodes))
}
