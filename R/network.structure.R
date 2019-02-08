network.structure <- function(data.nma) {

  trial <- quo(!! as.name(data.nma$varname.s))

  if ("n" %in% colnames(data.nma$raw.data)) {
    nodes <- data.nma$raw.data %>% select(-n) %>% count_(data.nma$varname.t) %>% rename(node.weight = n) %>%
      mutate(id = as.character(1:n())) %>% select(id, data.nma$varname.t, node.weight)
  }
    else {
    nodes <- data.nma$raw.data %>% count_(data.nma$varname.t) %>% rename(node.weight = n) %>%
    mutate(id = as.character(1:n())) %>% select(id, data.nma$varname.t, node.weight)
  }
  
  tmp1 <- data.nma$raw.data %>% 
    left_join(., nodes, by=data.nma$varname.t) %>%
    mutate(id2 = id) %>%
    group_by_(data.nma$varname.s) %>%
    expand(id, id2) %>%
    filter(id != id2) %>% 
    mutate(comparison = ifelse(id < id2,
                               paste(id, id2, sep = "@@@"), 
                               paste(id2, id, sep = "@@@"))) %>%
    distinct(comparison, .keep_all=TRUE) 
 
 study.label <- tmp1 %>% 
     ungroup() %>%
     group_by(comparison) %>%
     summarise(study = paste(!! trial, collapse =", \n"))
  
  edges <- tmp1 %>%
    ungroup() %>% 
    count(comparison) %>%
    rename(edge.weight = n) %>%
    left_join(.,study.label, by="comparison") %>%
    separate(comparison, into=c("from", "to"), sep = "@@@") %>%
    select(from, to, edge.weight, study)
 
  return(list("edges"=edges, "nodes"=nodes))
   
}