network.structure <- function(slr) {

  trial <- quo(!! as.name(slr$varname.s))

  if ("n" %in% colnames(slr$raw.data)) {
    nodes <- slr$raw.data %>% select(-n) %>% count_(slr$varname.t) %>% rename(node.weight = n) %>%
      mutate(id = as.character(1:n())) %>% select(id, slr$varname.t, node.weight)
  }
    else {
    nodes <- slr$raw.data %>% count_(slr$varname.t) %>% rename(node.weight = n) %>%
    mutate(id = as.character(1:n())) %>% select(id, slr$varname.t, node.weight)
  }
  
  tmp1 <- slr$raw.data %>% 
    left_join(., nodes, by=slr$varname.t) %>%
    mutate(id2 = id) %>%
    group_by_(slr$varname.s) %>%
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