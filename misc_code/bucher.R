###Non working code
#This code identifies loops of evidence in the data



# Bucher Method (2 arm studies only) --------------------------------------

#data.outcome.rct.2arms <- data.outcome.rct %>%
#  group_by(trial) %>%
#  mutate(n.arms.study = n()) %>%
#  ungroup() %>%
#  filter(n.arms.study == 2)

source("R\\network.plot.R")

#network.plot(data.outcome.rct.2arms, node.scale=5, edge.scale=2)

#data.outcome.rct.2arms <- data.outcome.rct %>%
#  filter(arm %in% c(1,2))

#network.plot(data.outcome.rct.2arms, node.scale=5, edge.scale=2)

source("R\\network.structure.R")

#edgesANDnodes <- network.structure(data.outcome.rct.2arms)
#edges <- edgesANDnodes[[1]]
#nodes <- edgesANDnodes[[2]]

#net <- graph_from_data_frame(d=edges, vertices=nodes, directed=F) 


edgesANDnodes <- network.structure(data.outcome.rct) #need to use slr
edges <- edgesANDnodes[[1]]
nodes <- edgesANDnodes[[2]]

net <- graph_from_data_frame(d=edges, vertices=nodes, directed=F) 
network.plot(data.outcome.rct, node.scale=5, edge.scale=2)




#start with all direct pairwise comparisons
#find all loops for those with all_simple_paths(net, from=, to=)
#sort
#unique
#now for each loop, check if it is only 2 arm studies, if not, drop loop

network <- data.outcome.rct

edgesANDnodes <- network.structure(network) #need to use slr
edges <- edgesANDnodes[[1]]
nodes <- edgesANDnodes[[2]]

net <- graph_from_data_frame(d=edges, vertices=nodes, directed=F) 

tmp1 <- network %>% 
  select(trial, trt) %>% 
  nest(trt)

pairs <- tmp1[1,"data"] %>% unlist %>% sort %>% combn(2)
for(i in 2:nrow(tmp1)){
  pairs <- tmp1[i,"data"] %>% unlist %>% sort %>% combn(2) %>% cbind(pairs)
}
direct.comparisons <- unique(pairs, MARGIN=2)

loops <- NULL

for(i in 1:ncol(direct.comparisons)){
  loops %<>% c(all_simple_paths(net, 
                                from=which(V(net)$trt==direct.comparisons[1,i]), 
                                to=which(V(net)$trt==direct.comparisons[2,i])))
}

loops %<>% lapply(function(x){names(x) %>% as.integer})


#loops.tibble <- tibble(loops) %>%
#  mutate(sorted.nodes = map(loops,sort))

#loops.tibble$sorted.nodes %>% duplicated(fromLast = TRUE)

#distinct(loops.tibble, sorted.nodes, .keep_all=TRUE)

#loops %<>% lapply(function(x){names(x) %>% as.integer %>% sort}) %>%
#  unique

is.loop <- lapply(loops, function(x){length(x)>=3}) %>% unlist
is.duplicate <- lapply(loops, sort) %>% duplicated

loops <- loops[which(is.loop & !is.duplicate)]


test <- loops[[1]]


for(i in 1:(length(test)-1)){
  print(test[i:(i+1)] %>% sort)
}


edges %<>%
  mutate(fromto = paste0(from, "--", to))

edges %>%
  filter(fromto %in% c("1--9","1--5")) %>%
  select(study)

#V(net)$trt[test] %>% combn(2)



loops.unique.studies <- 
  
  network.plot(network)

#ICDF <- n.contrasts - n.treatments + 1
