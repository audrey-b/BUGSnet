inconsistency.model <- function(inconsistency.out, nma.out, central.tdcy="median"){
  
  source("league.table.R")
  
  
  
  summary.incons <- summary(inconsistency.out$samples)
  
  
  
  
  
  tmp <- cbind(summary.incons[[2]][,"50%"] %>% exp(), #need to choose median/mean
               
               #summary.incons[[1]][,"SD"],  #need to change all those to exponential
               
               summary.incons[[2]][,"2.5%"] %>% exp(),
               
               summary.incons[[2]][,"97.5%"] %>% exp()) %>%
    
    as.data.frame() %>%
    
    tibble::rownames_to_column() %>%
    
    mutate(rowname=map_chr(rowname, function(x){paste0(inconsistency.out$trt.key[strsplit(x,"\\[|\\]|,")[[1]][3] %>% as.integer],
                                                       
                                                       " - ",
                                                       
                                                       inconsistency.out$trt.key[strsplit(x,"\\[|\\]|,")[[1]][2] %>% as.integer])}))
  
  #mutate(Treatment=map(rowname, function(x){inconsistency.out$trt.key[strsplit(x,"\\[|\\]|,")[[1]][2] %>% as.integer]}), 
  
  #       Comparator=map(rowname, function(x){inconsistency.out$trt.key[strsplit(x,"\\[|\\]|,")[[1]][3] %>% as.integer]}))
  
  #colnames(tmp) <- c("rowname", "mean", "SD", "lower", "upper","Treatment", "Comparator")
  
  #tmp %<>% select(Treatment, Comparator, mean, SD, lower, upper)
  
  colnames(tmp) <- c("Comparison", "median", "lci", "uci")
  
  tmp %<>% filter(uci < exp(10)) %>% ###ARBITRARY
    
    mutate(method="inconsistency model")
  
  
  
  binded <- leaguetable(nma.out, central.tdcy=central.tdcy, layout="long") %>%
    
    mutate(method="nma", Comparison = paste0(Treatment, " - ", Comparator)) %>%
    
    filter(Comparison %in% tmp$Comparison) %>%
    
    bind_rows(tmp)
  
  
  
  
  
  
  
  #Devon's code (adapted)
  
  #####Which is which???? factor/color
  
  # forest plot
  
  line.size=1
  
  f.plot <- ggplot(binded, aes(x=factor(Comparison), y=median, ymin=lci, ymax=uci, color=factor(method))) +
    
    geom_pointrange(size=0.4, alpha=0.5, position = position_dodge(width=0.5)) + 
    
    scale_colour_discrete(labels=c(paste0("Inconsistency Model,\nposterior ",central.tdcy," with 95% CrI"),
                                   
                                   paste0("Consistency Model,\nposterior ",central.tdcy," with 95% CrI"))) +
    
    coord_flip() +
    
    geom_hline(yintercept=1,lty=2) +
    
    xlab("Treatment Comparison") +
    
    labs(color="") +
    
    theme_classic() +
    
    theme(legend.position = "bottom") #+
  
  #scale_x_discrete(limits = sort(p.names, decreasing=TRUE))
  
  
  
  if(inconsistency.out$scale == "Risk Ratio" & nma.out$scale == "Risk Ratio"){
    
    f.plot <- f.plot + ylab("Risk Ratio")
    
  } else if(inconsistency.out$scale == "Odds Ratio" & nma.out$scale == "Odds Ratio"){
    
    f.plot <- f.plot + ylab("Odds Ratio")
    
  }
  
  
  
  return(f.plot)
  
  
}