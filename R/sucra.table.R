#' SUCRA Table and Plot
#' @description Produces a SUCRA (Surface Under the Cumulative Ranking Curve) plot and table. A Sucra table summarizes the probabilities
#' that each  treatment is the best, second best...worst treatment in the network.
#' 
#' @param jagsoutput Resulting dataset from running \code{nma.analysis()}.
#' @param largerbetter A boolean variable indicating whether a larger probability should indicate a more effective treatment (TRUE) or
#' if a smaller probability should indicate a more effective treatment (FALSE). 
#' @param line.thickness A number indicating the thickness of the lines in the SUCRA plot.
#' @param colour.set A string indicating the colour set from RcolorBrewer. "set1" is great, but may need a different colour set if 
#' there are many treatments in your network.
#' 
#' @return \code{s.table} - A dataset containing a SUCRA table.
#' @return \code{s.plot} - A plot showing the probability of each treatment being the nth best treatment.
#' 
#' @examples
#' 
#' #get sucra results
#' sucra_results <- sucra(jagsoutput = nma_results, largerbetter = TRUE)
#' 
#' #plot sucra results
#' sucra_results$s.plot


sucra <- function(jagsoutput, 
                  largerbetter, 
                  line.thickness = 1.0,
                  colour.set= "set1") {
  
  x <- do.call(rbind, jagsoutput$samples) %>% data.frame() %>% select(starts_with("d."))
  
  tmp.var <- vector(mode="character", length=ncol(x) - 1)
  for(i in 1:ncol(x)) {
    tmp.var[i] <- i
  }
  
colnames(x) <- jagsoutput$trt.key
  
  
x2 <- x
x3 <- x2 %>% 
  mutate(iteration = row_number()) %>% 
  gather(trt, value, -iteration) %>%
  group_by(iteration)

if(largerbetter){
  x3 %<>% arrange(iteration, desc(value))
  }else if(!largerbetter){
    x3 %<>% arrange(iteration, value)
  } 

x3 %<>% summarise(tmp.rank = paste(trt, collapse =",")) %>%
  select(-iteration) %>%
  separate_(., col="tmp.rank", into=tmp.var, sep=",") %>%
  gather(key=rank, value=trt) %>%
  mutate(rank=as.numeric(rank))


s.table <- x3 %>% 
  group_by(rank) %>% 
  count(trt) %>%
  mutate(percent = round(n / nrow(x)*100, 2)) %>%
  select(-n) %>%
  spread(trt, -rank) %>%
  mutate_all(funs(replace(., is.na(.), 0)))

x4 <- s.table %>% gather(trt, prob, -rank)

s.plot <- ggplot(data=x4, aes(x=factor(rank), y=prob, group=trt)) +
  geom_line(aes(color=trt), size=line.thickness) +
  geom_point(aes(color=trt)) +
  scale_color_manual(values = brewer.pal(n = length(jagsoutput$trt.key), name = colour.set))+
  theme_bw()


if(largerbetter==TRUE){
  s.plot <- s.plot +
    labs(x=paste0("Rank of Treatment",
                  "\n(Higher ranks associated with larger outcome values)"), 
         y="Probability (%)",
         color="Treatment")
    
}else if(largerbetter==FALSE){
  s.plot <- s.plot +
    labs(x=paste0("Rank of Treatment",
                  "\n(Higher ranks associated with smaller outcome values)"), 
         y="Probability (%)",
         color="Treatment")
}

x5 <- (list(s.table, s.plot))
names(x5) <- c("s.table", "s.plot")
return(x5)
}


