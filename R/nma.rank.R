#' Table and Plots of Treatment Rankings
#' @description Produces a SUCRA (Surface Under the Cumulative Ranking Curve) plot and table. A Sucra table summarizes the probabilities
#' that each treatment is the best, second best...worst treatment in the network.
#' 
#' @param nma Resulting output from running \code{nma.run()}.
#' @param largerbetter A boolean variable indicating whether a larger probability should indicate a more effective treatment (TRUE) or
#' if a smaller probability should indicate a more effective treatment (FALSE). 
#' @param sucra.lwd Line width relative to the default (default=1) in the SUCRA plot.
#' @param sucra.palette A string indicating the colour set from RcolorBrewer for the SUCRA plot. "set1" is great, but you may need a different colour set if 
#' there are many treatments in your network.
#' 
#' @return \code{ranktable} - A rank table showing the probability of each treatment being the nth best treatment.
#' @return \code{sucratable} - A table showing the probability of each treatment being the nth best treatment or better and an overall SUCRA value for each treatment.
#' @return \code{order} - A vector containing the order of efficacy of treatments (from best to worst) based on their SUCRA value. This vector 
#' is useful for ordering treatments when creating the league heat plot with \code{nma.league()}.
#' @return \code{longtable} - A long form table of ranking probabilities and SUCRA value.
#' @return \code{sucraplot} - A SUCRA plot showing the probability of each treatment being the nth best treatment or better.
#' @return \code{rankogram} - A rankogram showing the probability of each treatment being the nth best treatment.

#' 
#' @examples
#' 
#' #get sucra results
#' sucra_results <- nma.rank(nma = nma_results, largerbetter = TRUE)
#' 
#' #plot sucra results
#' sucra_results$sucra


nma.rank <- function(nma, 
                  largerbetter, 
                  sucra.lwd = 1.0,
                  sucra.palette= "Set1") {
  
  x <- do.call(rbind, nma$samples) %>% data.frame() %>% select(starts_with("d."))
  
  tmp.var <- vector(mode="character", length=ncol(x) - 1)
  for(i in 1:ncol(x)) {
    tmp.var[i] <- i
  }
  
colnames(x) <- nma$trt.key
  
  
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
  replace(., is.na(.), 0)

x4 <- s.table %>% gather(trt, prob, -rank)

x5 <- x4 %>% 
  ungroup() %>% 
  group_by(trt) %>%
  mutate(cumprob = cumsum(prob))

sucra.table <- x5 %>% 
  select(-prob) %>% 
  spread(trt, -rank)
sucra.table$rank = as.character(sucra.table$rank)
sucras <- colMeans(sucra.table[1:(dim(sucra.table)[1]-1),-1]) %>% round(2)
sucra.table %<>% rbind(c("SUCRA",sucras))

order <- sucra.table[dim(sucra.table)[1],] %>% 
  gather(key=rank, value=SUCRA) %>%
  mutate(SUCRA=as.numeric(SUCRA)) %>%
  arrange(by=desc(SUCRA)) %>%
  select(rank) %>%
  t() %>%
  as.vector

long.table <- x5

s.plot <- ggplot(data=x5, aes(x=factor(rank), y=cumprob, group=trt)) +
  geom_line(aes(color=trt), size=sucra.lwd) +
  geom_point(aes(color=trt)) +
  scale_color_manual(values = brewer.pal(n = length(nma$trt.key), name = sucra.palette))+
  theme_bw()

rankogram <- ggplot(data=x4, aes(y=prob, x=trt, fill=factor(rank)))+ 
  geom_bar(stat = "identity", width = 0.5) +
  scale_fill_brewer(palette="Blues")


if(largerbetter==TRUE){
  s.plot <- s.plot +
    labs(x=paste0("Rank of Treatment",
                  "\n(Higher ranks associated with larger outcome values)"), 
         y="Probability of that rank or higher (%)",
         color="Treatment")
  rankogram <- rankogram + labs(x=paste0("Treatment",
                                         "\n(Higher ranks associated with larger outcome values)"), 
                                y="Probability of rank (%)",
                                fill="Rank")
    
}else if(largerbetter==FALSE){
  s.plot <- s.plot +
    labs(x=paste0("Rank of Treatment",
                  "\n(Higher ranks associated with smaller outcome values)"), 
         y="Probability of that rank or higher (%)",
         color="Treatment")
  rankogram <- rankogram + labs(x=paste0("Treatment",
                                         "\n(Higher ranks associated with smaller outcome values)"), 
                                y="Probability of rank (%)",
                                fill="Rank")
}

#output sucra ranks to input into league table
# sucra.totals <- s.table %>%ungroup() %>% select(-rank) %>% t()
# sucra.ranks <- sucra.totals %*% c(1:ncol(sucra.totals)) %>% 
#   data.frame() %>% 
#   cbind(rownames(sucra.totals))
# 
# colnames(sucra.ranks) <- c("total", "treatment")
# 
# sucra.ranks <- sucra.ranks %>% arrange(total) %>% select(treatment)

x5 <- (list(s.table, sucra.table, order, long.table, s.plot, rankogram))
names(x5) <- c("ranktable", "sucratable", "order", "longtable", "sucraplot", "rankogram")
return(x5)
}


