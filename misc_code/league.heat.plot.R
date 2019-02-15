

league.heat.plot <- function(leaguetable,
                             sucra.ranks = NULL,
                             low.colour = "red", 
                             mid.colour = "white",
                             high.colour = "springgreen4",
                             midpoint = 0){
  
if (ncol(leaguetable) > 5){warning("leaguetable must be in 'long' format")}
  
  league.tmp <- leaguetable%>%filter(Treatment != Comparator)
  
  ggplot(data = league.tmp, aes(x=Treatment, y=Comparator, fill=log(median))) + 
    geom_tile()+
    coord_flip()+
    geom_text(aes(label = 
                    ifelse(((1<lci & 1< uci) | (1>lci & 1> uci)),
                           ifelse(Treatment!=Comparator, paste0("**", median, "**", "\n", "(",lci, ", ", uci,")"), " "),
                           ifelse(Treatment!=Comparator, paste0(median, "\n", "(",lci, ", ", uci,")"), " "))),
              size=3)+
    scale_fill_gradient2(low = low.colour, high = high.colour, mid = mid.colour, midpoint = midpoint)+
    theme_dark()+
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
          legend.position="none", panel.border=element_blank(),
          axis.ticks.x=element_blank(),
          axis.ticks.y=element_blank())+
    scale_x_discrete(limits = sucra.ranks, expand = c(0, 0)) +
    scale_y_discrete(limits = sucra.ranks, expand = c(0, 0))
}

