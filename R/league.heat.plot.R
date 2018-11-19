#' League Heat Plot
#' @description Creates and plots a heatmap for a league table produced by \code{leaguetable()}. The user inputs a colour indicating
#' a negative relative treatment effect, a null relative treatment effect, and a large relative treatment effect.
#' @param leaguetable Results from \code{leaguetable()}. Be sure to set layout="long" from \code{leaguetable()} (see example). 
#' @param low.colour A string indicating the colour of negative relative treatment effects (e.g relative risk of 0.5).
#' @param mid.colour A string indicating the colour of null relative treatment effects (e.g relative risk of ~1.0). 
#' @param high.colour A string indicating the colour of high relative treatment effects (e.g relative risk of ~2.0).
#' @param midpoint Point indicating a null treatment effect. THIS NEEDS TO BE UPDATED. CURRENTLY ONLY WORKS FOR BINOMIAL VARIABLES.
#' 
#' @examples
#' #Make league table, be sure to set layout="long".
#' lt <- leaguetable(jagsoutput=nma.results, central.tdcy="median", layout = "long")
#' 
#' #make plot using default colours
#' league.heat.plot(leaguetable=lt)

league.heat.plot <- function(leaguetable,
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
                           ifelse(Treatment!=Comparator, paste0("**", median," ", "(",lci, ", ", uci,")","**"), " "),
                           ifelse(Treatment!=Comparator, paste0(median," ", "(",lci, ", ", uci,")"), " "))),
              size=3)+
    scale_fill_gradient2(low = low.colour, high = high.colour, mid = mid.colour, midpoint = midpoint)+
    theme_dark()+
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
          legend.position="none", panel.border=element_blank(),
          axis.ticks.x=element_blank(),
          axis.ticks.y=element_blank())+
    scale_x_discrete(expand = c(0, 0)) +
    scale_y_discrete(expand = c(0, 0))
}

