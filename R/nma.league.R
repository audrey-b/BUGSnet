#' League Table and Heat Plot
#' @description Produces a league table and a league heat plot that contain point estimates of relative effectiveness 
#' for all possible pairs of treatments point estimates along with 95% credible intervals obtained with the quantile method.
#' @param nma Results from running \code{nma.analysis()}.
#' @param central.tdcy The statistic that you want to use in order to measure relative effectiveness. The options are "mean" and "median".
#' @param log.scale If TRUE, odds ratios, relative risk or hazard ratios are reported on the log scale. Default is FALSE.
#' @param order A vector of strings representing the order in which to display the treatments.
#' @param low.colour A string indicating the colour of low relative treatment effects for the heat plot (e.g relative risk of 0.5).
#' @param mid.colour A string indicating the colour of null relative treatment effects for the heat plot (e.g relative risk of ~1.0). 
#' @param high.colour A string indicating the colour of high relative treatment effects for the heat plot (e.g relative risk of ~2.0).
#' @param x  Must be specified for meta-regression. This is the value of the covariate for which to report the results.
#' 
#' @return \code{table} - League table.
#' @return \code{longtable} - League table in the long format.
#' @return \code{heatplot} - League heat plot, where a color scale is used to represent relative treatment effects and ** are used to highlight statistically significant differences.
#' 
#' @examples
#'
#' league_table <- leaguetable(nma=nma.results, central.tdcy="median")

nma.league <- function(nma, 
                       central.tdcy = "median",
                       log.scale = FALSE,
                       order = NULL,
                       low.colour = "darkgoldenrod1", 
                       mid.colour = "white",
                       high.colour = "cornflowerblue",
                       x=NULL) {
  
if(!is.null(nma$model$covariate) & is.null(x)) stop("x must be specified for meta-regression")
  
  dmat <- do.call(rbind, nma$samples) %>% data.frame() %>% select(starts_with("d."))
trt.names <- nma$trt.key
colnames(dmat) <- trt.names

if(!is.null(x)){#meta-regression
  
  betamat <- do.call(rbind, nma$samples) %>% data.frame() %>% select(starts_with("beta."))
  trt.names <- nma$trt.key
  colnames(betamat) <- trt.names
  
  dmat <- dmat+betamat*(x-nma$model$mean.x)
}

if(is.null(order)){
  order <- trt.names
}

dmat %<>% select(order)
trt.names <- order

colvals <- function(dmat, b.col=1, paste=TRUE) {
  
  base <- colnames(dmat)[b.col]
  
  dmat2 <- dmat
  new.vars <- paste0(colnames(dmat2), "-", b.col)
  for(i in 1:ncol(dmat)) {
    dmat2[[new.vars[i]]] <- dmat[, i] - dmat[, b.col]
  }
  dmat2 %<>% select(new.vars)
  colnames(dmat2) <- trt.names
  
  if(central.tdcy=="mean" & log.scale==FALSE & nma$link!="identity"){
    tmp.estimate <- dmat2 %>%  
      summarise_all(list(estimate = exp.mean)) %>% gather() %>%
      rename(trt = key, estimate = value) %>%
      mutate(trt = sub("_estimate", "", trt))
  } else if(central.tdcy=="mean"){
    tmp.estimate <- dmat2 %>%  
      summarise_all(list(estimate = id.mean)) %>% gather() %>%
      rename(trt = key, estimate = value) %>%
      mutate(trt = sub("_estimate", "", trt))}
  if(central.tdcy=="median" & log.scale==FALSE  & nma$link!="identity"){
    tmp.estimate <- dmat2 %>%  
      summarise_all(list(estimate = exp.median)) %>% gather() %>%
      rename(trt = key, estimate = value) %>%
      mutate(trt = sub("_estimate", "", trt))
  } else if(central.tdcy=="median"){
    tmp.estimate <- dmat2 %>%  
      summarise_all(list(estimate = id.median)) %>% gather() %>%
      rename(trt = key, estimate = value) %>%
      mutate(trt = sub("_estimate", "", trt))
  }
  
  if(log.scale==FALSE & nma$link!="identity"){
  tmp.lci <- dmat2 %>%  
    summarise_all(funs(lci = exp.lci)) %>% gather() %>%
    rename(trt = key, lci = value) %>%
    mutate(trt = sub("_lci", "", trt))
  
  tmp.uci <- dmat2 %>%  
    summarise_all(funs(uci = exp.uci)) %>% gather() %>%
    rename(trt = key, uci = value) %>%
    mutate(trt = sub("_uci", "", trt))
  
  null.value <- 1
  } else{
    tmp.lci <- dmat2 %>%  
      summarise_all(funs(lci = id.lci)) %>% gather() %>%
      rename(trt = key, lci = value) %>%
      mutate(trt = sub("_lci", "", trt))
    
    tmp.uci <- dmat2 %>%  
      summarise_all(funs(uci = id.uci)) %>% gather() %>%
      rename(trt = key, uci = value) %>%
      mutate(trt = sub("_uci", "", trt))
    
    null.value <- 0
  }
  
  if(paste){
    tmp1 <- left_join(tmp.estimate, tmp.lci, by = "trt") %>%
      left_join(., tmp.uci, by = "trt") %>%
      mutate(result = paste(format(estimate,drop0Trailing = F), 
                            " (", 
                            format(lci,drop0Trialing = F),
                            " to ", 
                            format(uci,drop0Triailing = F), 
                            ")", sep="")) %>%
      select(trt, result)
    
    colnames(tmp1)[2] <- as.character(tmp.estimate %>% filter(estimate==null.value) %>% select(trt))
  } else{
    tmp1 <- left_join(tmp.estimate, tmp.lci, by = "trt") %>%
      left_join(., tmp.uci, by = "trt")
    
    colnames(tmp1)[2] <- central.tdcy
  }
  
  return(tmp1)
}


#Default layout

tmp1.list <- list()
for(i in 1:ncol(dmat)) {
  tmp1.list[[i]] <- colvals(dmat, b.col=i)
}

tmp.df <- bind_cols(tmp1.list) %>%
  select(-starts_with("trt")) %>%
  t()

colnames(tmp.df) <- colnames(dmat)
rownames(tmp.df) <- colnames(dmat)

for(i in 1:dim(tmp.df)[1]){
  tmp.df[i,i] <- colnames(tmp.df)[i]
}

default <- tmp.df

#long layout

tmp2.list <- list()
for(i in 1:ncol(dmat)) {
  tmp2.list[[i]] <- colvals(dmat, b.col=i, paste=FALSE)
}

longtable <- tmp2.list %>% 
  bind_rows() %>%
  mutate(Treatment = trt,
         Comparator = rep(trt.names, each=length(trt.names))) %>%
  select(Treatment, Comparator, everything(), -trt)

if(log.scale==FALSE & nma$link!="identity"){
  null.value <- 1
} else{
  null.value <- 0
}



return(list("table"=default, "longtable"=longtable, "heatplot"=league.heat.plot(leaguetable=longtable,
                                                                                  central.tdcy=central.tdcy,
                                                                            order = order,
                                                                            low.colour = low.colour, 
                                                                            mid.colour = mid.colour,
                                                                            high.colour = high.colour,
                                                                            midpoint = null.value)))
}





league.heat.plot <- function(leaguetable,
                             central.tdcy,
                             order = NULL,
                             low.colour = "red", 
                             mid.colour = "white",
                             high.colour = "springgreen4",
                             midpoint){
  
  if (ncol(leaguetable) > 5){warning("leaguetable must be in 'long' format")}
  
  league.tmp <- leaguetable%>%filter(Treatment != Comparator)
  
  if(central.tdcy=="mean"){
    heatplot <- ggplot(data = league.tmp, aes(dmat=Treatment, y=Comparator, fill=mean)) + 
      geom_tile()+
      coord_flip()+
      geom_text(aes(label = 
                      ifelse(((midpoint<lci & midpoint< uci) | (midpoint>lci & midpoint> uci)),
                             ifelse(Treatment!=Comparator, paste0("**", formatC(mean, 2), "**", "\n", "(",formatC(lci, 2), ", ", formatC(uci, 2),")"), " "),
                             ifelse(Treatment!=Comparator, paste0(formatC(mean, 2), "\n", "(",formatC(lci, 2), ", ", formatC(uci, 2),")"), " "))),
                size=3)
  } else if(central.tdcy=="median"){
    heatplot <- ggplot(data = league.tmp, aes(x=Treatment, y=Comparator, fill=median)) + 
      geom_tile()+
      coord_flip()+
      geom_text(aes(label = 
                      ifelse(((midpoint<lci & midpoint< uci) | (midpoint>lci & midpoint> uci)),
                             ifelse(Treatment!=Comparator, paste0("**", formatC(median, 2), "**", "\n", "(",formatC(lci, 2), ", ", formatC(uci, 2),")"), " "),
                             ifelse(Treatment!=Comparator, paste0(formatC(median,2), "\n", "(",formatC(lci, 2), ", ", formatC(uci, 2),")"), " "))),
                size=3)
  }
  
  heatplot <- heatplot +
    scale_fill_gradient2(low = low.colour, high = high.colour, mid = mid.colour, midpoint = midpoint)+
    theme_dark()+
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
          legend.position="none", panel.border=element_blank(),
          axis.ticks.x=element_blank(),
          axis.ticks.y=element_blank())+
    scale_x_discrete(limits = order, expand = c(0, 0)) +
    scale_y_discrete(limits = order, expand = c(0, 0)) 
  
  return(heatplot)
}

