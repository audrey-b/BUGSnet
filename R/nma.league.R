#' League Table and Heat Plot
#' @description Produces a league table and a league heat plot that contain point estimates of relative effectiveness 
#' for all possible pairs of treatments point estimates along with 95% credible intervals obtained with the quantile method.
#' @param nma A \code{BUGSnetRun} object produced by running \code{nma.run()}.
#' @param central.tdcy The statistic that you want to use in order to measure relative effectiveness. The options are "mean" and "median".
#' @param log.scale If TRUE, odds ratios, relative risk or hazard ratios are reported on the log scale. Default is FALSE.
#' @param order A vector of strings representing the order in which to display the treatments.
#' @param low.colour A string indicating the colour of low relative treatment effects for the heat plot (e.g relative risk of 0.5).
#' @param mid.colour A string indicating the colour of null relative treatment effects for the heat plot (e.g relative risk of ~1.0). 
#' @param high.colour A string indicating the colour of high relative treatment effects for the heat plot (e.g relative risk of ~2.0).
#' @param cov.value  Must be specified for meta-regression. This is the value of the covariate for which to report the results.
#' @param digits Number of digits to display after the decimal point
#' 
#' @return \code{table} - A league table. Row names indicate comparator treatments.
#' @return \code{longtable} - League table in the long format.
#' @return \code{heatplot} - League heat plot, where a color scale is used to represent relative treatment effects and ** are used to highlight statistically significant differences.
#' 
#' @examples
#' data(diabetes.sim)
#' 
#' diabetes.slr <- data.prep(arm.data = diabetes.sim, 
#' varname.t = "Treatment", 
#' varname.s = "Study")
#' 
#' #Random effects, consistency model.
#' #Binomial family, cloglog link. This implies that the scale will be the Hazard Ratio.
#'diabetes.re.c <- nma.model(data = diabetes.slr,
#'        outcome = "diabetes", 
#'        N = "n",
#'        reference = "Placebo",
#'        family = "binomial",
#'        link = "cloglog",
#'        effects = "random",
#'        type="consistency",
#'        time="followup"
#'        )
#'  
#'diabetes.re.c.res <- nma.run(diabetes.re.c,
#'n.adapt=100,
#'n.burnin=0,
#'n.iter=100)
#'
#'  
#' league_table <- nma.league(nma=diabetes.re.c.res, central.tdcy="median")
#' league_table$heatplot
#' league_table$table
#' league_table$longtable
#' @export
#' @seealso \code{\link{nma.run}}, \code{\link{nma.rank}}, \code{\link{nma.forest}} 


nma.league <- function(nma, 
                       central.tdcy = "median",
                       log.scale = FALSE,
                       order = NULL,
                       low.colour = "darkgoldenrod1", 
                       mid.colour = "white",
                       high.colour = "cornflowerblue",
                       cov.value=NULL,
                       digits = 2) {
  
  # Bind variables to function
  trt <- NULL
  Treatment <- NULL
  Comparator <- NULL
  
  
  if (class(nma) != "BUGSnetRun")
    stop("\'nma\' must be a valid BUGSnetRun object created using the nma.run function.")
  
  if(!is.null(nma$model$covariate) & is.null(cov.value)) stop("cov.value must be specified for meta-regression")
  
  if(nma$model$type=="inconsistency") {
    stop('This function has not been implemented yet for the inconsistency model.')
  } 
  
  dmat <- do.call(rbind, nma$samples) %>% data.frame() %>% select(starts_with("d."))
  trt.names <- nma$trt.key
  colnames(dmat) <- trt.names
  
  if(!is.null(cov.value)){#meta-regression
    
    betamat <- do.call(rbind, nma$samples) %>% data.frame() %>% select(starts_with("beta."))
    trt.names <- nma$trt.key
    colnames(betamat) <- trt.names
    
    dmat <- dmat+betamat*(cov.value-nma$model$mean.cov)
  }
  
  if(is.null(order)){
    order <- trt.names
  }
  
  dmat %<>% select(order)
  trt.names <- order
  
  colvals <- function(dmat, b.col=1, paste=TRUE) {
    
    # Bind variables to function
    key <- NULL
    value <- NULL
    trt <- NULL
    estimate <- NULL
    lci <- NULL
    uci <- NULL
    result <- NULL
    
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
        summarise_all(list(lci = exp.lci)) %>% gather() %>%
        rename(trt = key, lci = value) %>%
        mutate(trt = sub("_lci", "", trt))
      
      tmp.uci <- dmat2 %>%  
        summarise_all(list(uci = exp.uci)) %>% gather() %>%
        rename(trt = key, uci = value) %>%
        mutate(trt = sub("_uci", "", trt))
      
      null.value <- 1
    } else{
      tmp.lci <- dmat2 %>%  
        summarise_all(list(lci = id.lci)) %>% gather() %>%
        rename(trt = key, lci = value) %>%
        mutate(trt = sub("_lci", "", trt))
      
      tmp.uci <- dmat2 %>%  
        summarise_all(list(uci = id.uci)) %>% gather() %>%
        rename(trt = key, uci = value) %>%
        mutate(trt = sub("_uci", "", trt))
      
      null.value <- 0
    }
    
    #create C-style formatting string from the digits parameter
    fmt <- paste0("%.", digits, "f")
    
    if(paste){
      tmp1 <- left_join(tmp.estimate, tmp.lci, by = "trt") %>%
        left_join(., tmp.uci, by = "trt") %>%
        mutate(result = paste(sprintf(fmt, estimate), 
                              " (", 
                              sprintf(fmt, lci),
                              " to ", 
                              sprintf(fmt, uci), 
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
  
  tmp.df <- suppressMessages(bind_cols(tmp1.list)) %>%
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
                                                                                  midpoint = null.value,
                                                                                  digits = digits)))
}





league.heat.plot <- function(leaguetable,
                             central.tdcy,
                             order = NULL,
                             low.colour = "red", 
                             mid.colour = "white",
                             high.colour = "springgreen4",
                             midpoint,
                             digits){
  
  #Bind Variables to function
  Treatment <- NULL
  Comparator <- NULL
  ct.stat <- NULL
  lci <- NULL
  uci <- NULL
  
  #if (ncol(leaguetable) > 5){warning("leaguetable must be in 'long' format")}
  
  league.tmp <- leaguetable%>%filter(Treatment != Comparator)
  
  #rename central tendancy statistic to ct.stat
  names(league.tmp)[names(league.tmp) == central.tdcy] <- "ct.stat"
  #create C-style formatting string from the digits parameter
  fmt <- paste0("%.", digits, "f")
  
  heatplot <- ggplot(data = league.tmp, aes(x=Treatment, y=Comparator, fill=ct.stat)) + 
    geom_tile()+
    #coord_flip()+
    geom_text(aes(label = 
                    ifelse(((midpoint<lci & midpoint< uci) | (midpoint>lci & midpoint> uci)),
                           ifelse(Treatment!=Comparator, paste0("**", sprintf(fmt, ct.stat), "**", "\n", "(",sprintf(fmt, lci), ", ", sprintf(fmt, uci),")"), " "),
                           ifelse(Treatment!=Comparator, paste0(sprintf(fmt, ct.stat), "\n", "(",sprintf(fmt, lci), ", ", sprintf(fmt, uci),")"), " "))),
              size=3)
  
  heatplot <- heatplot +
    scale_fill_gradient2(low = low.colour, high = high.colour, mid = mid.colour, midpoint = midpoint)+
    theme_dark()+
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
          legend.position="none", panel.border=element_blank(),
          axis.ticks.x=element_blank(),
          axis.ticks.y=element_blank())+
    scale_x_discrete(limits = order, expand = c(0, 0), position="top") +
    scale_y_discrete(limits = rev(order), expand = c(0, 0)) 
  
  return(heatplot)
}

