#' Table and Plots of Treatment Rankings
#' @description Produces a SUCRA (Surface Under the Cumulative Ranking Curve) plot and table. A Sucra table summarizes the probabilities
#' that each treatment is the best, second best...worst treatment in the network.
#' 
#' @param nma A \code{BUGSnetRun} object produced by running \code{nma.run()}.
#' @param largerbetter A boolean variable indicating whether a larger probability should indicate a more effective treatment (TRUE) or
#' if a smaller probability should indicate a more effective treatment (FALSE). 
#' @param sucra.lwd Line width relative to the default (default=1) in the SUCRA plot.
#' @param sucra.palette A string indicating the colour set from RcolorBrewer for the SUCRA plot. "Set1" is used by default and is automatically extended if 
#' there are many treatments in your network.
#' @param ranko.palette A string indicating the colour set from RcolorBrewer for the rankogram. "Blues" is used by default and is automatically extended if 
#' there are many treatments in your network.
#' @param cov.value  Must be specified for meta-regression. This is the value of the covariate for which to report the results.
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
#' #get sucra results
#' sucra_results <- nma.rank(nma = diabetes.re.c.res, largerbetter = FALSE)
#' 
#' #plot sucra results
#' sucra_results$sucraplot
#' @export
#' @seealso \code{\link{nma.run}}, \code{\link{nma.league}}, \code{\link{nma.forest}} 




nma.rank <- function(nma, 
                     largerbetter, 
                     sucra.lwd = 1.0,
                     sucra.palette= "Set1",
                     ranko.palette="Blues",
                     cov.value=NULL) {
  
  #Bind variables to function
  trt <- NULL
  value <- NULL
  iteration <- NULL
  prob <- NULL
  SUCRA <- NULL
  cumprob <- NULL
  
  if (class(nma) != "BUGSnetRun")
    stop("\'nma\' must be a valid BUGSnetRun object created using the nma.run function.")
  
  if(!is.null(nma$model$covariate) & is.null(cov.value)){
    if(nma$model$prior.beta!="EQUAL") stop("cov.value must be specified for meta-regression")
  }
  if(is.null(nma$model$covariate) & !is.null(cov.value)) stop("cov.value cannot be specified outside of meta-regression")
  
  if(nma$model$type=="inconsistency") {
    stop('This function has not been implemented yet for the inconsistency model.')
  } 
  
  x <- do.call(rbind, nma$samples) %>% data.frame() %>% select(starts_with("d."))
  colnames(x) <- nma$trt.key
  
  if(!is.null(cov.value)){#meta-regression
    
    betamat <- do.call(rbind, nma$samples) %>% data.frame() %>% select(starts_with("beta."))
    colnames(betamat) <- nma$trt.key
    
    x <- x+betamat*(cov.value-nma$model$mean.cov)
  }
  
  
  #tmp.var <- vector(mode="character", length=ncol(x) - 1)
  #for(i in 1:ncol(x)) {
  #  tmp.var[i] <- i
  #}
  
  tmp.var <- 1:ncol(x) %>% as.character()
  
  
  
  
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
  
  ##Extend sucra & ranko color palettes if necessary
  n.trts <- length(nma$trt.key)
  max.colors.sucra <- brewer.pal.info[sucra.palette,]$maxcolors
  if (max.colors.sucra < n.trts){
    tmp.colors <- brewer.pal(max.colors.sucra, sucra.palette)
    sucra.colors <- colorRampPalette(tmp.colors)(n.trts)
  } else{
    sucra.colors <- brewer.pal(n.trts, sucra.palette)
  }
  
  max.colors.ranko <- brewer.pal.info[ranko.palette,]$maxcolors
  if (max.colors.ranko < n.trts){
    tmp.colors <- brewer.pal(max.colors.ranko, ranko.palette)
    ranko.colors <- colorRampPalette(tmp.colors)(n.trts)
  } else{
    ranko.colors <- brewer.pal(n.trts, ranko.palette)
  }
  
  s.plot <- ggplot(data=x5, aes(x=factor(rank), y=cumprob, group=trt)) +
    geom_line(aes(color=trt), size=sucra.lwd) +
    geom_point(aes(color=trt)) +
    scale_color_manual(values = sucra.colors)+
    theme_bw()
  
  rankogram <- ggplot(data=x4, aes(y=prob, x=trt, fill=factor(rank)))+ 
    geom_bar(stat = "identity", width = 0.5) +
    scale_fill_manual(values=ranko.colors)
  
  if(!is.null(cov.value)){#meta-regression
    cov.str <- paste0(" when ",nma$model$covariate,"=",cov.value)
  }else cov.str <- ""
  
  if(largerbetter==TRUE){
    s.plot <- s.plot +
      labs(x=paste0("Ranking of Treatment",
                    "\n(Higher rankings associated with larger outcome values)"), 
           y=paste0("Probability of ranking or better (%)", cov.str),
           color="Treatment")
    rankogram <- rankogram + labs(x=paste0("Treatment",
                                           "\n(Higher rankings associated with larger outcome values)"), 
                                  y=paste0("Probability of ranking (%)", cov.str),
                                  fill="Ranking")
    
  }else if(largerbetter==FALSE){
    s.plot <- s.plot +
      labs(x=paste0("Ranking of Treatment",
                    "\n(Higher rankings associated with smaller outcome values)"), 
           y=paste0("Probability of ranking or better (%)", cov.str),
           color="Treatment")
    rankogram <- rankogram + labs(x=paste0("Treatment",
                                           "\n(Higher rankings associated with smaller outcome values)"), 
                                  y=paste0("Probability of ranking (%)", cov.str),
                                  fill="Ranking")
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


