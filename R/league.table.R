#' League Table
#' @description Produces league table that contains all information about relative effectiveness and their 
#' uncertainty for all possible pairs of treatments.
#' @param jagsoutput Results from running \code{nma.analysis()}.
#' @param central.tdcy The statistic that you want to use in order to measure relative effectiveness. The options are "mean" and "median".
#' @param layout Indicates the form of the results. The options are "default" and "long". Use "default" if you want an nxn matrix
#' comparing all n treatments. Use "long" if you want 1 line per treatment comparison. Note that if you want to use the function 
#' \code{league.heat.plot}, set this parameter to "long".
#' @seealso \code{league.heat.plot()}
#' 
#' @examples
#'
#' league_table <- leaguetable(jagsoutput=nma.results, central.tdcy="median", layout = "long")

leaguetable <- function(jagsoutput, 
                        central.tdcy = "median",
                        layout="default") {
  
x <- do.call(rbind, jagsoutput$samples) %>% data.frame() %>% select(starts_with("d."))
trt.names <- jagsoutput$trt.key
colnames(x) <- trt.names
  
colvals <- function(x, b.col=1, paste=TRUE) {
  
  base <- colnames(x)[b.col]

  x2 <- x
  new.vars <- paste0(colnames(x2), "-", b.col)
 for(i in 1:ncol(x)) {
  x2[[new.vars[i]]] <- x[, i] - x[, b.col]
          }
  x2 %<>% select(new.vars)
  colnames(x2) <- trt.names

 if(central.tdcy=="mean"){
 tmp.estimate <- x2 %>%  
    summarise_all(funs(estimate = e.mean.round)) %>% gather() %>%
    rename(trt = key, estimate = value) %>%
    mutate(trt = sub("_estimate", "", trt))
 }else if(central.tdcy=="median"){
 tmp.estimate <- x2 %>%  
   summarise_all(funs(estimate = e.median.round)) %>% gather() %>%
   rename(trt = key, estimate = value) %>%
   mutate(trt = sub("_estimate", "", trt))
 }
 
 tmp.lci <- x2 %>%  
   summarise_all(funs(lci = e.lci.round)) %>% gather() %>%
   rename(trt = key, lci = value) %>%
   mutate(trt = sub("_lci", "", trt))
 
 tmp.uci <- x2 %>%  
   summarise_all(funs(uci = e.uci.round)) %>% gather() %>%
   rename(trt = key, uci = value) %>%
   mutate(trt = sub("_uci", "", trt))
 
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
 
colnames(tmp1)[2] <- as.character(tmp.estimate %>% filter(estimate==1) %>% select(trt))
 } else{
   tmp1 <- left_join(tmp.estimate, tmp.lci, by = "trt") %>%
     left_join(., tmp.uci, by = "trt")
   
   colnames(tmp1)[2] <- central.tdcy
 }

return(tmp1)
}


if(layout=="default"){
  
  tmp1.list <- list()
  for(i in 1:ncol(x)) {
    tmp1.list[[i]] <- colvals(x, b.col=i)
  }
 
  tmp.df <- bind_cols(tmp1.list) %>%
    select(-starts_with("trt")) %>%
    t()
  
  colnames(tmp.df) <- colnames(x)
  rownames(tmp.df) <- colnames(x)

for(i in 1:dim(tmp.df)[1]){
  tmp.df[i,i] <- colnames(tmp.df)[i]
}

return(tmp.df) 
} else if(layout=="long"){
  
  tmp1.list <- list()
  for(i in 1:ncol(x)) {
    tmp1.list[[i]] <- colvals(x, b.col=i, paste=FALSE)
  }
  
  tble <- tmp1.list %>% 
    bind_rows() %>%
    mutate(Treatment = trt,
           Comparator = rep(trt.names, each=length(trt.names))) %>%
    select(Treatment, Comparator, everything(), -trt)

  return(tble)
}


}





 
 
