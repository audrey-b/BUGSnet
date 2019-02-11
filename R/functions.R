extract.names <- function(string){
  sapply(string, function(x) strsplit(x, "\\.")[[1]][2]) %>%
    as.character
}

calc.report <- function(x, fct="identity", arg=NULL, trans="identity", digits=16){
  if(is.null(arg)){
    round(eval(call(fct,(call(trans, x)))), digits)
    } else 
      round(eval(call(fct,(call(trans, x)), arg)), digits)
}

e.mean <- function(x) calc.report(x, "mean", trans="exp")
e.median <- function(x) calc.report(x, "median", trans="exp")
e.sd <- function(x) calc.report(x, "sd", trans="exp")
e.lci <- function(x) calc.report(x, "quantile", arg=0.025, trans="exp")
e.uci <- function(x) calc.report(x, "quantile", arg=0.975, trans="exp")

e.mean.round <- function(x) calc.report(x, "mean", trans="exp", digits=2)
e.median.round <- function(x) calc.report(x, "median", trans="exp", digits=2)
e.sd.round <- function(x) calc.report(x, "sd", trans="exp", digits=2)
e.lci.round <- function(x) calc.report(x, "quantile", arg=0.025, trans="exp", digits=2)
e.uci.round <- function(x) calc.report(x, "quantile", arg=0.975, trans="exp", digits=2)
