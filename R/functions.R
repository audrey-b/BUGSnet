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

exp.mean <- function(x) calc.report(x, "mean", trans="exp")
exp.median <- function(x) calc.report(x, "median", trans="exp")
exp.sd <- function(x) calc.report(x, "sd", trans="exp")
exp.lci <- function(x) calc.report(x, "quantile", arg=0.025, trans="exp")
exp.uci <- function(x) calc.report(x, "quantile", arg=0.975, trans="exp")

exp.mean.round <- function(x) calc.report(x, "mean", trans="exp", digits=2)
exp.median.round <- function(x) calc.report(x, "median", trans="exp", digits=2)
exp.sd.round <- function(x) calc.report(x, "sd", trans="exp", digits=2)
exp.lci.round <- function(x) calc.report(x, "quantile", arg=0.025, trans="exp", digits=2)
exp.uci.round <- function(x) calc.report(x, "quantile", arg=0.975, trans="exp", digits=2)

mean.round <- function(x) calc.report(x, "mean", trans="identity", digits=2)
median.round <- function(x) calc.report(x, "median", trans="identity", digits=2)
sd.round <- function(x) calc.report(x, "sd", trans="identity", digits=2)
lci.round <- function(x) calc.report(x, "quantile", arg=0.025, trans="identity", digits=2)
uci.round <- function(x) calc.report(x, "quantile", arg=0.975, trans="identity", digits=2)