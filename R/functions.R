extract.names <- function(string){
  sapply(string, function(x) strsplit(x, "\\.")[[1]][2]) %>%
    as.character
}

calc.report <- function(x, fct="identity", arg=NULL, trans="identity"){
  if(is.null(arg)){
    eval(call(fct,(call(trans, x))))
    } else 
      eval(call(fct,(call(trans, x)), arg))
}

exp.mean <- function(x) calc.report(x, "mean", trans="exp")
exp.median <- function(x) calc.report(x, "median", trans="exp")
exp.sd <- function(x) calc.report(x, "sd", trans="exp")
exp.lci <- function(x) calc.report(x, "quantile", arg=0.025, trans="exp")
exp.uci <- function(x) calc.report(x, "quantile", arg=0.975, trans="exp")

id.mean <- function(x) calc.report(x, "mean", trans="identity")
id.median <- function(x) calc.report(x, "median", trans="identity")
id.sd <- function(x) calc.report(x, "sd", trans="identity")
id.lci <- function(x) calc.report(x, "quantile", arg=0.025, trans="identity")
id.uci <- function(x) calc.report(x, "quantile", arg=0.975, trans="identity")
