#' @noRd
extract.names <- function(string){
  sapply(string, function(x) strsplit(x, "\\.")[[1]][2]) %>%
    as.character
}

#' @noRd
calc.report <- function(x, fct="identity", arg=NULL, trans="identity"){
  if(is.null(arg)){
    eval(call(fct,(call(trans, x))))
    } else 
      eval(call(fct,(call(trans, x)), arg))
}


#' @noRd
exp.mean <- function(x) calc.report(x, "mean", trans="exp")

#' @noRd
exp.median <- function(x) calc.report(x, "median", trans="exp")

#' @noRd
exp.sd <- function(x) calc.report(x, "sd", trans="exp")

#' @noRd
exp.lci <- function(x) calc.report(x, "quantile", arg=0.025, trans="exp")

#' @noRd
exp.uci <- function(x) calc.report(x, "quantile", arg=0.975, trans="exp")


#' @noRd
id.mean <- function(x) calc.report(x, "mean", trans="identity")

#' @noRd
id.median <- function(x) calc.report(x, "median", trans="identity")

#' @noRd
id.sd <- function(x) calc.report(x, "sd", trans="identity")

#' @noRd
id.lci <- function(x) calc.report(x, "quantile", arg=0.025, trans="identity")

#' @noRd
id.uci <- function(x) calc.report(x, "quantile", arg=0.975, trans="identity")
