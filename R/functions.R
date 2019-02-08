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

#' Data preparation
#' @description Puts data into appropriate form for subsequent analysis using BUGSnet.
#' @param raw.data Data with 1 row for each study arm.
#' @param patient.data An optional separate data object containing patient characteristics (e.g comorbidities).
#' Data should be in same format as \code{raw.data} (1 row for each arm). Patient characteristic data should 
#' be plotted and compared when assessing the feasibility of an NMA.
#'
#' @param varname.t A string indicating the name of the treatment variable.
#' @param varname.s A string indicating the name of the study variable.
#' @param N A string indicating the name of the variable containing the number of participants in each arm
#' 
#' @return \code{raw.data} - A tibble containing the arm level study data
#' @return \code{patient.data} - A tibble containing the arm level patient characteristics (if applicable)
#' @return \code{treatments} - A list of all treatments in the network
#' @return \code{studies} - A list of all studies in the network
#' @return \code{n.arms} - A tibble containing the number of arms for each study
#' @return \code{varname.t} - A string containing the name of the treatment variable
#' @return \code{varname.s} - A string containing the name of the study variable
#' 
#' @examples
#' # Example 1, no patient data
#' my.data <- read.csv("diabetes_NMAdata.csv")
#' diabetes.slr <- data.prep(raw.data = my.data, varname.t = "treatment", varname.s = "studyName", N="n_patients")
#' 
#' # Example 2, contains patient data
#' NMA.data <- read.csv("MI_NMAdata.csv")
#' patient.data <- read.csv("patient_covariates.csv")
#' myocardInf.slr <- data.prep(raw.data = NMA.data, varname.t = "trt", varname.s = "stud", N="N")

data.prep <- function(raw.data,
                     patient.data = NULL,
                     varname.t,
                     varname.s,
                     N
                     ){
  return(list(raw.data=raw.data,
              patient.data = patient.data,
              treatments=raw.data %>% select_(varname.t) %>% unique %>% arrange_(varname.t),
              studies=raw.data %>% select_(varname.s) %>% unique %>% arrange_(varname.s),
              n.arms= raw.data %>% group_by_(varname.s) %>% summarize(n.arms = n()) %>% ungroup() %>% arrange,
              varname.t = varname.t,
              varname.s = varname.s)
         )
}
