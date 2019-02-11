
#' Data preparation
#' @description Puts data into appropriate form for subsequent analysis using BUGSnet.
#' @param arm.data Data with 1 row for each study arm.
#' @param patient.data An optional separate data object containing patient characteristics (e.g comorbidities).
#' Data should be in same format as \code{arm.data} (1 row for each arm). Patient characteristic data should 
#' be plotted and compared when assessing the feasibility of an NMA.
#'
#' @param varname.t A string indicating the name of the treatment variable.
#' @param varname.s A string indicating the name of the study variable.
#' @param N A string indicating the name of the variable containing the number of participants in each arm
#' 
#' @return \code{arm.data} - A tibble containing the arm level study data
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
#' diabetes.slr <- data.prep(arm.data = my.data, varname.t = "treatment", varname.s = "studyName", N="n_patients")
#' 
#' # Example 2, contains patient data
#' NMA.data <- read.csv("MI_NMAdata.csv")
#' patient.data <- read.csv("patient_covariates.csv")
#' myocardInf.slr <- data.prep(arm.data = NMA.data, varname.t = "trt", varname.s = "stud", N="N")

data.prep <- function(arm.data,
                      patient.data = NULL,
                      varname.t,
                      varname.s,
                      N
){
  return(list(arm.data=arm.data,
              patient.data = patient.data,
              treatments=arm.data %>% select_(varname.t) %>% unique %>% arrange_(varname.t),
              studies=arm.data %>% select_(varname.s) %>% unique %>% arrange_(varname.s),
              n.arms= arm.data %>% group_by_(varname.s) %>% summarize(n.arms = n()) %>% ungroup() %>% arrange,
              varname.t = varname.t,
              varname.s = varname.s)
  )
}
