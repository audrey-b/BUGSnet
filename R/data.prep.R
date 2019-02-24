
#' Data preparation
#' @description Puts data into appropriate form for subsequent analysis using BUGSnet.
#' @param arm.data Data with 1 row for each study arm.
#' @param varname.t A string indicating the name of the treatment variable.
#' @param varname.s A string indicating the name of the study variable.
#' @param N A string indicating the name of the variable containing the number of participants in each arm
#' 
#' @return \code{arm.data} - A tibble containing the arm level study data
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

data.prep <- function(arm.data,
                      #patient.data = NULL,
                      varname.t,
                      varname.s,
                      N
){
  return(list(arm.data=arm.data %>% mutate_if(is.factor, as.character),
              #patient.data = patient.data,
              treatments=arm.data %>% select_(varname.t) %>% unique %>% arrange_(varname.t),
              studies=arm.data %>% select_(varname.s) %>% unique %>% arrange_(varname.s),
              n.arms= arm.data %>% group_by_(varname.s) %>% summarize(n.arms = n()) %>% ungroup() %>% arrange,
              varname.t = varname.t,
              varname.s = varname.s)
  )
}
