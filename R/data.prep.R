#' Data preparation
#' @description Puts data into appropriate form for subsequent analysis using BUGSnet.
#' @param arm.data Data with 1 row for each study arm.
#' @param varname.t A string indicating the name of the treatment variable.
#' @param varname.s A string indicating the name of the study variable.
#' 
#' @return \code{data.prep} returns an object of class \code{BUGSnetData} which is a list containing the following components:
#' @return \code{arm.data} - A tibble containing the arm level study data
#' @return \code{treatments} - A list of all treatments in the network
#' @return \code{studies} - A list of all studies in the network
#' @return \code{n.arms} - A tibble containing the number of arms for each study
#' @return \code{varname.t} - A string containing the name of the treatment variable
#' @return \code{varname.s} - A string containing the name of the study variable
#' 
#' @examples 
#' data(diabetes.sim)
#' 
#' diabetes.slr <- data.prep(arm.data = diabetes.sim, 
#' varname.t = "Treatment", 
#' varname.s = "Study")
#' 
#' diabetes.slr$arm.data
#' diabetes.slr$treatments
#' diabetes.slr$studies
#' @export
#' @seealso \code{\link{net.plot}}, \code{\link{net.tab}}, \code{\link{data.plot}}, \code{\link{nma.model}}


data.prep <- function(arm.data,
                      #patient.data = NULL,
                      varname.t,
                      varname.s
){
  
  varname.t.quo <- rlang::quo(!! as.name(varname.t))
  varname.s.quo <- rlang::quo(!! as.name(varname.s))
  
  
  arm.data %<>% mutate(!! varname.t.quo := as.character(!! varname.t.quo),
                       !! varname.s.quo := as.character(!! varname.s.quo))
  
  bdata <- structure(list(arm.data=arm.data,
                          #patient.data = patient.data,
                          treatments=arm.data %>% select(varname.t) %>% unique %>% arrange(!! varname.t.quo),
                          studies=as.data.frame(arm.data %>% select(varname.s) %>% unique %>% arrange(!! varname.s.quo)),
                          n.arms= arm.data %>% group_by(!! varname.s.quo) %>% summarize(n.arms = n()) %>% ungroup() %>% arrange,
                          varname.t = varname.t,
                          varname.s = varname.s),
                     class = "BUGSnetData")
  
  return(bdata)
}
