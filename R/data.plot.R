#' Patient Characteristic Plot
#' @description Plots a particular patient characteristic by study or by treatment. Useful for assessing
#' differences in potential effect modifiers.
#' 
#' @param patients.data A dataset containing patient level data. \code{patients.data = my.slr$patient.data} is a
#' convenient way to input this parameter
#' @param treatment.var A string indicating the name of the treatment variable
#' @param trial.var A string indicating the name of the study variable
#' @param var.name A string indicating the name of the patient characteristic to be plotted
#' @param fill.str An optional string indicating the variable to categorize measurements. For instance,
#' some studies report the mean treatment effect and others may report the median treatment effect. If there 
#' is a variable in \code{patients.data} called "type.measure" indicating whether the mean or median is 
#' reported, setting fill.str="type.measure" would colour all studies reporting the mean as red, and all
#' the studies reporting the median as turquoise.
#' @param errorbar If TRUE, will add error bars to your plot.
#' @param errorbar.min Mandatory if errorbar=TRUE. A string indicating how to calculate the minimum of the 
#' error bar (see Example 2).
#' @param errorbar.max Mandatory if errorbar=TRUE. A string indicating how to calculate the maximum of the 
#' error bar (see Example 2).
#' @param overall.avg If true, adds overall average line to plot.
#' @param y.lab y-label of the plot.
#' @param caption Caption of plot.
#' @param by If by="trial" then data from arms will be grouped by study/trial. If by="treatment" then
#' bar graph is grouped by treatment.
#' 
#' @examples
#' # Example containing a fill.str, an overall average, and no error
#' 
#'data.plot(patients.data = my.slr$patient.data, 
#'              treatment.var = "trt",
#'              trial.var = "studyName",
#'              var.name = "Age", 
#'              fill.str="type_measure",
#'              y.lab = "Age (Years)", 
#'              caption = "Johnson study did not report age",
#'              by = "trial")
#'              
#'# Example containing no fill.str, no overall average, but contains errorbars
#'
#'data.plot(patients.data = my.slr$patient.data, 
#'              treatment.var = "trt",
#'              trial.var = "studyName",
#'              var.name = "Age", 
#'              errorbar = TRUE,
#'              errorbar.min = "Age - Age_sd", 
#'              errorbar.max = "Age + Age_sd",
#'              overall.avg=FALSE,
#'              y.lab = "Age (Years)", 
#'              caption = "Error bars: Mean +/- sd",
#'              by = "trial")

data.plot <- function(patients.data, 
                          treatment.var,
                          trial.var,
                          var.name, 
                          errorbar = FALSE,
                          errorbar.min = NULL, 
                          errorbar.max = NULL,
                          overall.avg=TRUE,
                          fill.str=NULL,
                          y.lab, 
                          caption,
                          by = "trial"){
  
  if (errorbar == TRUE && (is.null(errorbar.min) | is.null(errorbar.max))){
    stop('When errorbar = TRUE, errorbar.min and errorbar.max must be specified')
    }
  
  if (errorbar){eb = geom_errorbar(aes_string(ymin=errorbar.min, ymax=errorbar.max))}
  else {eb = NULL}
  
  if (overall.avg){
    overall.mean <- patients.data %>% select(var.name)%>%colMeans( na.rm=TRUE) %>% as.numeric()
    om <- geom_hline(yintercept = overall.mean, color = "red", linetype=2)}
  else{om=NULL}
  
  patients.data <- patients.data %>%
    mutate(trt =(!! as.name(treatment.var)), trial =(!! as.name(trial.var)))
  
  if (by == "trial"){
  if(is.null(fill.str)) {p <- ggplot(patients.data, aes_string(x=treatment.var, y=var.name))}
  else p <- ggplot(patients.data, aes_string(x=treatment.var, y=var.name, fill=fill.str))
  
  p <- p + geom_bar(stat = "identity") +
    eb +
    om + 
    theme(axis.text.x = element_text(angle = 90, hjust = 1),
          legend.title=element_blank(),
          plot.margin=unit(c(25, 5.5, 5.5, 50), "points"),
          strip.text.x = element_text(angle = 90)) +
    labs(y = y.lab, x="", caption = caption)+
    facet_grid(. ~ trial,  space="free_x", scales="free_x")
  
  }
  else if (by == "treatment"){
    
  if(is.null(fill.str)) {p <- ggplot(patients.data, aes_string(x=trial.var, y=var.name))}
  else p <- ggplot(patients.data, aes_string(x=trial.var, y=var.name, fill=fill.str))
  
  p <- p + geom_bar(stat = "identity") +
    eb +
    om + 
    theme(axis.text.x = element_text(angle = 90, hjust = 1),
          legend.title=element_blank(),
         plot.margin=unit(c(10, 4, 4, 20), "points"),
          strip.text.x = element_text(angle = 90)) +
    labs(y = y.lab, x="", caption = caption)+
    facet_grid(. ~ trt,  space="free_x", scales="free_x")
  
  }
  p
}