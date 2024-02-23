#' @title 
#' Patient Characteristic Plot
#' 
#' @description
#' Plots a particular patient characteristic by study or by treatment. Useful for assessing
#' differences in potential effect modifiers.
#' 
#' @param data A \code{BUGSnetData} object produced by \code{data.prep()}
#' @param covariate A string indicating the name of the patient characteristic to be plotted
#' @param covariate.label An optional string indicating the label of the patient characteristic to be displayed. 
#' Default is to display the covariate parameter as the axis label.
#' @param half.length A string indicating how to calculate the half-length of error bars (optional)
#' @param by If by="study" then data from arms will be grouped by study/trial. If by="treatment" then
#' bar graph is grouped by treatment.
#' @param fill.str An optional string indicating the variable to categorize measurements. For instance,
#' some studies report the mean treatment effect and others may report the median treatment effect. If there 
#' is a variable in \code{data} called "type.measure" indicating whether the mean or median is 
#' reported, setting fill.str="type.measure" would colour all studies reporting the mean as red, and all
#' the studies reporting the median as turquoise.
#' @param avg.hline If TRUE, adds overall average line to plot. Default is TRUE.
#' @param text.size Font size of the text. Default is 20.
#' @param orientation An optional string indicating the orientation of the plot. Default is "landscape".
#' (The plot expands width-ways with additional studies). 
#' If "portrait", the plot is rotated 90° and expands length-ways with additional studies.
#'              
#' @seealso
#' \code{\link{data.prep}}
#' 
#' @importFrom dplyr mutate select
#' @importFrom ggplot2 aes_string element_blank element_text facet_grid ggplot geom_errorbar geom_hline geom_point labs unit
#' @importFrom magrittr %>%
#' @importFrom rlang !!
#' 
#' @examples
#' data(diabetes.sim)
#' 
#' diabetes.slr <- data.prep(arm.data = diabetes.sim, 
#' varname.t = "Treatment", 
#' varname.s = "Study")
#'
#' # Example containing a fill.str, an overall average, and no error
#' 
#'data.plot(data = diabetes.slr,
#'              covariate = "age", 
#'              fill.str="age_type",
#'              by = "study")
#'              
#'# Example containing no fill.str, no overall average, but contains errorbars
#'
#'data.plot(data = diabetes.slr,
#'              covariate = "age", 
#'              half.length = "age_SD",
#'              avg.hline=FALSE,
#'              by = "study")

#' @export
data.plot <- function(
  data,
  covariate, 
  covariate.label = covariate,
  half.length = NULL,
  by = "study",
  avg.hline = TRUE,
  fill.str = NULL,
  text.size = 20,
  orientation = "landscape"
){
  
  if (!inherits(data, 'BUGSnetData'))
    stop("\'data\' must be a valid BUGSnetData object created using the data.prep function.")
  
  patients.data <- data$arm.data
  treatment.var <- data$varname.t
  trial.var <- data$varname.s
  axis.lab <- covariate.label
  
  if(is.null(half.length)){
    errorbar = FALSE
    errorbar.min = NULL
    errorbar.max = NULL
    caption <- ""
  } else{
    errorbar = TRUE
    errorbar.min <- paste0(covariate,"-",half.length)
    errorbar.max <- paste0(covariate,"+",half.length)
    caption <- paste0("Error bars: ", covariate," +/- ", half.length)
  }
  
  
  if (errorbar == TRUE && (is.null(errorbar.min) | is.null(errorbar.max))){
    stop('When errorbar = TRUE, errorbar.min and errorbar.max must be specified')
    }
  
  if (orientation == "landscape") {
  
    if (errorbar){eb = geom_errorbar(aes_string(ymin=errorbar.min, ymax=errorbar.max))}
    else {eb = NULL}
    
    if (avg.hline){
      overall.mean <- patients.data %>% select(covariate) %>% colMeans(na.rm=TRUE) %>% as.numeric()
      om <- geom_hline(yintercept = overall.mean, color = "red", linetype=2)}
    else{om=NULL}
    
    patients.data <- patients.data %>%
      mutate(trt =(!! as.name(treatment.var)), trial =(!! as.name(trial.var)))
    
    if (by == "study"){
    if(is.null(fill.str)) {p <- ggplot(patients.data, aes_string(x=treatment.var, y=covariate))}
    else p <- ggplot(patients.data, aes_string(x=treatment.var, y=covariate, color=fill.str))
    
    p <- p + geom_point() +
      eb +
      om + 
      theme(axis.text.x = element_text(angle = 90, hjust = 1),
            legend.title=element_blank(),
            plot.margin=unit(c(25, 5.5, 5.5, 50), "points"),
            strip.text.x = element_text(angle = 90),
            text = element_text(size=text.size)) +
      labs(y = axis.lab, x="", caption = caption)+
      facet_grid(. ~ trial,  space="free_x", scales="free_x")
    
    }
    else if (by == "treatment"){
      
    if(is.null(fill.str)) {p <- ggplot(patients.data, aes_string(x=trial.var, y=covariate))}
    else p <- ggplot(patients.data, aes_string(x=trial.var, y=covariate, color=fill.str))
    
    p <- p + geom_point() +
      eb +
      om + 
      theme(axis.text.x = element_text(angle = 90, hjust = 1),
            legend.title=element_blank(),
           plot.margin=unit(c(10, 4, 4, 20), "points"),
            strip.text.x = element_text(angle = 90),
           text = element_text(size=text.size)) +
      labs(y = axis.lab, x="", caption = caption)+
      facet_grid(. ~ trt,  space="free_x", scales="free_x")
    
    }
    
  } else if (orientation == "portrait") {
    
    if (errorbar){eb = geom_errorbar(aes_string(xmin=errorbar.min, xmax=errorbar.max))}
    else {eb = NULL}
    
    if (avg.hline){
      overall.mean <- patients.data %>% select(covariate) %>% colMeans(na.rm=TRUE) %>% as.numeric()
      om <- geom_vline(xintercept = overall.mean, color = "red", linetype=2)}
    else{om=NULL}
    
    patients.data <- patients.data %>%
      mutate(trt =(!! as.name(treatment.var)), trial =(!! as.name(trial.var)))
    
    if (by == "study"){
      if(is.null(fill.str)) {p <- ggplot(patients.data, aes_string(x=covariate, y=treatment.var))}
      else p <- ggplot(patients.data, aes_string(x=covariate, y=treatment.var, color=fill.str))
      
      p <- p + geom_point() +
        eb +
        om + 
        theme(axis.text.x = element_text(angle = 90, hjust = 1),
              legend.title=element_blank(),
              # plot.margin=unit(c(25, 5.5, 5.5, 50), "points"),
              plot.margin=unit(c(25, 5.5, 50, 5.5), "points"),
              strip.text.y = element_text(angle = 90),
              text = element_text(size=text.size)) +
        labs(x = axis.lab, y="", caption = caption)+
        facet_grid(trial ~ .,  space="free_y", scales="free_y") +
        theme(strip.text.y = element_text(angle = 0)) # Rotates facet label
      
    }
    else if (by == "treatment"){
      
      if(is.null(fill.str)) {p <- ggplot(patients.data, aes_string(x=covariate, y=trial.var))}
      else p <- ggplot(patients.data, aes_string(x=covariate, y=trial.var, color=fill.str))
      
      p <- p + geom_point() +
        eb +
        om + 
        theme(axis.text.x = element_text(angle = 90, hjust = 1),
              legend.title=element_blank(),
              plot.margin=unit(c(10, 4, 20, 4), "points"),
              strip.text.y = element_text(angle = 90),
              text = element_text(size=text.size)) +
        labs(x = axis.lab, y="", caption = caption)+
        facet_grid(trt ~ .,  space="free_y", scales="free_y") +
        theme(strip.text.y = element_text(angle = 0)) # Rotates facet label
      
    }
  }
  p
}
