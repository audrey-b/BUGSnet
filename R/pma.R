#no compatibility for data type: event rates with different study follow up times

#' Pairwise meta-analysis
#' @description implements pairwise meta-analysis via the package \code{meta}
#' @param data A BUGSnetData object produced by \code{data.prep()}
#' @param name.trt1 A string indicating the name of the comparator treatment (often Placebo)
#' @param name.trt2 A string indicating the name of the experimental treatment
#' @param outcome A string indicating the name of your outcome variable
#' @param N A string indicating the name of the variable containing the number of participants in each arm
#' @param type.outcome A string. Options are: "binomial", "continuous", "rate" (e.g # of events and # person-years reported)
#' @param sm A character string indicating which summary measure ("RR", "OR", "RD", or
#' ASD") is to be used for pooling of studies.
#' @param sd A string (only required when type.outcome="continuous") indicating variable name
#' of the standard deviation of the outcome
#' @param time A string required when type.outcome = "rate". Name of variable indicating person-time followup (e.g person years).
#' @param method A character string indicating what type of test was used. For more info, see \code{meta}'s
#' documentation.
#' @param method.tau A character string indicating which method is used to estimate the between study variance. 
#' Either "DL", "PM", "REML", "ML", "HS", "SJ", "HE", or "EB", can be abbreviated.
#' 
#' @return A forest plot as produced by the package \code{meta}
#' @return \code{raw} - dataset containing summary statistics of meta-analysis (effect estimates,
#'  confidence bounds, I-squared, Q-statistic)
#'  
#' @note Depending on the the value of \code{type.outcome}, this function will implement the functions
#' \code{metabin} (dichotomous outcomes), \code{metacont} (continuos outcomes), or \code{metainc} (rate outcomes)
#'  from the package \code{meta}. 
#' 
#' 
#' @examples
#' 
#' data(diabetes.sim)
#' 
#' diabetes.slr <- data.prep(arm.data = diabetes.sim, 
#' varname.t = "Treatment", 
#' varname.s = "Study")
#' 
#' pma(data = diabetes.slr,
#' type.outcome="binomial",
#' sm="OR",
#' name.trt1 = "Placebo", 
#' name.trt2 = "Diuretic", 
#' outcome = "diabetes",
#' N = "n")
#'                            
#' @export
#' @seealso \code{\link{data.prep}}

pma <- function(data,
                     name.trt1, 
                     name.trt2, 
                     outcome,
                     N,
                     sd = NULL,
                     time = NULL,
                     type.outcome,
                     method = "MH",
                     method.tau="DL",
                     sm){
  
  #Bind variables to function
  comparison <- NULL
  trt.e <- NULL
  trt.c <- NULL
  
  
  if(class(data) != "BUGSnetData")
    stop("\'data\' must be a valid BUGSnetData object created using the data.prep function.")
  
  if(type.outcome=="continuous" & is.null(sd)) stop("sd must be specified for continuous outcomes")
  
  data.nma <- data
  
  # Check if number of events and participants is integer
  tmp.check1 <- data.nma$arm.data %>% select(outcome)
  tmp.check2 <- data.nma$arm.data %>% select(N)
  
  if(type.outcome %in% c("binomial","rate") && all(tmp.check1%%1!=0)) {
    stop('The "outcome" variable (number of events) must be an integer')
  } else if(all(tmp.check2%%1!=0)) {
    stop('The "N" variable (number of participants) must be an integer')
  }
  if(type.outcome=="continuous" && !is.null(method) && method !="Inverse"){
    warning('For meta-analysis with continuous outcomes, inverse variance 
            weighting was used for pooling.')
  }
  
  if(type.outcome %in% c("bin","binom","binomial","binary")){ 
    pairwise.dat <- by.comparison(data.nma=data.nma, outcome=outcome, type.outcome="binomial", N=N)
    
    pairwise.dat.with.c <- pairwise.dat %>% 
      filter(grepl(name.trt1, comparison)) %>%
      filter((trt.e == name.trt2 | trt.c == name.trt2))
    
    #names(pairwise.dat.with.c)[names(pairwise.dat.with.c) == data.nma$varname.s] <- "study.id"
    
    if(dim(pairwise.dat.with.c)[1]!=0){
      meta1 <- metabin(pairwise.dat.with.c[,paste0(outcome,".e")] %>% t() %>% as.vector,
                       pairwise.dat.with.c[,paste0(N,".e")] %>% t() %>% as.vector,
                       pairwise.dat.with.c[,paste0(outcome,".c")] %>% t() %>% as.vector,
                       pairwise.dat.with.c[,paste0(N,".c")] %>% t() %>% as.vector,
                       method = method,
                       method.tau=method.tau,
                       sm=sm)
      tmp.tab <- summary(meta1)
      
      est <- data.frame(
        comparison = paste(name.trt1, " vs. ", name.trt2, sep=""),
        n.studies = tmp.tab$k,
        i.squared = tmp.tab$I2$TE,
        re.estimate = tmp.tab$random$TE %>% exp,
        re.lci = tmp.tab$random$lower %>% exp,
        re.uci = tmp.tab$random$upper %>% exp,
        fe.estimate = tmp.tab$fixed$TE %>% exp,
        fe.lci = tmp.tab$fixed$lower %>% exp,
        fe.uci = tmp.tab$fixed$upper %>% exp)
      
      return(list("summary" = est, "forest" = forest(meta1, studlab = pairwise.dat.with.c$trial),"raw"=meta1))
      
    }
    else(return("No direct information for this comparison"))
    
  } else if (type.outcome %in% c("cont", "continuous")){
    pairwise.dat <- by.comparison(data.nma=data.nma, outcome=outcome, type.outcome="continuous", N=N, sd=sd)
    
    pairwise.dat.with.c <- pairwise.dat %>% 
      filter(grepl(name.trt1, comparison)) %>%
      filter((trt.e == name.trt2 | trt.c == name.trt2))
    
    #names(pairwise.dat.with.c)[names(pairwise.dat.with.c) == data.nma$varname.s] <- "study.id"
    
    if(dim(pairwise.dat.with.c)[1]!=0){
      meta1 <- metacont(pairwise.dat.with.c[,paste0(N,".e")] %>% t() %>% as.vector,
                        pairwise.dat.with.c[,paste0(outcome,".e")] %>% t() %>% as.vector,
                        pairwise.dat.with.c[,paste0(sd,".e")] %>% t() %>% as.vector,
                        pairwise.dat.with.c[,paste0(N,".c")] %>% t() %>% as.vector,
                        pairwise.dat.with.c[,paste0(outcome,".c")] %>% t() %>% as.vector,
                        pairwise.dat.with.c[,paste0(sd,".c")] %>% t() %>% as.vector,
                        method.tau=method.tau,
                        sm=sm)
      tmp.tab <- summary(meta1)
      
      est <- data.frame(
        comparison = paste(name.trt1, " vs. ", name.trt2, sep=""),
        n.studies = tmp.tab$k,
        i.squared = tmp.tab$I2$TE,
        re.estimate = tmp.tab$random$TE %>% exp,
        re.lci = tmp.tab$random$lower %>% exp,
        re.uci = tmp.tab$random$upper %>% exp,
        fe.estimate = tmp.tab$fixed$TE %>% exp,
        fe.lci = tmp.tab$fixed$lower %>% exp,
        fe.uci = tmp.tab$fixed$upper %>% exp)
      
      return(list("summary" = est, "forest" = forest(meta1, studlab = pairwise.dat.with.c$trial),"raw"=meta1))
      
    }
    else(return("No direct information for this comparison"))
    
  } else if (type.outcome %in% c("rate")){
    pairwise.dat <- by.comparison(data.nma=data.nma, outcome=outcome, type.outcome="rate", N=N, time=time)
    
    pairwise.dat.with.c <- pairwise.dat %>% 
      filter(grepl(name.trt1, comparison)) %>%
      filter((trt.e == name.trt2 | trt.c == name.trt2))
    
    #names(pairwise.dat.with.c)[names(pairwise.dat.with.c) == data.nma$varname.s] <- "study.id"
    
    if(dim(pairwise.dat.with.c)[1]!=0){
      meta1 <- metainc(event.e = pairwise.dat.with.c[,paste0(outcome,".e")] %>% t() %>% as.vector,
                       time.e = pairwise.dat.with.c[,paste0(time,".e")] %>% t() %>% as.vector,
                       event.c= pairwise.dat.with.c[,paste0(outcome,".c")] %>% t() %>% as.vector,
                       time.c = pairwise.dat.with.c[,paste0(time,".c")] %>% t() %>% as.vector,
                       method.tau=method.tau,
                       sm=sm)
      tmp.tab <- summary(meta1)
      
      est <- data.frame(
        comparison = paste(name.trt1, " vs. ", name.trt2, sep=""),
        n.studies = tmp.tab$k,
        i.squared = tmp.tab$I2$TE,
        re.estimate = tmp.tab$random$TE %>% exp,
        re.lci = tmp.tab$random$lower %>% exp,
        re.uci = tmp.tab$random$upper %>% exp,
        fe.estimate = tmp.tab$fixed$TE %>% exp,
        fe.lci = tmp.tab$fixed$lower %>% exp,
        fe.uci = tmp.tab$fixed$upper %>% exp)
      
      return(list("summary" = est, "forest" = forest(meta1),"raw"=meta1))
      
    }
    else(return("No direct information for this comparison"))
  }
  
}

