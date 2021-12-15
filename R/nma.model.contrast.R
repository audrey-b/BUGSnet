#' Create Bugs Model for contrast-Level data
#' @description Creates BUGS code which can be ran through \code{nma.run()}.
#' 
#' @param data_contrast A \code{BUGSnetData} object containing the data from contrast-based trials produced by \code{data.prep()}
#' @param differences A string indicating the name of the differences for contrast-based studies
#' @param se.diffs A string indicating the variable name of the standard errors of the differences
#' @param var.arm1 A string (only required for networks with multi-arm trials) indicating the variable name of the variance of the treatment in arm 1 of each study
#' @param reference A string for the treatment that will be seen as the 'referent' comparator and labeled as treatment 1 in the BUGS code. This is often
#' a placebo or control drug of some kind.  
#' @param effects A string indicating the type of treatment effect relative to baseline. Options are "fixed" or "random".
#' @param prior.mu A string of BUGS code that defines priors on the baseline treatment effects. By default, independent normal priors are used with mean 0 and standard deviation 15u, where u is the largest maximum likelihood estimator in single trials \insertCite{@see @gemtc}{BUGSnet}.
#' @param prior.d A string of BUGS code that defines define priors on relative treatment effects. By default, independent normal priors are used with mean 0 and standard deviation 15u, where u is the largest maximum likelihood estimator in single trials \insertCite{@see @gemtc}{BUGSnet}.
#' @param prior.sigma A string of BUGS code that defines the prior on the variance of relative treatment effects. By default, a uniform distribution with range 0 to u is used, where u is the largest maximum likelihood estimator in single trials \insertCite{@see @gemtc}{BUGSnet}.
#' @param type If type="inconsistency", an inconsistency model will be built. By default, type="consistency" and a consistency model is built.
#' will be built.
#' @param scale A string indicating the scale of the data, such as "Mean Difference" or "Log-Odds Ratio",
#' @return \code{nma.model} returns an object of class \code{BUGSnetModel} which is a list containing the following components:
#' @return \code{bugs} - A long character string containing BUGS code that will be run in \code{jags}.
#' @return \code{data} - The data used in the BUGS code.
#' @return \code{scale} - The scale of the outcome, based on the chosen family and link function
#' examples are "Risk Ratio" (relative risk), "Odds Ratio", "Mean Difference", "Hazard Ratio"
#' @return \code{trt.key} - Treatments mapped to integer numbers, used to run BUGS code.
#' @return ...
#' 
#' @examples
#' 
#' 
#' data(diabetes.sim)
#' 
#' diabetes.slr <- data.prep(arm.data = diabetes.sim, 
#' varname.t = "Treatment", 
#' varname.s = "Study")
#' 
#' #Random effects, consistency model.
#' #Binomial family, cloglog link. This implies that the scale will be the Hazard Ratio.
#'diabetes.re.c <- nma.model(data = diabetes.slr,
#'        outcome = "diabetes", 
#'        N = "n",
#'        reference = "Placebo",
#'        family = "binomial",
#'        link = "cloglog",
#'        effects = "random",
#'        type="consistency",
#'        time="followup"
#'        )
#'        
#' #Fixed effects, consistency model.
#' #Binomial family, cloglog link. This implies that the scale will be the Hazard Ratio.
#'diabetes.fe.c <- nma.model(data = diabetes.slr,
#'        outcome = "diabetes", 
#'        N = "n",
#'        reference = "Placebo",
#'        family = "binomial",
#'        link = "cloglog",
#'        effects = "fixed",
#'        type="consistency",
#'        time="followup"
#'        )
#'         
#' @export
#'
#' @references
#' \insertRef{gemtc}{BUGSnet}
#' 
#' \insertRef{TSD3}{BUGSnet}
#' 
#' 
#' @seealso \code{\link{data.prep}}, \code{\link{nma.run}}


nma.model.contrast <- function(data_contrast = NULL,
                      differences,
                      se.diffs,
                      var.arm1 = NULL,
                      reference,
                      type="consistency",
                      effects,
                      scale,
                      prior.mu = "DEFAULT",
                      prior.d = "DEFAULT",
                      prior.sigma = "DEFAULT"){
  
  contrast <- TRUE
  time <- NULL
  covariate <- NULL
  
  # Bind variables to function
  trt.ini <- NULL
  trt <- NULL
  trial <- NULL
  trt.jags <- NULL
  arm <- NULL
  value <- NULL
  variable <- NULL
  n.arms <- NULL
  
  # if(is.null(data_arm)) {arm <- F}
  # if(is.null(data_contrast)) {contrast <- F}
  
  if(!is.null(data_contrast) && class(data_contrast) != "BUGSnetData")
    stop("\'data_contrast\' must be a valid BUGSnetData object created using the data.prep function.")
  

    # set scale and outcome to dummies
    family <- ""
    link <- "identity"
    outcome <- "none"

  
  #pull relevant fields from the data and apply naming convention
  cvarlist <- c(trt = data_contrast$varname.t, trial = data_contrast$varname.s, r1 = differences, se.diffs = se.diffs, var.arm1 = var.arm1, covariate = NULL) #se.diffs = se.diffs, var.arm1 = var.arm1
  cdata <- data_contrast$arm.data[, cvarlist]
  names(cdata) <- names(cvarlist)
  
  trts <- cdata$trt
  
  if(!(reference %in% trts)) {
    
    stop("Reference is not present in network - pick a reference treatment that is in the network.")
    
  }
  
  #Trt mapping
  trt.key <- trts %>% unique %>% sort %>% tibble(trt.ini=.) %>%
    filter(trt.ini!=reference) %>% add_row(trt.ini=reference, .before=1) %>%
    mutate(trt.jags = 1:dim(.)[1])
  
  ctrt <- trt.key[trt.key$trt.ini %in% cdata$trt,]
  #add treatment mapping to data
  cdata %<>% mutate(trt.jags=mapvalues(trt,
                                       from=ctrt$trt.ini,
                                       to=ctrt$trt.jags) %>% as.integer)
  
  mean.cov <- NULL
  
  # determine number of treatments in dataset
  nt <- length(unique(trts))
  
  # generate BUGS data object for cb data
  bugstemp2 <- cdata %>% arrange(trial) %>% group_by(trial) %>% mutate(arm = row_number()) %>% # changed
    ungroup() %>% select(-trt) %>% gather("variable", "value", -trial, -arm) %>% spread(arm, value)
  bugsdata3 <- list()
  
  for (v in unique(bugstemp2$variable))
    bugsdata3[[v]] <- as.matrix(bugstemp2 %>% filter(variable == v) %>% select(-trial, -variable))
  
  # modify BUGS contrast object and add treatment index, response, covariate, studies, number of arms to data
  names(bugsdata3)[names(bugsdata3) == "trt.jags"] <- "t_c"  
  names(bugsdata3)[names(bugsdata3) == "r1"] <- "y_c"
  names(bugsdata3)[names(bugsdata3) == "covariate"] <- "x_c"
  bugsdata3 <- bugsdata3[names(bugsdata3) %in% c("ns_c", "nt", "na_c", "y_c", "se.diffs", "var.arm1", "t_c", "x_c")]
  
  bugsdata3$ns_c <- data_contrast$studies %>% nrow()
  bugsdata3$na_c <- data_contrast$n.arms %>% select(n.arms) %>% t() %>% as.vector
  
  # Data checking for contrast input
  
  # check that first arm difference is NA, change them to 0
  check_nas <-bugsdata3$y_c[,1]
  
  if(!all(bugsdata3$y_c[,1] %in% c(NA, "NA"))) {
    
    stop("The response in the first arm of each contrast-based trial should be NA.")
    
  } else {
    #convert to 0's
    bugsdata3$y_c[,1] <- 0
  }
  
  # check that first arm se is NA
  if(!all(bugsdata3$se.diffs[,1] %in% c(NA, "NA"))) {
    
    stop("The standard errors in the first arm of each contrast-based trial should be NA.")
    
  }
  
  # Check that se's are positive
  if(!all(bugsdata3$se.diffs[,-c(1)] > 0, na.rm = T)) {
    
    stop("The standard errors of contrasts should be positive")
    
  }
  
  # check if there are multi-arm trials, and if there are, check that var.arm1 is specified for the first arm
  if(!all(bugsdata3$na_c == 2)) {
    
    if(is.null(var.arm1)) { 
      stop("var.arm1 must be specified if there are multi-arm contrast-based trials.") 
    } else {
      
      #check that only the first arm is specified
      if(!all(c(bugsdata3$var.arm1[,-c(1)]) %in% c(NA, "NA"), na.rm = T)) {
        
        stop("Only the observed variances for the control arms for contrast-based trials should be included, all other arms should be NA")
        
      }  else if(!all(is.numeric(bugsdata3$var.arm1[bugsdata3$na_c >2,1]))) { # make sure the control arms for all multi arm trials is numeric
        
        stop("Control arm observed variances for multi-arm contrast-based trials must be numeric")
        
      } else if (!all(bugsdata3$var.arm1[bugsdata3$na_c >2,1] > 0)) { # make sure variances are positiveF
        
        stop("COntrol arm observed variances for multi-arm contrast-based trials must be positive")
        
      }
      # set the variances to 0 in the first column to avoid compilation error
      bugsdata3$var.arm1[bugsdata3$na_c <3,1] <- 0
      
    }
      
  } else {
    
    if(!is.null(var.arm1)) { # if var.arm1 is specified when there are no multi-armed trials, set them all to 0 and print warning
      
      message("Control arm variances are not used for contrast-based trials with two arms.")
      
    }
    
    bugsdata3$var.arm1 <- matrix(0, nrow = bugsdata3$ns_c, ncol = 2) # set to zero to avoid compilation error
    
  }
  
  
  
  # make legend for treatment names and numbering in jags program
  
  add.to.model <- trt.key %>%
    transmute(Treatments = paste0("# ", trt.jags, ": ", trt.ini, "\n")) %>%
    t() %>%
    paste0() %>%
    paste(collapse="") %>%
    paste0("\n\n# Treatment key\n",.)
  
  ############
  ###Priors###
  ############

  max.delta <- paste0(nma.prior(data_arm=NULL, data_contrast = data_contrast, outcome=outcome,
                                differences = differences, scale=scale, N=NULL, sd=NULL, time = NULL))
  # BASELINE EFFECTS PRIOR
  if (prior.mu == "DEFAULT"){
    prior.mu.str <- sprintf("for(i in 1:ns_a){
                               mu[i] ~ dnorm(0,(%s*15)^(-2))
                               }", max.delta)
  } else {
    prior.mu.str <- sprintf("for(i in 1:ns_a){
                               mu[i] ~ %s
  }", prior.mu)
  }
  
  # RELATIVE EFFECTS PRIOR
  if (prior.d =="DEFAULT"){
    if(type=="consistency"){
      prior.d.str <- sprintf("for(k in 2:nt){
                             d[k] ~ dnorm(0,(%s*15)^(-2))
    }", max.delta)
      
    } else if(type=="inconsistency"){
      prior.d.str <- sprintf("for (c in 1:(nt-1)) {
                           for (k in (c+1):nt)  {
                           d[c,k] ~ dnorm(0,(%s*15)^(-2))
                           d[k,c] <- d[c,k]
                           }
    }", max.delta)
    }
  } else {
    if(type=="consistency"){
      prior.d.str <- sprintf("for(k in 2:nt){
                             d[k] ~ %s
    }", prior.d)
      
    }else if(type=="inconsistency"){
      prior.d.str <- sprintf("for (c in 1:(nt-1)) {
                           for (k in (c+1):nt)  {
                           d[c,k] ~ %s
                           d[k,c] <- d[c,k]
                           }
  }", prior.d)
    }
  }
  
  # RANDOM EFFECTS VARIANCE PRIOR
  if (prior.sigma == "DEFAULT"){
    prior.sigma2.str <-  sprintf("sigma ~ dunif(0,%s)
                                 sigma2 <- sigma^2", max.delta)
  } else {
    prior.sigma2.str <-  sprintf("sigma ~ %s
                                 sigma2 <- sigma^2", prior.sigma)
  }
  
  ###hot fix for binomial family with log link.
  if(scale == "Risk Ratio"){
    prior.mu.str <-  "for(i in 1:ns_a){
     mu[i] <- log(p.base[i])           #priors for all trial baselines
     p.base[i] ~ dunif(0, 1)
   }"
  }
  
  # # #meta regression string
  # if (!is.null(covariate)){
  # 
  #   if (prior.beta=="UNRELATED"){
  #     prior.meta.reg <- sprintf("beta[1]<-0
  #   for (k in 2:nt){
  #     beta[k] ~ dt(0, (%s)^(-2), 1)
  #   }", max.delta)
  #   }else if(prior.beta=="EXCHANGEABLE"){
  #     prior.meta.reg <- sprintf("beta[1]<-0
  #   for (k in 2:nt){
  #     beta[k] ~ dnorm(b, gamma^(-2))
  #   }
  #   b~dt(0, %s^(-2), 1)
  #   gamma~dunif(0, %s)", max.delta, max.delta)
  #   }else if(prior.beta=="EQUAL"){
  #     prior.meta.reg <- sprintf("beta[1]<-0
  #   for (k in 2:nt){
  #     beta[k] <- B
  #   }
  #   B~dt(0, %s^(-2), 1)", max.delta)
  #   }else {
  #     prior.meta.reg <- prior.beta
  #   }
  # }else prior.meta.reg <- ""
  
  prior.meta.reg <- ""
  
  # #remove covariate from bugsdata2 if unused
  # if (is.null(covariate)){bugsdata2 <- bugsdata2[names(bugsdata2)!="x"]}
  
  # make the code for the model
  model <- makeBUGScode(family=family,
                        link=link,
                        effects=effects,
                        inconsistency=(type=="inconsistency"),
                        prior.mu.str,
                        prior.d.str,
                        prior.sigma2.str,
                        meta.covariate = NULL,
                        prior.meta.reg,
                        auto = FALSE, # for compatibility with auto-run function - can change this if the feature is added
                        arm = FALSE,
                        contrast = contrast) %>%
    paste0(add.to.model)
  
  if(effects == "random") { # for random effects, need to specify ns_a=0
    
    bugsdata <- c(bugsdata3, nt=nt, ns_a = 0)
    
  } else {
    
    bugsdata <- c(bugsdata3, nt=nt) # fixed effects - don't specify ns_a=0 to avoid warning about unused variable from JAGS
    
  }
  
  
  bmodel <- structure(list(bugs=model,
                           data=bugsdata,
                           scale=scale,
                           trt.key=trt.key,
                           family=family,
                           link=link,
                           type=type,
                           effects=effects,
                           covariate=NULL,
                           prior.mu=prior.mu,
                           prior.d=prior.d,
                           prior.sigma=prior.sigma,
                           prior.beta=NULL,
                           reference=reference,
                           time=NULL,
                           outcome=outcome,
                           N=NULL,
                           sd=sd,
                           mean.cov=mean.cov,
                           study_c = data_contrast$studies),
                      class = "BUGSnetModel")
  return(bmodel)
}

