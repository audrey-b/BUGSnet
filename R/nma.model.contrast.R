#' Create Bugs Model for contrast-based data
#' @description Creates BUGS code which can be ran through \code{nma.run()}.
#' 
#' @param data_contrast A \code{BUGSnetData} object containing the data from contrast-based trials produced by \code{data.prep()}
#' @param differences A string indicating the name of the differences for contrast-based studies
#' @param se.c A string indicating the variable name of the standard errors of the differences.
#' @param var.ref A string (only required for networks with multi-arm trials) indicating the variable name of the variance of the reference treatment in each study
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
#' @details 
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
                      se.diffs = NULL,
                      var.ref = NULL,
                      reference,
                      type="consistency",
                      time=NULL,
                      effects,
                      scale,
                      prior.mu = "DEFAULT",
                      prior.d = "DEFAULT",
                      prior.sigma = "DEFAULT"){
  
  arm <- FALSE
  contrast <- TRUE
  # 
  # if(is.null(data_arm)) {arm <- F}
  # if(is.null(data_contrast)) {contrast <- F}
  
  if(!is.null(data_contrast) && class(data_contrast) != "BUGSnetData")
    stop("\'data_contrast\' must be a valid BUGSnetData object created using the data.prep function.")
  

    # set scale and outcome to dummies
    family <- ""
    link <- "identity"
    outcome <- "none"

  
  #pull relevant fields from the data and apply naming convention
  cvarlist <- c(trt = data_contrast$varname.t, trial = data_contrast$varname.s, r1 = differences, se.diffs = se.diffs, var.ref = var.ref, covariate = covariate) #se.diffs = se.diffs, var.ref = var.ref
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
  
  # # TODO test of this implementation works for arm and cb data
  # #pre-process the covariate if specified
  # if (!is.null(covariate)) {
  #   covariates <- c(adata$covariate, cdata$covariate)
  #   if (is.numeric(covariates) == TRUE){
  #     #issue warning if covariate appears to be categorical
  #     if (length(unique(covariates)) < 5) {
  #       warning(paste0("The specified covariate is being treated as continuous. Ignore this warning if this is the intended behaviour. ",
  #                      "For the covariate to be treated as a categorical variable it must be converted to a factor in the data that is ",
  #                      "passed to data.prep."))
  #     }
  #     
  #     #de-mean covariate if continuous
  #     mean.cov <- mean(covariates, na.rm=TRUE)
  #     adata$covariate <- adata$covariate-mean.cov
  #     cdata$covariate <- cdata$covariate-mean.cov
  #   } else if (is.factor(covariates) == TRUE || is.character(covariates) == TRUE) {
  #     #check that covariate has fewer than 3 levels and convert strings and factors to binary covariates
  #     if (length(unique(covariates)) > 2)
  #       stop("BUGSnet does not currently support meta-regression with categorical variables that have more than two levels.")
  #     if (is.character(covariates) == TRUE) {
  #       adata$covariate <- factor(adata$covariate, levels = unique(covariates))
  #       cdata$covariate <- factor(cdata$covariate, levels = unique(covariates))
  #       adata$covariate <- as.numeric(adata$covariate != levels(adata$covariate)[1])
  #       cdata$covariate <- as.numeric(cdata$covariate != levels(adata$covariate)[1])
  #     }
  #   }
  #   
  # } else{mean.cov <- NULL}
  mean.cov <- NULL
  
  # determine number of treatments in dataset
  nt <- length(unique(trts))
  
  #generate BUGS data object for arm-based data
  
  # if(arm) {
  #   
  #   bugstemp <- adata %>% arrange(trial, trt.jags) %>% group_by(trial) %>% mutate(arm = row_number()) %>%
  #     ungroup() %>% select(-trt) %>% gather("variable", "value", -trial, -arm) %>% spread(arm, value)
  #   bugsdata2 <- list()
  #   for (v in unique(bugstemp$variable))
  #     bugsdata2[[v]] <- as.matrix(bugstemp %>% filter(variable == v) %>% select(-trial, -variable))
  #   #modify BUGS arm object for the various family/link combinations
  #   names(bugsdata2)[names(bugsdata2) == "trt.jags"] <- "t_a"
  #   names(bugsdata2)[names(bugsdata2) == "N"] <- "n"
  #   names(bugsdata2)[names(bugsdata2) == "covariate"] <- "x_a"
  #   if (family == "binomial" && link %in% c("log","logit")){
  #     names(bugsdata2)[names(bugsdata2) == "r1"] <- "r"
  #     bugsdata2 <- bugsdata2[names(bugsdata2) %in% c("ns_a", "nt", "na_a", "r", "n", "t_a", "x_a")]
  #   } else if (family == "normal" && link=="identity"){
  #     names(bugsdata2)[names(bugsdata2) == "r1"] <- "y"
  #     bugsdata2$se <- bugsdata2$sd / sqrt(bugsdata2$n)
  #     bugsdata2 <- bugsdata2[names(bugsdata2) %in% c("ns_a", "nt", "na_a", "y", "se", "t_a", "x_a")]
  #   } else if (family == "poisson" && link == "log"){
  #     names(bugsdata2)[names(bugsdata2) == "r1"] <- "r"
  #     names(bugsdata2)[names(bugsdata2) == "timevar"] <- "E"
  #     bugsdata2 <- bugsdata2[names(bugsdata2) %in% c("ns_a", "nt", "na_a", "r", "E", "t_a", "x_a")]
  #   } else if (family == "binomial" && link == "cloglog"){
  #     names(bugsdata2)[names(bugsdata2) == "r1"] <- "r"
  #     names(bugsdata2)[names(bugsdata2) == "timevar"] <- "time"
  #     bugsdata2 <- bugsdata2[names(bugsdata2) %in% c("ns_a", "nt", "na_a", "r", "n", "t_a", "x_a", "time")]
  #   }
  #   
  #   bugsdata2$ns_a <- data_arm$studies %>% nrow()
  #   bugsdata2$na_a <- data_arm$n.arms %>% select(n.arms) %>% t() %>% as.vector
  #   
  # } else (bugsdata2 <- data.frame())
  
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
  bugsdata3 <- bugsdata3[names(bugsdata3) %in% c("ns_c", "nt", "na_c", "y_c", "se.diffs", "var.ref", "t_c", "x_c")]
  
  bugsdata3$ns_c <- data_contrast$studies %>% nrow()
  bugsdata3$na_c <- data_contrast$n.arms %>% select(n.arms) %>% t() %>% as.vector
  
  # Data checking for contrast input
  
  # check that first arm difference is NA, change them to 0
  if(!all(is.na(bugsdata3$y_c[,1]))) {
    
    stop("The response in the first arm of each contrast-based trial should be NA.")
    
  } else {
    #convert to 0's
    bugsdata3$y_c[,1] <- 0
  }
  
  # check that first arm se is NA
  if(!all(is.na(bugsdata3$se.diffs[,1]))) {
    
    stop("The standard errors in the first arm of each contrast-based trial should be NA.")
    
  }
  
  # check if there are multi-arm trials, and if there are, check that var.ref is specified for the first arm
  if(!all(bugsdata3$na_c == 2)) {
    
    if(is.null(var.ref)) { 
      stop("var.ref must be specified if there are multi-arm contrast-based trials.") 
    } else {
      
      #check that only the first arm is specified
      if(!all(is.na(c(bugsdata3$var.ref[,-c(1)])))) {
        
        stop("Only the observed variances for the control arms for contrast-based trials should be included, all other arms should be NA")
        
      }  else if(!all(is.numeric(bugsdata3$var.ref[bugsdata3$na_c >2,1]))) { # make sure the control arms for all multi arm trials is numeric
        
        stop("Control arm observed variances for multi-arm conrtast-based trials must be numeric")
        
        bugsdata3$var.ref <- matrix(0, nrow = bugsdata3$ns_c, ncol = 2)
        
      }
      
    }
  } else {
    
    if(!is.null(var.ref)) { # if var.ref is specified when there are no multi-armed trials, set them all to 0 and print warning
      
      message("Control arm variances are not used for contrast-based trials with two arms.")
      
    }
    
    bugsdata3$var.ref <- matrix(0, nrow = bugsdata3$ns_c, ncol = 2) # set to zero to avoid compilation error
    
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

  max.delta <- paste0(nma.prior(data_arm=NULL, data_contrast, outcome=outcome, differences = differences, scale=scale, N=N, sd=sd.a, time = NULL))
  
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
  
  # #meta regression string
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
  #FIXME
  model <- makeBUGScode(family=family,
                        link=link,
                        effects=effects,
                        inconsistency=(type=="inconsistency"),
                        prior.mu.str,
                        prior.d.str,
                        prior.sigma2.str,
                        covariate,
                        prior.meta.reg,
                        auto = FALSE, # for compatibility with auto-run function - can change this if the feature is added
                        arm = arm,
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
                           covariate=covariate,
                           prior.mu=prior.mu,
                           prior.d=prior.d,
                           prior.sigma=prior.sigma,
                           prior.beta=prior.beta,
                           reference=reference,
                           time=NULL,
                           outcome=outcome,
                           N=NULL,
                           sd=sd,
                           mean.cov=mean.cov),
                      class = "BUGSnetModel")
  return(bmodel)
}

