#' Create Bugs Model for Arm-Based data
#' @description Creates BUGS code which can be ran through \code{nma.run()}.
#' 
#' @param data_arm A \code{BUGSnetData} object containing the data from arm-based trials produced by \code{data.prep()}
#' @param outcome A string indicating the name of your outcome variable for arm-based studies.
#' @param differences A string indicating the name of the differences for contrast-based studies
#' @param N A string indicating the name of the variable containing the number of participants in each arm for arm-based data - required.
#' @param sd.a A string (only required for continuous outcomes with arm-level data) indicating variable name
#' of the standard deviation of the outcome. Standard errors should be converted to standard deviation by multiplying by the square root of the sample size prior to using this function.
#' @param reference A string for the treatment that will be seen as the 'referent' comparator and labeled as treatment 1 in the BUGS code. This is often
#' a placebo or control drug of some kind.  
#' @param family A string indicating the family of the distribution of the outcome for arm-based trials. Options are:
#' "binomial", "normal", "poisson" 
#' @param link The link function for the nma model for arm-based models. Options are "logit" (binomial family), "log" (binomial family), "cloglog" (poisson family), "identity" (normal family).
#' @param time A string (only required for binomial-cloglog or poisson-log models) indicating the name of variable 
#'   indicating person-time followup (e.g person years) or study followup time.
#' @param effects A string indicating the type of treatment effect relative to baseline. Options are "fixed" or "random".
#' @param prior.mu A string of BUGS code that defines priors on the baseline treatment effects. By default, independent normal priors are used with mean 0 and standard deviation 15u, where u is the largest maximum likelihood estimator in single trials \insertCite{@see @gemtc}{BUGSnet}.
#' @param prior.d A string of BUGS code that defines define priors on relative treatment effects. By default, independent normal priors are used with mean 0 and standard deviation 15u, where u is the largest maximum likelihood estimator in single trials \insertCite{@see @gemtc}{BUGSnet}.
#' @param prior.sigma A string of BUGS code that defines the prior on the variance of relative treatment effects. By default, a uniform distribution with range 0 to u is used, where u is the largest maximum likelihood estimator in single trials \insertCite{@see @gemtc}{BUGSnet}.
#' @param prior.beta Optional string that defines the prior on the meta-regression coefficients. Options are "UNRELATED", "EXCHANGEABLE", "EQUAL" \insertCite{@TSD3}{BUGSnet} or a string of BUGS code.
#' @param covariate Optional string indicating the name of the variable in your data set that you would like to
#' adjust for via meta regression - only implemented for arm-based models. By default, covariate=NULL and no covariate adjustment is applied. If the specified covariate is numeric then
#' it will be centered for the analysis. If it is a character or factor then it will be treated as categorical. Currently only categorical variables
#' with fewer than 3 levels are supported. The name of the covariate variable must be the same in the arm and contrast-based data (if both are used)
#' @param type If type="inconsistency", an inconsistency model will be built. By default, type="consistency" and a consistency model is built.
#' will be built.
#' @return \code{nma.model} returns an object of class \code{BUGSnetModel} which is a list containing the following components:
#' @return \code{bugs} - A long character string containing BUGS code that will be run in \code{jags}.
#' @return \code{data} - The data used in the BUGS code.
#' @return \code{scale} - The scale of the outcome, based on the chosen family and link function
#' examples are "Risk Ratio" (relative risk), "Odds Ratio", "Mean Difference", "Hazard Ratio"
#' @return \code{trt.key} - Treatments mapped to integer numbers, used to run BUGS code.
#' @return ...
#' 
#' @details 
#' For meta-regression, the prespecified prior choices for the regression coefficients \eqn{\beta_{(1,2)},…,\beta_{(1,K)}} are
#' \describe{
#'   \item{Unrelated:}{\deqn{iid t(0, u^2, 1)}}
#'   \item{Exchangeable:}{\deqn{iid N(b, \gamma^2), b ~ t(0, u^2, 1), \gamma ~ U(0,u)}}
#'   \item{Equal:}{\deqn{\beta_2=...=\beta_T=B, B ~ t(0, u^2, 1)}}
#' }
#' where \eqn{u} is the largest maximum likelihood estimator in single trials \insertCite{@see @gemtc}{BUGSnet}.
#' 
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


nma.model <- function(data_arm = NULL,
                          outcome,
                          N = NULL,
                          sd.a=NULL,
                          reference,
                          type="consistency",
                          time=NULL,
                          family = NULL,
                          link = NULL,
                          effects,
                          prior.mu = "DEFAULT",
                          prior.d = "DEFAULT",
                          prior.sigma = "DEFAULT",
                          prior.beta = NULL,
                          covariate = NULL){
  
  arm <- TRUE
  contrast <- FALSE
  # 
  # if(is.null(data_arm)) {arm <- F}
  # if(is.null(data_contrast)) {contrast <- F}
  
  if((!is.null(data_arm) && class(data_arm) != "BUGSnetData"))
    stop("\'data_arm\' must be a valid BUGSnetData object created using the data.prep function.")
  
  if(!is.null(covariate) & is.null(prior.beta))stop("prior.beta must be specified when covariate is specified")
  if(is.null(covariate) & !is.null(prior.beta))stop("covariate must be specified when prior.beta is specified")
  if(!is.null(prior.beta)){
    if(!(prior.beta %in% c("UNRELATED","EQUAL","EXCHANGEABLE"))){
      stop("prior.beta must be either UNRELATED, EQUAL, or EXCHANGEABLE")
    }
  }
  
  # Warnings for different families and data requirements
  if(family=="normal" & is.null(sd.a)) stop("sd.a must be specified for continuous outcomes")
  if(family=="normal" & !(link%in% c("identity"))) stop("This combination of family and link is currently not supported in BUGSnet.")
  if(family=="poisson" & link!="log") stop("This combination of family and link is currently not supported in BUGSnet.")
  if(family=="binomial" & !(link %in% c("log","logit", "cloglog"))) stop("This combination of family and link is currently not supported in BUGSnet.")
  
  # Set up measurement scale
  if(link=="logit" & family %in% c("binomial", "binary", "bin", "binom")){
    scale <- "Odds Ratio"
  }else if(link=="log" & family %in% c("binomial", "binary", "bin", "binom")){
    scale <- "Risk Ratio"
  }else if(link== "identity" & family =="normal"){
    scale <- "Mean Difference"
  }else if(link =="cloglog" & family %in% c("binomial", "binary", "bin", "binom")){
    if(is.null(time)) stop("time must be specified when using a binomial family with the cloglog link")
    scale <- "Hazard Ratio"
  }else if(link == "log" & family =="poisson"){
    if(is.null(time)) stop("time must be specified when using a poisson family with the log link")
    scale <- "Rate Ratio"
  } 
  
  #pull relevant fields from the data and apply naming convention
  avarlist <- c(trt = data_arm$varname.t, trial = data_arm$varname.s, r1 = outcome, N = N, sd.a = sd.a, timevar = time, covariate = covariate) #se.diffs = se.diffs, var.ref = var.ref
  adata <- data_arm$arm.data[, avarlist]
  names(adata) <- names(avarlist)
  trts <- adata$trt
  
  if(!(reference %in% trts)) {
    
    stop("Reference is not present in network - pick a reference treatment that is in the network.")
    
  }
  
  #Trt mapping
  trt.key <- trts %>% unique %>% sort %>% tibble(trt.ini=.) %>%
    filter(trt.ini!=reference) %>% add_row(trt.ini=reference, .before=1) %>%
    mutate(trt.jags = 1:dim(.)[1])
  
  atrt <- trt.key[trt.key$trt.ini %in% adata$trt,]
  #add treatment mapping to data
  adata %<>% mutate(trt.jags=mapvalues(trt,
                                       from=atrt$trt.ini,
                                       to=atrt$trt.jags) %>% as.integer)
  
  # TODO test of this implementation works for arm and cb data
  #pre-process the covariate if specified
  if (!is.null(covariate)) {
    covariates <- adata$covariate
    if (is.numeric(covariates) == TRUE){
      #issue warning if covariate appears to be categorical
      if (length(unique(covariates)) < 5) {
        warning(paste0("The specified covariate is being treated as continuous. Ignore this warning if this is the intended behaviour. ",
                       "For the covariate to be treated as a categorical variable it must be converted to a factor in the data that is ",
                       "passed to data.prep."))
      }
      
      #de-mean covariate if continuous
      mean.cov <- mean(covariates, na.rm=TRUE)
      adata$covariate <- adata$covariate-mean.cov
    } else if (is.factor(covariates) == TRUE || is.character(covariates) == TRUE) {
      #check that covariate has fewer than 3 levels and convert strings and factors to binary covariates
      if (length(unique(covariates)) > 2)
        stop("BUGSnet does not currently support meta-regression with categorical variables that have more than two levels.")
      if (is.character(covariates) == TRUE) {
        adata$covariate <- factor(adata$covariate, levels = unique(covariates))
        adata$covariate <- as.numeric(adata$covariate != levels(adata$covariate)[1])
      }
    }
    
  } else{mean.cov <- NULL}
  
  # determine number of treatments in dataset
  nt <- length(unique(trts))
  
  #generate BUGS data object for arm-based data
  
  bugstemp <- adata %>% arrange(trial, trt.jags) %>% group_by(trial) %>% mutate(arm = row_number()) %>%
    ungroup() %>% select(-trt) %>% gather("variable", "value", -trial, -arm) %>% spread(arm, value)
  bugsdata2 <- list()
  for (v in unique(bugstemp$variable))
    bugsdata2[[v]] <- as.matrix(bugstemp %>% filter(variable == v) %>% select(-trial, -variable))
  #modify BUGS arm object for the various family/link combinations
  names(bugsdata2)[names(bugsdata2) == "trt.jags"] <- "t_a"
  names(bugsdata2)[names(bugsdata2) == "N"] <- "n"
  names(bugsdata2)[names(bugsdata2) == "covariate"] <- "x_a"
  if (family == "binomial" && link %in% c("log","logit")){
    names(bugsdata2)[names(bugsdata2) == "r1"] <- "r"
    bugsdata2 <- bugsdata2[names(bugsdata2) %in% c("ns_a", "nt", "na_a", "r", "n", "t_a", "x_a")]
  } else if (family == "normal" && link=="identity"){
    names(bugsdata2)[names(bugsdata2) == "r1"] <- "y"
    bugsdata2$se <- bugsdata2$sd / sqrt(bugsdata2$n)
    bugsdata2 <- bugsdata2[names(bugsdata2) %in% c("ns_a", "nt", "na_a", "y", "se", "t_a", "x_a")]
  } else if (family == "poisson" && link == "log"){
    names(bugsdata2)[names(bugsdata2) == "r1"] <- "r"
    names(bugsdata2)[names(bugsdata2) == "timevar"] <- "E"
    bugsdata2 <- bugsdata2[names(bugsdata2) %in% c("ns_a", "nt", "na_a", "r", "E", "t_a", "x_a")]
  } else if (family == "binomial" && link == "cloglog"){
    names(bugsdata2)[names(bugsdata2) == "r1"] <- "r"
    names(bugsdata2)[names(bugsdata2) == "timevar"] <- "time"
    bugsdata2 <- bugsdata2[names(bugsdata2) %in% c("ns_a", "nt", "na_a", "r", "n", "t_a", "x_a", "time")]
  }
  
  bugsdata2$ns_a <- data_arm$studies %>% nrow()
  bugsdata2$na_a <- data_arm$n.arms %>% select(n.arms) %>% t() %>% as.vector
  
  # # generate BUGS data object for cb data
  # 
  # if(contrast) {
  #   
  #   bugstemp2 <- cdata %>% arrange(trial) %>% group_by(trial) %>% mutate(arm = row_number()) %>% # changed
  #     ungroup() %>% select(-trt) %>% gather("variable", "value", -trial, -arm) %>% spread(arm, value)
  #   bugsdata3 <- list()
  #   
  #   for (v in unique(bugstemp2$variable))
  #     bugsdata3[[v]] <- as.matrix(bugstemp2 %>% filter(variable == v) %>% select(-trial, -variable))
  #   
  #   # modify BUGS contrast object and add treatment index, response, covariate, studies, number of arms to data
  #   
  #   names(bugsdata3)[names(bugsdata3) == "trt.jags"] <- "t_c"  
  #   names(bugsdata3)[names(bugsdata3) == "r1"] <- "y_c"
  #   names(bugsdata3)[names(bugsdata3) == "covariate"] <- "x_c"
  #   bugsdata3 <- bugsdata3[names(bugsdata3) %in% c("ns_c", "nt", "na_c", "y_c", "se.diffs", "var.ref", "t_c", "x_c")]
  #   
  #   bugsdata3$ns_c <- data_contrast$studies %>% nrow()
  #   bugsdata3$na_c <- data_contrast$n.arms %>% select(n.arms) %>% t() %>% as.vector
  #   
  #   # Data checking for contrast input
  #   
  #   # check that first arm difference is NA, change them to 0
  #   if(!all(is.na(bugsdata3$y_c[,1]))) {
  #     
  #     stop("The response in the first arm of each contrast-based trial should be NA.")
  #     
  #   } else {
  #     #convert to 0's
  #     bugsdata3$y_c[,1] <- 0
  #   }
  #   
  #   # check that first arm se is NA
  #   if(!all(is.na(bugsdata3$se.diffs[,1]))) {
  #     
  #     stop("The standard errors in the first arm of each contrast-based trial should be NA.")
  #     
  #   }
  #   
  #   # check if there are multi-arm trials, and if there are, check that var.ref is specified for the first arm
  #   if(!all(bugsdata3$na_c == 2)) {
  #     
  #     message("there are multi-arm trials")
  #     if(is.null(var.ref)) { 
  #       stop("var.ref must be specified if there are multi-arm contrast-based trials.") 
  #     } else {
  #       
  #       #check that only the first arm is specified
  #       if(!all(is.na(c(bugsdata3$var.ref[,-c(1)])))) {
  #         
  #         stop("Only the observed variances for the control arms for contrast-based trials should be included, all other arms should be NA")
  #         
  #       }  else if(!all(is.numeric(bugsdata3$var.ref[bugsdata3$na_c >2,1]))) { # make sure the control arms for all multi arm trials is numeric
  #         
  #         stop("Control arm observed variances for multi-arm conrtast-based trials must be numeric")
  #         
  #         bugsdata3$var.ref <- matrix(0, nrow = bugsdata3$ns_c, ncol = 2)
  #         
  #       }
  #       
  #     }
  #   } else {
  #     
  #     if(!is.null(var.ref)) { # if var.ref is specified when there are no multi-armed trials, set them all to 0 and print warning
  #       
  #       message("Control arm variances are not used for contrast-based trials with two arms.")
  #       
  #     }
  #     
  #     bugsdata3$var.ref <- matrix(0, nrow = bugsdata3$ns_c, ncol = 2) # set to zero to avoid compilation error
  #     
  #   }
  #   
  # } else {bugsdata3 <- data.frame()}
  
  
  
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
  
  max.delta <- paste0(nma.prior(data_arm, data_contrast=NULL, outcome=outcome, differences = differences, scale=scale, N=N, sd=sd.a, time = time))
  
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
  
  #meta regression string
  if (!is.null(covariate)){
    
    if (prior.beta=="UNRELATED"){
      prior.meta.reg <- sprintf("beta[1]<-0
    for (k in 2:nt){
      beta[k] ~ dt(0, (%s)^(-2), 1)
    }", max.delta)
    }else if(prior.beta=="EXCHANGEABLE"){
      prior.meta.reg <- sprintf("beta[1]<-0
    for (k in 2:nt){
      beta[k] ~ dnorm(b, gamma^(-2))
    }
    b~dt(0, %s^(-2), 1)
    gamma~dunif(0, %s)", max.delta, max.delta)
    }else if(prior.beta=="EQUAL"){
      prior.meta.reg <- sprintf("beta[1]<-0
    for (k in 2:nt){
      beta[k] <- B
    }
    B~dt(0, %s^(-2), 1)", max.delta)
    }else {
      prior.meta.reg <- prior.beta
    }
  }else prior.meta.reg <- ""
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
                        covariate,
                        prior.meta.reg,
                        auto = FALSE, # for compatibility with auto-run function - can change this if the feature is added
                        arm = arm,
                        contrast = contrast) %>%
    paste0(add.to.model)
  
  # if(!arm) {
  #   
  #   if(effects == "random") { # for random effects, need to specify ns_a=0
  #     
  #     bugsdata <- c(bugsdata3, nt=nt, ns_a = 0)
  #     
  #   } else {
  #     
  #     bugsdata <- c(bugsdata3, nt=nt) # fixed effects - don't specify ns_a=0 to avoid warning about unused variable from JAGS
  #     
  #   }
  #   
  #   
  # } else { # otherwise (just arm or both arm and contrast)
  #   
  #   bugsdata <- c(bugsdata2,bugsdata3, nt = nt)
  #   
  # }
  
  bugsdata <- c(bugsdata2, nt=nt)
  
  
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
                           time=time,
                           outcome=outcome,
                           N=N,
                           sd=sd,
                           mean.cov=mean.cov),
                      class = "BUGSnetModel")
  return(bmodel)
}

