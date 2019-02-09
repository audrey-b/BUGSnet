#' Make Bugs Model
#' @description Creates bugs BUGS code which will can be ran through \code{nma.analysis()}.
#' 
#' @param data.prep A data object produced by \code{data.prep()}
#' @param outcome A string indicating the name of your outcome variable
#' @param N A string indicating the name of the variable containing the number of participants in each arm
#' @param reference A string for the treatment that will be seen as the 'referent' comparator and labeled as treatment 1 in the BUGS code. This is often
#' a placebo or control drug of some kind.  
#' @param sd A string (only required for continuous outcomes) indicating variable name
#' of the standard deviation of the outcome
#' @param type If type="inconsistency", an inconsistency model will be built. Otherwise a consistency model
#' will be built
#' @param exposure.time A string (only required for rate or rate2 outcomes) indicating the name of variable 
#'   indicating person-time followup (e.g person years) or study followup time.
#' @param family A string indicating the family of the distribution of the outcome. Options are:
#' "binomial", "normal", "poisson".
#' @param link The link function for the nma model. Options are "logit", "log", "cloglog", "identity".
#' @param effects A string indicating the type of treatment effect relative to baseline. Options are "fixed" or "random".
#' @param prior.mu A string of BUGS code to define priors on the baseline treatment effects. By default, independent normal priors are used with mean 0 and standard deviation 15u, where u is the largest maximum likelihood estimator in single trials (as suggested by ).
#' @param prior.d A string of BUGS code to define priors on relative treatment effects. By default, independent normal priors are used with mean 0 and standard deviation 15u, where u is the largest maximum likelihood estimator in single trials (as suggested by ).
#' @param prior.sigma A string of BUGS code for the priors on the variance of relative treatment effects. By default, a uniform distribution with range 0 to u is used, where u is the largest maximum likelihood estimator in single trials (as suggested by ).
#' @param meta.covariate Optional string indicating the name of the variable in your data set that you would like to
#' adjust for via meta regression.
#' 
#' @return \code{model.str} - A long character string containing BUGS code that will be run in \code{rjags}.
#' @return \code{bugsdata2} - Data accompanying BUGS code.
#' @return \code{scale} - The scale of the outcome, based on the chosen family and link function
#' examples are "RR", "OR", "SMD", "HR"
#' @return \code{trt.map.table} - Treatments mapped to integer numbers, used to run BUGS code.
#' 
#' @details
#' 
#' @examples
#' #Example 1
#' #Fixed effects, consistency model.
#' #Binomial family, log link. This implies that the relative risk (RR) will be the scale.
#' 
#'nma.model(data.nma = my.slr,
#'        outcome = "n_died", 
#'        N = "n",
#'        reference = "plbo",
#'        family = "binomial",
#'        link = "log",
#'        effects = "fixed")
#'        
#' #Example 2
#' #Random Effects, inconsistency model
#' #Normal family, identity link.        
#'  
#'nma.model(data.nma = my.slr,
#'         outcome = "weight_change", 
#'         N = "n",
#'         reference = "control",
#'         type = "inconsistency",
#'         sd = "std_dev",
#'         family = "normal",
#'         link = "identity",
#'         effects = "random")
#'         
#'# Example 3
#'# Random Effects, consistency model,
#'# Poisson family, cloglog link.
#'
#'nma.model(data.nma = my.slr,
#'         outcome = "n_died", 
#'         N = "n",
#'         reference = "plbo",
#'         exposure.time="personyears",
#'         family = "poisson",
#'         link = "cloglog",
#'         effects = "random")

nma.model <- function(data.nma,
                     outcome, 
                     N,
                     reference,
                     type=NULL,
                     sd=NULL,
                     exposure.time=NULL,
                     family,
                     link,
                     effects,
                     prior.mu = "DEFAULT",
                     prior.d = "DEFAULT",
                     prior.sigma = "DEFAULT",
                     meta.covariate = NULL){
  if(link=="logit" & family %in% c("binomial", "binary", "bin", "binom")){
    scale <- "OR"
  }else if(link=="log" & family %in% c("binomial", "binary", "bin", "binom")){
    scale <- "RR"
  }else if(link== "identity" & family =="normal"){
    scale <- "SMD"
  }else if(link =="cloglog" & family %in% c("binomial", "binary", "bin", "binom")){
    scale <- "lograte"
  }else if(link == "log" & family =="poisson"){
    scale <- "HR"
  }
  
  data <- data.nma$raw.data
  
  # rename variables as appropriate
  names(data)[names(data) == data.nma$varname.t] <- "trt"
  names(data)[names(data) == data.nma$varname.s] <- "trial"
  names(data)[names(data) == outcome] <- "r1"
  names(data)[names(data) == N] <- "N"
  if (family == "normal"){names(data)[names(data) == sd] <- "sd"}
  if (!is.null(exposure.time)){names(data)[names(data) == exposure.time] <- "timevar"}
  if (!is.null(meta.covariate)){names(data)[names(data) == meta.covariate] <- "covariate"}
  
  #Trt mapping
  trt.map.table <- data$trt %>%
    unique %>% 
    sort  %>%
    tibble(trt.ini=.) %>%
    filter(trt.ini!=reference) %>%
    add_row(trt.ini=reference, .before=1) %>%
    mutate(trt.jags = 1:dim(.)[1])
  
  data %<>% mutate(trt.jags=mapvalues(trt,
                                      from=trt.map.table$trt.ini,
                                      to=trt.map.table$trt.jags) %>% as.integer)
  
  nt = data.nma$treatments %>% nrow()
  ns = data.nma$studies %>% nrow()
  na = data.nma$n.arms %>% select(n.arms) %>% t() %>% as.vector
  
  
  
  sorted.data <- data 
  
  # rename variables as appropriate
  names(sorted.data)[names(sorted.data) == data.nma$varname.t] <- "trt"
  names(sorted.data)[names(sorted.data) == data.nma$varname.s] <- "trial"
  names(sorted.data)[names(sorted.data) == outcome] <- "r1"
  names(sorted.data)[names(sorted.data) == N] <- "N"
  if (family == "normal"){names(sorted.data)[names(sorted.data) == sd] <- "sd"}
  if (!is.null(exposure.time)){names(sorted.data)[names(sorted.data) == exposure.time] <- "timevar"}
  if (!is.null(meta.covariate)){names(sorted.data)[names(sorted.data) == meta.covariate] <- "covariate"}
    
  sorted.data %<>% arrange(trial, trt.jags)
  
  if (family == "binomial" && link %in% c("log","logit")){
    r <- matrix(NA, ns, max(na))
    n <- matrix(NA, ns, max(na))
    t <- matrix(NA, ns, max(na))
    x <- matrix(NA, ns, max(na))
    
    line <- 1
    
    for(i in 1:ns){
      for(a in 1:na[i]){
        r[i,a] <- sorted.data %>% select(r1) %>% slice(line+a-1) %>% as.numeric()
        n[i,a] <- sorted.data %>% select(N) %>% slice(line+a-1) %>% as.numeric()
        t[i,a] <- sorted.data %>% select(trt.jags) %>% slice(line+a-1) %>% as.numeric()
        if (!is.null(meta.covariate)){x[i,a] <- sorted.data %>% select(covariate) %>% slice(line+a-1) %>% as.numeric()} #fix
      }
      line <- line + na[i]
    }
    
    if (!is.null(meta.covariate)) {x <- x-mean(x, na.rm=TRUE)}
    
    bugsdata2 <- list(ns=ns,
                      nt=nt,
                      #meanA=0,
                      #precA=0.5,
                      r=r,
                      n=n,
                      t=t,
                      na=na,
                      x=x) 
  } else if (family == "normal" && link=="identity"){
    y <- matrix(NA, ns, max(na))
    n <- matrix(NA, ns, max(na))
    sd<- matrix(NA, ns, max(na))
    se<- matrix(NA, ns, max(na))
    t <- matrix(NA, ns, max(na))
    x <- matrix(NA, ns, max(na))
    
    line <- 1
    
    for(i in 1:ns){
      for(a in 1:na[i]){
        y[i,a] <- sorted.data %>% select(r1) %>% slice(line+a-1) %>% as.numeric()
        sd[i,a]<- sorted.data %>% select(sd) %>% slice(line+a-1) %>% as.numeric()
        n[i,a] <- sorted.data %>% select(N) %>% slice(line+a-1) %>% as.numeric()
        t[i,a] <- sorted.data %>% select(trt.jags) %>% slice(line+a-1) %>% as.numeric()
        if (!is.null(meta.covariate)){x[i,a] <- sorted.data %>% select(covariate) %>% slice(line+a-1) %>% as.numeric()}        
      }
      line <- line + na[i]
    }
    se <- sd/sqrt(n)
    if (!is.null(meta.covariate)) {x <- x-mean(x, na.rm=TRUE)}
    
    bugsdata2 <- list(ns=ns,
                      nt=nt,
                      #meanA=0,
                      #precA=0.5,
                      y=y,
                      se=se,
                      t=t,
                      na=na,
                      x = x)
    
  } else if (family == "poisson" && link == "log"){
    E <- matrix(NA, ns, max(na))
    r <- matrix(NA, ns, max(na))
    t <- matrix(NA, ns, max(na))
    x <- matrix(NA, ns, max(na))
    
    line <- 1
    
    for(i in 1:ns){
      for(a in 1:na[i]){
        E[i,a] <- sorted.data %>% select(timevar) %>% slice(line+a-1) %>% as.numeric()
        r[i,a] <- sorted.data %>% select(r1) %>% slice(line+a-1) %>% as.numeric()
        t[i,a] <- sorted.data %>% select(trt.jags) %>% slice(line+a-1) %>% as.numeric()
        if (!is.null(meta.covariate)){x[i,a] <- sorted.data %>% select(covariate) %>% slice(line+a-1) %>% as.numeric()}
      }
      line <- line + na[i]
    }
    if (!is.null(meta.covariate)) {x <- x-mean(x, na.rm=TRUE)}
    
    bugsdata2 <- list(ns=ns,
                      nt=nt,
                      #meanA=0,
                      #precA=0.5,
                      E=E,
                      r=r,
                      t=t,
                      na=na,
                      x=x) 
  } else if (family == "binomial" && link == "cloglog"){
    time <- matrix(NA, ns, max(na))
    n <- matrix(NA, ns, max(na))
    r <- matrix(NA, ns, max(na))
    t <- matrix(NA, ns, max(na))
    x <- matrix(NA, ns, max(na))
    
    line <- 1
    
    for(i in 1:ns){
      for(a in 1:na[i]){
        time[i,a] <- sorted.data %>% select(timevar) %>% slice(line+a-1) %>% as.numeric()
        n[i,a] <- sorted.data %>% select(N) %>% slice(line+a-1) %>% as.numeric()
        r[i,a] <- sorted.data %>% select(r1) %>% slice(line+a-1) %>% as.numeric()
        t[i,a] <- sorted.data %>% select(trt.jags) %>% slice(line+a-1) %>% as.numeric()
        if (!is.null(meta.covariate)){x[i,a] <- sorted.data %>% select(covariate) %>% slice(line+a-1) %>% as.numeric()}        
      }
      line <- line + na[i]
    }
    if (!is.null(meta.covariate)) {x <- x-mean(x, na.rm=TRUE)}
    
    bugsdata2 <- list(ns=ns,
                      nt=nt,
                      #meanA=0,
                      #precA=0.5,
                      time=time,
                      n=n,
                      r=r,
                      t=t,
                      na=na,
                      x=x) 
  } 
  add.to.model <- trt.map.table %>% 
    transmute(Treatments = paste0("# ", trt.jags, ": ", trt.ini, "\n")) %>% 
    t() %>% 
    paste0() %>% 
    paste(collapse="") %>%
    paste0("\n\n# Treatment key\n",.)
  
  ############
  ###Priors###
  ############
  
  max.delta <- paste0(compute.prior(data.nma=data.nma, outcome=outcome, scale=scale, N=N, sd=sd, time = exposure.time))
  
  # BASELINE EFFECTS PRIOR
  if (prior.mu == "DEFAULT"){
    prior.mu.str <- sprintf("for(i in 1:ns){
                               mu[i] ~ dnorm(0,(%s*15)^(-2))
                               }", max.delta)
  } else {
    prior.mu.str <- sprintf("for(i in 1:ns){
                               mu[i] ~ %s
  }", prior.mu)
	}
  
  # RELATIVE EFFECTS PRIOR
  if (prior.d =="DEFAULT"){
    if(is.null(type)){
      prior.d.str <- sprintf("for(k in 2:nt){
                             d[k] ~ dnorm(0,(%s*15)^(-2))
    }", max.delta)

  } else if(type=="inconsistency"){
    prior.d.str <- sprintf("for (c in 1:(nt-1)) {
                           for (k in (c+1):nt)  { 
                           d[c,k] ~ dnorm(0,(%s*15)^(-2))
                           } 
    }", max.delta)
		}
  } else {
    if(is.null(type)){
      prior.d.str <- sprintf("for(k in 2:nt){
                             d[k] ~ %s)
    }", prior.d)

  }else if(type=="inconsistency"){
    prior.d.str <- sprintf("for (c in 1:(nt-1)) {
                           for (k in (c+1):nt)  { 
                           d[c,k] ~ %s
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
   if(scale == "RR"){
     prior.mu.str <-  "for(i in 1:ns){
     mu[i] <- log(p.base[i])           #priors for all trial baselines
     p.base[i] ~ dunif(0, 1)
   }"
   }

  #meta regression string
  if (!is.null(meta.covariate)){
    prior.meta.reg <- "beta[1]<-0
    for (k in 2:nt){
      beta[k] <- B
    }
    B~dnorm(0, .0001)"
  } else {
    prior.meta.reg <- ""
  }
  
  #remove covariate from bugsdata2 if unused
  if (is.null(meta.covariate)){bugsdata2 <- bugsdata2[names(bugsdata2)!="x"]}
  

  
  if(!is.null(type)){
    if(type=="inconsistency"){     
      model.str <- makeBUGScode(family=family,
                                link=link,
                                effects=effects,
                                inconsistency=TRUE,
                                prior.mu.str,
                                prior.d.str,
                                prior.sigma2.str,
                                meta.covariate,
                                prior.meta.reg) %>%
        paste0(add.to.model)
    }
  } else if(is.null(type)){
    model.str <- makeBUGScode(family=family, 
                              link=link, 
                              effects=effects, 
                              inconsistency=FALSE, 
                              prior.mu.str,
                              prior.d.str, 
                              prior.sigma2.str, 
                              meta.covariate, 
                              prior.meta.reg) %>%
      paste0(add.to.model)
  }
  
  return(bugs=list(model.str=model.str,
                   bugsdata2=bugsdata2, 
                   scale=scale, 
                   trt.map.table=trt.map.table, 
                   family=family, 
                   link=link))
}

