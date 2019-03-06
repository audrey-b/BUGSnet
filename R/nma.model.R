#' Create Bugs Model
#' @description Creates BUGS code which can be ran through \code{nma.run()}.
#' 
#' @param data A data object produced by \code{data.prep()}
#' @param outcome A string indicating the name of your outcome variable.
#' @param N A string indicating the name of the variable containing the number of participants in each arm
#' @param sd A string (only required for continuous outcomes) indicating variable name
#' of the standard deviation of the outcome. Standard errors should be converted to standard deviation by multiplying by the sample size prior to using this function.
#' @param reference A string for the treatment that will be seen as the 'referent' comparator and labeled as treatment 1 in the BUGS code. This is often
#' a placebo or control drug of some kind.  
#' @param family A string indicating the family of the distribution of the outcome. Options are:
#' "binomial", "normal", "poisson".
#' @param link The link function for the nma model. Options are "logit", "log", "cloglog", "identity".
#' @param time A string (only required for binomial-cloglog or poisson-log models) indicating the name of variable 
#'   indicating person-time followup (e.g person years) or study followup time.
#' @param effects A string indicating the type of treatment effect relative to baseline. Options are "fixed" or "random".
#' @param prior.mu A string of BUGS code that defines priors on the baseline treatment effects. By default, independent normal priors are used with mean 0 and standard deviation 15u, where u is the largest maximum likelihood estimator in single trials \insertCite{@see @gemtc}{BUGSnet}.
#' @param prior.d A string of BUGS code that defines define priors on relative treatment effects. By default, independent normal priors are used with mean 0 and standard deviation 15u, where u is the largest maximum likelihood estimator in single trials \insertCite{@see @gemtc}{BUGSnet}.
#' @param prior.sigma A string of BUGS code that defines the prior on the variance of relative treatment effects. By default, a uniform distribution with range 0 to u is used, where u is the largest maximum likelihood estimator in single trials \insertCite{@see @gemtc}{BUGSnet}.
#' @param prior.beta Optional string that defines the prior on the meta-regression coefficients. Options are "UNRELATED", "EXCHANGEABLE", "EQUAL" \insertCite{@TSD3}{BUGSnet} or a string of BUGS code.
#' @param covariate Optional string indicating the name of the variable in your data set that you would like to
#' adjust for via meta regression. By default, covariate=NULL and no covariate adjustment is applied. The covariate will be centered for the analysis.
#' @param type If type="inconsistency", an inconsistency model will be built. By default, type="consistency" and a consistency model is built.
#' will be built.
#' 
#' @return \code{model} - A long character string containing BUGS code that will be run in \code{jags}.
#' @return \code{data} - The data used in the BUGS code.
#' @return \code{scale} - The scale of the outcome, based on the chosen family and link function
#' examples are "Risk Ratio" (relative risk), "Odds Ratio", "Mean Difference", "Hazard Ratio"
#' @return \code{trt.key} - Treatments mapped to integer numbers, used to run BUGS code.
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
#' #Example 1
#' #Random effects, consistency model.
#' #Binomial family, cloglog link. This implies that the scale will be the Hazard Ratio.
#' 
#' data(diabetes.sim)
#' diabetes.slr <- data.prep(arm.data = diabetes.sim, 
#' varname.t = "Treatment", 
#' varname.s = "Study")
#' 
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
#' #Example 2
#' #Fixed effects, consistency model.
#' #Binomial family, cloglog link. This implies that the scale will be the Hazard Ratio.
#' 
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


nma.model <- function(data,
                      outcome, 
                      N,
                      sd=NULL,
                      reference,
                      type="consistency",
                      time=NULL,
                      family,
                      link,
                      effects,
                      prior.mu = "DEFAULT",
                      prior.d = "DEFAULT",
                      prior.sigma = "DEFAULT",
                      prior.beta = NULL,
                      covariate = NULL){
  
  if(!is.null(covariate) & is.null(prior.beta))stop("prior.beta must be specified when covariate is specified")
  if(is.null(covariate) & !is.null(prior.beta))stop("covariate must be specified when prior.beta is specified")
  if(!is.null(prior.beta)){
    if(!(prior.beta %in% c("UNRELATED","EQUAL","EXCHANGEABLE"))){
      stop("prior.beta must be either UNRELATED, EQUAL, or EXCHANGEABLE")
    }
  }
  if(family=="normal" & is.null(sd)) stop("sd must be specified for continuous outcomes")
  if(family=="normal" & link!="identity") stop("This combination of family and link is currently not supported in BUGSnet.")
  if(family=="poisson" & link!="log") stop("This combination of family and link is currently not supported in BUGSnet.")
  if(family=="binomial" & !(link %in% c("log","logit", "cloglog"))) stop("This combination of family and link is currently not supported in BUGSnet.")
  
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
  
  data1 <- data$arm.data
  
  # rename variables as appropriate
  names(data1)[names(data1) == data$varname.t] <- "trt"
  names(data1)[names(data1) == data$varname.s] <- "trial"
  names(data1)[names(data1) == outcome] <- "r1"
  names(data1)[names(data1) == N] <- "N"
  if (family == "normal"){names(data1)[names(data1) == sd] <- "sd"}
  if (!is.null(time)){names(data1)[names(data1) == time] <- "timevar"}
  if (!is.null(covariate)){names(data1)[names(data1) == covariate] <- "covariate"}
  
  #Trt mapping
  trt.key <- data1$trt %>%
    unique %>% 
    sort  %>%
    tibble(trt.ini=.) %>%
    filter(trt.ini!=reference) %>%
    add_row(trt.ini=reference, .before=1) %>%
    mutate(trt.jags = 1:dim(.)[1])
  
  data1 %<>% mutate(trt.jags=mapvalues(trt,
                                       from=trt.key$trt.ini,
                                       to=trt.key$trt.jags) %>% as.integer)
  
  nt = data$treatments %>% nrow()
  ns = data$studies %>% nrow()
  na = data$n.arms %>% select(n.arms) %>% t() %>% as.vector
  
  
  
  sorted.data <- data1 
  
  # rename variables as appropriate
  names(sorted.data)[names(sorted.data) == data$varname.t] <- "trt"
  names(sorted.data)[names(sorted.data) == data$varname.s] <- "trial"
  names(sorted.data)[names(sorted.data) == outcome] <- "r1"
  names(sorted.data)[names(sorted.data) == N] <- "N"
  if (family == "normal"){names(sorted.data)[names(sorted.data) == sd] <- "sd"}
  if (!is.null(time)){names(sorted.data)[names(sorted.data) == time] <- "timevar"}
  if (!is.null(covariate)){names(sorted.data)[names(sorted.data) == covariate] <- "covariate"}
  
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
        if (!is.null(covariate)){x[i,a] <- sorted.data %>% select(covariate) %>% slice(line+a-1) %>% as.numeric()} #fix
      }
      line <- line + na[i]
    }
    
    if (!is.null(covariate)) {
      mean.cov <- mean(x[,1], na.rm=TRUE)
      x <- x-mean.cov
    } else{mean.cov <- NULL}
    
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
        if (!is.null(covariate)){x[i,a] <- sorted.data %>% select(covariate) %>% slice(line+a-1) %>% as.numeric()}        
      }
      line <- line + na[i]
    }
    se <- sd/sqrt(n)
    if (!is.null(covariate)) {
      mean.cov <- mean(x, na.rm=TRUE)
      x <- x-mean.cov
    } else{mean.cov <- NULL}
    
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
        if (!is.null(covariate)){x[i,a] <- sorted.data %>% select(covariate) %>% slice(line+a-1) %>% as.numeric()}
      }
      line <- line + na[i]
    }
    if (!is.null(covariate)) {
      mean.cov <- mean(x, na.rm=TRUE)
      x <- x-mean.cov
    } else{mean.cov <- NULL}
    
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
    time.mat <- matrix(NA, ns, max(na))
    n <- matrix(NA, ns, max(na))
    r <- matrix(NA, ns, max(na))
    t <- matrix(NA, ns, max(na))
    x <- matrix(NA, ns, max(na))
    
    line <- 1
    
    for(i in 1:ns){
      for(a in 1:na[i]){
        time.mat[i,a] <- sorted.data %>% select(timevar) %>% slice(line+a-1) %>% as.numeric()
        n[i,a] <- sorted.data %>% select(N) %>% slice(line+a-1) %>% as.numeric()
        r[i,a] <- sorted.data %>% select(r1) %>% slice(line+a-1) %>% as.numeric()
        t[i,a] <- sorted.data %>% select(trt.jags) %>% slice(line+a-1) %>% as.numeric()
        if (!is.null(covariate)){x[i,a] <- sorted.data %>% select(covariate) %>% slice(line+a-1) %>% as.numeric()}        
      }
      line <- line + na[i]
    }
    if (!is.null(covariate)) {
      mean.cov <- mean(x, na.rm=TRUE)
      x <- x-mean.cov
    } else{mean.cov <- NULL}
    
    bugsdata2 <- list(ns=ns,
                      nt=nt,
                      #meanA=0,
                      #precA=0.5,
                      "time"=time.mat,
                      n=n,
                      r=r,
                      t=t,
                      na=na,
                      x=x) 
  } 
  add.to.model <- trt.key %>% 
    transmute(Treatments = paste0("# ", trt.jags, ": ", trt.ini, "\n")) %>% 
    t() %>% 
    paste0() %>% 
    paste(collapse="") %>%
    paste0("\n\n# Treatment key\n",.)
  
  ############
  ###Priors###
  ############
  
  max.delta <- paste0(nma.prior(data, outcome=outcome, scale=scale, N=N, sd=sd, time = time))
  
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
    if(type=="consistency"){
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
    if(type=="consistency"){
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
  if(scale == "Risk Ratio"){
    prior.mu.str <-  "for(i in 1:ns){
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
  #remove covariate from bugsdata2 if unused
  if (is.null(covariate)){bugsdata2 <- bugsdata2[names(bugsdata2)!="x"]}
  
  
  
  
  if(type=="inconsistency"){     
    model <- makeBUGScode(family=family,
                          link=link,
                          effects=effects,
                          inconsistency=TRUE,
                          prior.mu.str,
                          prior.d.str,
                          prior.sigma2.str,
                          covariate,
                          prior.meta.reg) %>%
      paste0(add.to.model)
  } else if(type=="consistency"){
    model <- makeBUGScode(family=family, 
                          link=link, 
                          effects=effects, 
                          inconsistency=FALSE, 
                          prior.mu.str,
                          prior.d.str, 
                          prior.sigma2.str, 
                          covariate, 
                          prior.meta.reg) %>%
      paste0(add.to.model)
  }
  
  return(model=list(bugs=model,
                    data=bugsdata2, 
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
                    mean.cov=mean.cov))
  
  
}

