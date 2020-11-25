#add warning messages for incompatible link and family

makeBUGScode <- function(family, link, effects, inconsistency, prior.mu.str, prior.d.str, prior.sigma2.str, meta.covariate, prior.meta.reg){
  
  if (family=="binomial"){
    family.str <- "r[i,k] ~ dbin(p[i,k],n[i,k]) # binomial likelihood"
    monitor.str <- "rhat[i,k] <- p[i,k] * n[i,k] # expected value of the numerators
    dev[i,k] <- 2 * (r[i,k] * (log(r[i,k])-log(rhat[i,k])) #Deviance contribution
    + (n[i,k]-r[i,k]) * (log(n[i,k]-r[i,k]) - log(n[i,k]-rhat[i,k])))"
  }
  
  if (family=="normal"){
    family.str <- "y[i,k] ~ dnorm(theta[i,k],prec[i,k])"
    monitor.str <- "prec[i,k] <- pow(se[i,k],-2)
    dev[i,k] <- (y[i,k]-theta[i,k])*(y[i,k]-theta[i,k])*prec[i,k] #Deviance contribution"
  }
  
  if (family=="poisson"){
    family.str <- "r[i,k] ~ dpois(theta[i,k]) # Poisson likelihood"
    monitor.str <- "theta[i,k] <- lambda[i,k]*E[i,k] # failure rate * exposure
    dev[i,k] <- 2*((theta[i,k]-r[i,k]) + r[i,k]*log(r[i,k]/theta[i,k])) #Deviance contribution"
  }
  
  if (!is.null(meta.covariate)) {
    metareg.str <- "+ (beta[t[i,k]]-beta[t[i,1]])*(x[i,k])"
    } else {metareg.str <- ""}
  
  
  if (effects == "fixed"){
    
    if (family == "binomial" && link=="logit"){
      link.str <- "logit(p[i,k]) <- mu[i] + d[t[i,k]] - d[t[i,1]]"
    }  else if (family == "binomial" && link=="log"){
      link.str <- "log(p[i,k]) <- mu[i] + d[t[i,k]] - d[t[i,1]]"
    } else if (family == "normal" && link == "identity"){
      link.str <- "theta[i,k] <- mu[i] + d[t[i,k]] - d[t[i,1]] # model for linear predictor"
    } else if (family == "poisson" && link=="log"){
      link.str <- "log(lambda[i,k]) <- mu[i] + d[t[i,k]] - d[t[i,1]] # model for linear predictor"
    } else if (family== "binomial" && link=="cloglog"){
      link.str <- "cloglog(p[i,k]) <- log(time[i,k]) + mu[i] + d[t[i,k]] - d[t[i,1]] # model for linear predictor"
    }

    link.str <- paste0(link.str, metareg.str)
    
    if(!inconsistency){
      
      code.str <- sprintf("#This code is adapted from
    #Dias, S., Welton, N.J., Sutton, A.J. & Ades, A.E. NICE DSU Technical Support Document 2: 
    #A Generalised Linear Modelling Framework for Pairwise and Network Meta-Analysis of Randomised
    #Controlled Trials. 2011; last updated September 2016 (available from http:
    #//www.nicedsu.org.uk).
                          
    # fixed effects model for multi-arm trials
                          
    model{                               # *** PROGRAM STARTS

      for(i in 1:ns){                      # LOOP THROUGH STUDIES
      for (k in 1:na[i]) {             # LOOP THROUGH ARMS
        %s
        %s
        %s
      }
      resdev[i] <- sum(dev[i,1:na[i]])

    }   
      totresdev <- sum(resdev[])
      d[1]<-0
      %s
                          
      %s
      
      %s               
    }", monitor.str,
        link.str,
        family.str,
        prior.mu.str,
        prior.d.str,
        prior.meta.reg)
    }
    
    if(inconsistency){
      
      if (family == "binomial" && link=="logit"){
        link.str <- "logit(p[i,k]) <- mu[i] + d[t[i,1],t[i,k]]"
      }  else if (family == "binomial" && link=="log"){
        link.str <- "log(p[i,k]) <- mu[i] + d[t[i,1],t[i,k]]"
      } else if (family == "normal" && link == "identity"){
        link.str <- "theta[i,k] <- mu[i] + d[t[i,1],t[i,k]]"
      } else if (family == "poisson" && link=="log"){
        link.str <- "log(lambda[i,k]) <- mu[i] + d[t[i,1],t[i,k]]"
      } else if (family== "binomial" && link=="cloglog"){
        link.str <- "cloglog(p[i,k]) <- log(time[i,k]) + mu[i] + d[t[i,1],t[i,k]]"
      }
      
      link.str <- paste0(link.str, metareg.str)
      
      code.str <- sprintf("# inconsistency model
    # fixed effects model
    model{                      # *** PROGRAM STARTS
                          
      for(i in 1:ns){             # LOOP THROUGH STUDIES
        for (k in 1:na[i])  {   # LOOP THROUGH ARMS
          %s
          %s
          %s
          }
          resdev[i]<-sum(dev[i,1:na[i]]) 
      }
        totresdev<-sum(resdev[])
        for (k in 1:nt){d[k,k]<-0}  #set effects of k vs k to zero
        %s

        %s  
        
        %s

      }
      ", monitor.str,
        link.str,
        family.str,
        prior.mu.str,
        prior.d.str,
        prior.meta.reg)
    } 
  } 
  
  if (effects == "random"){
    
    if (family == "binomial" && link=="logit"){
      link.str <- "logit(p[i,k]) <- mu[i] + delta[i,k]"
    } else if (family == "binomial" && link=="log"){
      link.str <- "log(p[i,k]) <- mu[i] + delta[i,k]"
    } else if (family == "normal"){
      link.str <- "theta[i,k] <- mu[i] + delta[i,k]"
    } else if (family == "poisson" && link=="log"){
      link.str <- "log(lambda[i,k]) <- mu[i] + delta[i,k]"
    } else if (family == "binomial" && link=="cloglog"){
      link.str <- "cloglog(p[i,k]) <- log(time[i,k]) + mu[i] + delta[i,k]"
    }
    
    link.str <- paste0(link.str, metareg.str)
    
    if(!inconsistency){
      
      code.str <- sprintf("
     
      # Random effects model for multi-arm trials
    
      model{                               # *** PROGRAM STARTS
  
      for(i in 1:ns){                      # LOOP THROUGH STUDIES
        w[i,1] <- 0    # adjustment for multi-arm trials is zero for control arm
        delta[i,1] <- 0             # treatment effect is zero for control arm
        for (k in 1:na[i]) {             # LOOP THROUGH ARMS
          %s
          # model for linear predictor
          %s
          %s
        }
        resdev[i] <- sum(dev[i,1:na[i]]) #JS
        for (k in 2:na[i]) {             # LOOP THROUGH ARMS
          # trial-specific LOR distributions
          delta[i,k] ~ dnorm(md[i,k],taud[i,k])
          # mean of LOR distributions, with multi-arm trial correction
          md[i,k] <-  d[t[i,k]] - d[t[i,1]] + sw[i,k]
          # precision of LOR distributions (with multi-arm trial correction)
          taud[i,k] <- pow(sigma2,-1) *2*(k-1)/k
          # adjustment, multi-arm RCTs
          w[i,k] <- (delta[i,k] - d[t[i,k]] + d[t[i,1]])
          # cumulative adjustment for multi-arm trials
          sw[i,k] <- sum(w[i,1:(k-1)])/(k-1)
        }
      }   
      totresdev <- sum(resdev[])
      d[1]<-0       # treatment effect is zero for reference treatment
      %s
      %s
      %s
      tau <- pow(sigma,-2)
      %s
    }", monitor.str,
        link.str,
        family.str,
        prior.mu.str,
        prior.d.str,
        prior.sigma2.str,
        prior.meta.reg)
      
    }
    
    if(inconsistency){
      
      code.str <- sprintf("# Binomial likelihood, inconsistency model
      # Random effects model
      model{                      # *** PROGRAM STARTS
    
      for(i in 1:ns){             # LOOP THROUGH STUDIES
        delta[i,1]<-0           # treatment effect is zero in control arm
        for (k in 1:na[i])  {   # LOOP THROUGH ARMS
          %s
          %s
          %s
        }
        resdev[i] <- sum(dev[i,1:na[i]])
        for (k in 2:na[i]) {  # LOOP THROUGH ARMS
          delta[i,k] ~ dnorm(d[t[i,1], t[i,k]] , pow(sigma2,-1)) # trial-specific LOR distributions
        }
      }
      totresdev <- sum(resdev[])
      %s
      %s  
      %s
      tau <- 1/sigma2
      %s

    }
    ",monitor.str,
      link.str,
      family.str,
      prior.mu.str,
      prior.d.str,
      prior.sigma2.str,
      prior.meta.reg)
    }
  }
return(code.str)

}