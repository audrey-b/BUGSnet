#add warning messages for incompatible link and family

makeBUGScode <- function(family, link, effects, inconsistency, prior.mu.str, prior.d.str, prior.sigma2.str, meta.covariate, prior.meta.reg, auto){
  
  # Set up family and monitor strings for arm-based reporting trials
  
  if (family=="binomial"){
    family.str <- "r[i,k] ~ dbin(p[i,k],n[i,k]) # binomial likelihood"
    monitor.str <- "rhat[i,k] <- p[i,k] * n[i,k] # expected value of the numerators
    dev_a[i,k] <- 2 * (r[i,k] * (log(r[i,k])-log(rhat[i,k])) #Deviance contribution
    + (n[i,k]-r[i,k]) * (log(n[i,k]-r[i,k]) - log(n[i,k]-rhat[i,k])))"
  }
  
  if (family=="normal"){
    family.str <- "y[i,k] ~ dnorm(theta_a[i,k],prec[i,k])"
    monitor.str <- "prec[i,k] <- pow(se[i,k],-2)
    dev_a[i,k] <- (y[i,k]-theta_a[i,k])*(y[i,k]-theta_a[i,k])*prec[i,k] #Deviance contribution"
  }
  
  if (family=="poisson"){
    family.str <- "r[i,k] ~ dpois(theta_a[i,k]) # Poisson likelihood"
    monitor.str <- "theta_a[i,k] <- lambda[i,k]*E[i,k] # failure rate * exposure
    dev_a[i,k] <- 2*((theta_a[i,k]-r[i,k]) + r[i,k]*log(r[i,k]/theta_a[i,k])) #Deviance contribution"
  }
  
  # Set up family and monitor strings for contrast-based reporting trials
    
    family.str.c <- "for (k in 1:(na_c[i]-1)) {
                     for (j in 1:(na_c[i]-1)) {
                       Sigma[i,j,k] <- se.diffs[i,k+1]^2*(equals(j,k)) + var.ref[i,1]*(1-equals(j,k))
                     }
    }
    Omega[i, 1:(na_c[i]-1), 1:(na_c[i]-1)] <- inverse(Sigma[i,1:(na_c[i]-1), 1:(na_c[i]-1)])
    y_c[i,2:na_c[i]] ~ dmnorm(theta_c[i,2:na_c[i]], Omega[i, 1:(na_c[i]-1), 1:(na_c[i]-1)])"
    monitor.str.c <- "for(k in 1:(na_c[i]-1)) {
    ydiff[i,k] <- y_c[i,(k+1)]-theta_c[i,(k+1)]
    }"
  
    
    # Deviance strings
    dev.str.c <- "dev_c[i] <- t(ydiff[i,1:(na_c[i]-1)])%*%Omega[i,1:(na_c[i]-1), 1:(na_c[i]-1)]%*%ydiff[i,1:(na_c[i]-1)]
    resdev_c[i] <- dev_c[i]"
    
    dev.str <- "resdev_a[i] <- sum(dev_a[i,1:na_a[i]])"
  
    # TODO figure out metaregression thing
  # if (!is.null(meta.covariate)) {
  #   metareg.str <- "+ (beta[t[i,k]]-beta[t[i,1]])*(x[i,k])"
  #   } else {metareg.str <- ""}
  
  if (effects == "fixed"){
    
    # Set up link for arm-based reporting trials
    if (family == "binomial" && link=="logit"){
      link.str <- "logit(p[i,k]) <- mu[i] + d[t_a[i,k]] - d[t_a[i,1]]"
    }  else if (family == "binomial" && link=="log"){
      link.str <- "log(p[i,k]) <- mu[i] + d[t_a[i,k]] - d[t_a[i,1]]"
    } else if (family == "normal" && link == "identity"){
      link.str <- "theta_a[i,k] <- mu[i] + d[t_a[i,k]] - d[t_a[i,1]] # model for linear predictor"
    } else if (family == "poisson" && link=="log"){
      link.str <- "log(lambda[i,k]) <- mu[i] + d[t_a[i,k]] - d[t_a[i,1]] # model for linear predictor"
    } else if (family== "binomial" && link=="cloglog"){
      link.str <- "cloglog(p[i,k]) <- log(time[i,k]) + mu[i] + d[t_a[i,k]] - d[t_a[i,1]] # model for linear predictor"
    } 
    
    # Set up link for contrast-based reporting trials
    link.str.c <- "theta_c[i,k] <- d[t_c[i,k]] - d[t_c[i,1]]"

    # TODO metareg
      # link.str <- paste0(link.str, metareg.str)
    
    # Fixed Effects Consistency Model
    if(!inconsistency){
      
      code.str <- sprintf("#This code is adapted from
    #Dias, S., Welton, N.J., Sutton, A.J. & Ades, A.E. NICE DSU Technical Support Document 2: 
    #A Generalised Linear Modelling Framework for Pairwise and Network Meta-Analysis of Randomised
    #Controlled Trials. 2011; last updated September 2016 (available from http:
    #//www.nicedsu.org.uk).
                          
    # fixed effects model for multi-arm trials
                          
    %s
    
      for(i in 1:ns_a){                      # LOOP THROUGH ARM-BASED STUDIES
      
        for (k in 1:na_a[i]) {             # LOOP THROUGH ARMS
          %s
          %s
          %s
        }
      %s
      }
      
      for(i in 1:ns_c) {                  # LOOP THROUGH CONTRAST-BASED STUDIES
      
        %s
        %s
      
        for(k in 1:na_c[i]) {
      
          %s
      
        }
        
        %s
      
      }
      
      totresdev <- sum(resdev_a[], resdev_c[])
      d[1]<-0
      %s
                          
      %s
      
      %s               
    %s", paste0(ifelse(auto, "", "model{                               # *** PROGRAM STARTS")),
        family.str,
        monitor.str,
        link.str,
        dev.str,
        family.str.c,
        monitor.str.c,
        link.str.c,
        dev.str.c,
        prior.mu.str,
        prior.d.str,
        prior.meta.reg,
        paste0(ifelse(auto, "", "}")))
    }
    
    # if(inconsistency){
    #   
    #   if (family == "binomial" && link=="logit"){
    #     link.str <- "logit(p[i,k]) <- mu[i] + d[t[i,1],t[i,k]]"
    #   }  else if (family == "binomial" && link=="log"){
    #     link.str <- "log(p[i,k]) <- mu[i] + d[t[i,1],t[i,k]]"
    #   } else if (family == "normal" && link == "identity"){
    #     link.str <- "theta[i,k] <- mu[i] + d[t[i,1],t[i,k]]"
    #   } else if (family == "poisson" && link=="log"){
    #     link.str <- "log(lambda[i,k]) <- mu[i] + d[t[i,1],t[i,k]]"
    #   } else if (family== "binomial" && link=="cloglog"){
    #     link.str <- "cloglog(p[i,k]) <- log(time[i,k]) + mu[i] + d[t[i,1],t[i,k]]"
    #   } else if (family == "contrast" && link == "identity") {
    #     link.str <- "theta[i,k] <- d[t[i,1],t[i,k]]"
    #   }
    #   
    #   link.str <- paste0(link.str, metareg.str)
    #   
    #   code.str <- sprintf("# inconsistency model
    # # fixed effects model
    # 
    # %s
    #                       
    #   for(i in 1:ns){             # LOOP THROUGH STUDIES
    #     %s
    #     %s
    #     for (k in 1:na[i])  {   # LOOP THROUGH ARMS
    #       %s
    #       %s
    #       %s
    #       }
    #       %s 
    #   }
    #     totresdev<-sum(resdev[])
    #     for (k in 1:nt){d[k,k]<-0}  #set effects of k vs k to zero
    #     %s
    # 
    #     %s  
    # 
    #   %s
    #   ", paste0(ifelse(auto, "", "model{                      # *** PROGRAM STARTS")),
    #      paste0(ifelse(family == "contrast", family.str, "")), # If contrast-based, the trial likelihood is multivariate
    #      paste0(ifelse(family == "contrast", monitor.str, "")),
    #      paste0(ifelse(family == "contrast", "", monitor.str)),
    #     link.str,
    #     paste0(ifelse(family == "contrast", "", family.str)),
    #     dev.str,
    #     prior.mu.str,
    #     prior.d.str,
    #     paste0(ifelse(auto, "", "}")))
    # } 
  } 
  
  # if (effects == "random"){
  #   
  #   if (family == "binomial" && link=="logit"){
  #     link.str <- "logit(p[i,k]) <- mu[i] + delta[i,k]"
  #   } else if (family == "binomial" && link=="log"){
  #     link.str <- "log(p[i,k]) <- mu[i] + delta[i,k]"
  #   } else if (family == "normal"){
  #     link.str <- "theta[i,k] <- mu[i] + delta[i,k]"
  #   } else if (family == "poisson" && link=="log"){
  #     link.str <- "log(lambda[i,k]) <- mu[i] + delta[i,k]"
  #   } else if (family == "binomial" && link=="cloglog"){
  #     link.str <- "cloglog(p[i,k]) <- log(time[i,k]) + mu[i] + delta[i,k]"
  #   } else if (family == "contrast" && link == "identity") {
  #     link.str <- "theta[i,k] <- delta[i,k]"
  #   }
  #   
  #   link.str <- paste0(link.str, metareg.str)
  #   
  #   if(!inconsistency){
  #     
  #     code.str <- sprintf("
  #    
  #     # Random effects model for multi-arm trials
  #   
  #     %s
  # 
  #     for(i in 1:ns){                      # LOOP THROUGH STUDIES
  #       %s
  #       %s
  #       w[i,1] <- 0    # adjustment for multi-arm trials is zero for control arm
  #       delta[i,1] <- 0             # treatment effect is zero for control arm
  #       for (k in 1:na[i]) {             # LOOP THROUGH ARMS
  #         %s
  #         # model for linear predictor
  #         %s
  #         %s
  #       }
  #       %s
  #       for (k in 2:na[i]) {             # LOOP THROUGH ARMS
  #         # trial-specific LOR distributions
  #         delta[i,k] ~ dnorm(md[i,k],taud[i,k])
  #         # mean of LOR distributions, with multi-arm trial correction
  #         md[i,k] <-  d[t[i,k]] - d[t[i,1]] + sw[i,k]
  #         # precision of LOR distributions (with multi-arm trial correction)
  #         taud[i,k] <- pow(sigma2,-1) *2*(k-1)/k
  #         # adjustment, multi-arm RCTs
  #         w[i,k] <- (delta[i,k] - d[t[i,k]] + d[t[i,1]])
  #         # cumulative adjustment for multi-arm trials
  #         sw[i,k] <- sum(w[i,1:(k-1)])/(k-1)
  #       }
  #     }   
  #     totresdev <- sum(resdev[])
  #     d[1]<-0       # treatment effect is zero for reference treatment
  #     %s
  #     %s
  #     %s
  #     tau <- pow(sigma,-2)
  #     %s
  #   %s", paste0(ifelse(auto, "", "model{                               # *** PROGRAM STARTS")),
  #       paste0(ifelse(family == "contrast", family.str, "")),
  #       paste0(ifelse(family == "contrast", monitor.str, "")),
  #       paste0(ifelse(family == "contrast", "", monitor.str)),
  #       link.str,
  #       paste0(ifelse(family == "contrast", "", family.str)),
  #       dev.str,
  #       prior.mu.str,
  #       prior.d.str,
  #       prior.sigma2.str,
  #       prior.meta.reg,
  #       paste0(ifelse(auto, "", "}")))
  #     
  #   }
  #   
  #   if(inconsistency){
  #     
  #     code.str <- sprintf("# Binomial likelihood, inconsistency model
  #     # Random effects model
  #     %s
  #   
  #     for(i in 1:ns){             # LOOP THROUGH STUDIES
  #       delta[i,1]<-0           # treatment effect is zero in control arm
  #       %s
  #       %s
  #       for (k in 1:na[i])  {   # LOOP THROUGH ARMS
  #         %s
  #         %s
  #         %s
  #       }
  #       %s
  #       for (k in 2:na[i]) {  # LOOP THROUGH ARMS
  #         delta[i,k] ~ dnorm(d[t[i,1], t[i,k]] , pow(sigma2,-1)) # trial-specific LOR distributions
  #       }
  #     }
  #     totresdev <- sum(resdev[])
  #     %s
  #     %s  
  #     %s
  #     %s
  #     tau <- 1/sigma2
  # 
  #   %s
  #   ", paste0(ifelse(auto, "", "model{                      # *** PROGRAM STARTS")),
  #                         paste0(ifelse(family == "contrast", family.str, "")),
  #                         paste0(ifelse(family == "contrast", monitor.str, "")),
  #                         paste0(ifelse(family == "contrast", "", monitor.str)),
  #                         link.str,
  #                         paste0(ifelse(family == "contrast", "", family.str)),
  #                         dev.str,
  #                         prior.mu.str,
  #                         prior.d.str,
  #                         prior.sigma2.str,
  #                         prior.meta.reg,
  #                         paste0(ifelse(auto, "", "}")))
  #   }
  # }
return(code.str)

}