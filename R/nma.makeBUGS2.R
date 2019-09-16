#add warning messages for incompatible link and family

makeBUGScode2 <- function(max.arm, family, link, effects, inconsistency, prior.mu.str, prior.d.str, prior.sigma2.str, meta.covariate, prior.meta.reg){
  
  if(max.arm == 2 & effects == "fixed"){
 
  if (!is.null(meta.covariate)) {
    metareg.str <- "+ (beta[t[i,2]]-beta[t[i,1]])*(x[i,2])"
    } else {metareg.str <- ""}
    
    link.str <- "theta[i,2] <- d[t[i,2]] - d[t[i,1]]"
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
      
         y[i,2] ~ dnorm(theta[i,2],prec[i,2]) # normal likelihood for 2-arm trials
         dev[i] <- (y[i,2]-theta[i,2])*(y[i,2]-theta[i,2])*prec[i,2]
         prec[i,2]<- pow(se[i,2],-2)
         %s
      }

      totresdev <- sum(dev[])
      
      d[1]<-0
      %s
      %s
      
    }", link.str,
        prior.d.str,
        prior.meta.reg)
    }
  }
 
  if(max.arm == 3 & effects == "fixed"){
    
    if (!is.null(meta.covariate)) {
      metareg.str <- "+ (beta[t[i,k]]-beta[t[i,1]])*(x[i,k])"
    } else {metareg.str <- ""}
    
    link.str <- "theta[i,k] <- d[t[i,k]] - d[t[i,1]]"
    link.str <- paste0(link.str, metareg.str)
    
    if(!inconsistency){
      
      code.str <- sprintf("#This code is adapted from
                          #Dias, S., Welton, N.J., Sutton, A.J. & Ades, A.E. NICE DSU Technical Support Document 2: 
                          #A Generalised Linear Modelling Framework for Pairwise and Network Meta-Analysis of Randomised
                          #Controlled Trials. 2011; last updated September 2016 (available from http:
                          #//www.nicedsu.org.uk).
                          
                          # fixed effects model for multi-arm trials
                          
                          model{                               # *** PROGRAM STARTS
                          
                          for(i in 1:ns2){                      # LOOP THROUGH STUDIES
                          
                          y[i,2] ~ dnorm(theta[i,2],prec[i,2]) # normal likelihood for 2-arm trials
                          dev[i] <- (y[i,2]-theta[i,2])*(y[i,2]-theta[i,2])*prec[i,2]
                         
                          }

                          for(i in (ns2+1):(ns2+ns3)){
                              for(k in 1:(na[i]-1)) {
                                  for(j in 1:(na[i]-1)){
                                  Sigma[i,j,k] <- V[i,1]*(1-equals(j,k)) + Variance[i,k+1]*equals(j,k)
                              }
                          }
                          # precision matrix
                          Omega[i,1:(na[i]-1),1:(na[i]-1)] <- inverse(Sigma[i,,])
                          
                          # multivariate normal likelihood for 3-arm studies  
                          y[i,2:na[i]]~dmnorm(theta[i,2:na[i]],Omega[i,1:(na[i]-1),1:(na[i]-1)])
                          
                          # deviance contribution
                          for(k in 1:(na[i]-1)){
                                ydiff[i,k] <- y[i,(k+1)] - theta[i, (k+1)]
                                z[i,k] <- inprod(Omega[i,k,1:(na[i]-1)], ydiff[i,1:(na[i]-1)])
                          }
                          dev[i] <- inprod(ydiff[i,1:(na[i]-1)], z[i,1:(na[i]-1)])
                          }


                          for(i in 1:(ns2+ns3)) {
                              for(k in 2:na[i]){
                                 prec[i,k]<- pow(se[i,k],-2)
                                 Variance[i,k] <- pow(se[i,k],2)
                                 %s
                              }
                          }
                          totresdev <- sum(dev[])
                          d[1]<-0

                          %s
                          %s
                          
                          }", link.str,
        prior.d.str,
        prior.meta.reg)
    }
  }
 
  
  
return(code.str)
}