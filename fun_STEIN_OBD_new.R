#--------------------------OBD function for STEIN design-----------------------#

# Dependency to perform isotonic regression
# require(Iso)

fun_STEIN_OBD_new = function(doseDT, pT, psi1, w1 = 0.33, w2 = 1.09, size = 1000, pgrad = 0.1, 
                             alphaTO = 0.5, betaTO = 0.5, alphaEO = 0.5, betaEO = 0.5) {
  
  DTremain <- doseDT[doseDT$keep == 1 & doseDT$n > 0, ]  # need to select doses applied to patients
  rN <- nrow(DTremain)
  
  u_brench <- psi1 - w1 * pT 
    
  if(rN == 0)
  {
    return(99)
  } else
  {
    alphaTpost <- alphaTO + DTremain$x
    betaTpost <- betaTO + DTremain$n - DTremain$x
    alphaEpost <- alphaEO + DTremain$y
    betaEpost <- betaEO + DTremain$n - DTremain$y
    
    meanTpost <- alphaTpost/(alphaTpost + betaTpost)
    varTpost <- meanTpost * (1-meanTpost)/(1+alphaTpost+betaTpost)
    
    if(rN > 1)
    {
      # Obtain the isotonically transformed values for toxicity probs
      PT = (doseDT$x + 0.05)/(doseDT$n + 0.1)
      phat = pava_rcpp(PT, doseDT$n+0.1)+0.001*seq(1,nrow(doseDT))
      phat[which(doseDT$n == 0)] = Inf
      phat[which(doseDT$keep == 0)] = Inf
      
      # Obtain the unimodal-isotonically transformed values for efficacy probs
      qhat = peestimate_new(doseDT$y, doseDT$n, doseDT$y[doseDT$n > 0]/doseDT$n[doseDT$n > 0])
      qhat[which(doseDT$keep == 0)] = -Inf
      qhat[which(doseDT$n == 0)] = -Inf
      
      u_hat <- NULL
      for(i in 1:nrow(doseDT))
      {
        u_hat[i] <- qhat[i] - w1*phat[i] - ifelse(phat[i] > pT, w2 * phat[i], 0)
      }
      opt <- doseDT$id[which.max(u_hat)]
      
      D <- length(alphaTpost)
      ptilde <- sapply(1:D, FUN = function(x) rbeta(size, alphaTpost[x], betaTpost[x])) %>%
        apply(., MARGIN = 1, FUN = pava_rcpp, w = 1/varTpost) %>%
        t() # size * D
      
      # Obtain the unimodal-isotonically transformed values for efficacy probs
      sampleE <- sapply(1:D, FUN = function(x) rbeta(size, alphaEpost[x], betaEpost[x])) # size * D
      qtilde <- sapply(1:size, FUN = function(x) peestimate_new(doseDT$y, doseDT$n, sampleE[x, ])) %>% t()
      qtilde <- qtilde[,DTremain$id]
      
      u <- qtilde - w1 * ptilde - w2 * ptilde * (ptilde > pT)
      ubar <- colMeans(u)
      u_best <- u[,which.max(ubar)]

      # opt <- DTremain$id[which.max(ubar)]
      pin <- mean((ptilde[,which.max(ubar)] <= pT) * (qtilde[,which.max(ubar)] >= psi1), na.rm = T)
    } else   # rN = 1
    {
      opt <- DTremain$id
      ptilde <- rbeta(size, alphaTpost, betaTpost)
      qtilde <- rbeta(size, alphaEpost, betaEpost)
      
      u_best <- qtilde - w1 * ptilde - w2 * ptilde * (ptilde > pT)
    }
    
    if(mean(u_best > u_brench) < pgrad)
    {
      return(99)
    } else
    {
      return(opt)
    }
    
  }
  
  
  
}
#------------------------------------------------------------------------------#

#--------------------------Efficacy estimation function------------------------#

peestimate_new <- function(yE,n, qhat) {
  ndose = length(yE)
  lik = rep(0,ndose)
  pe = (yE+0.005)/(n+0.01)
  pe[n > 0] <- qhat
  p.e = matrix(NA,ndose,ndose)
  
  for (i in 1:ndose) {
    
    if (i==1) {x=seq(ndose,1,by=-1)} else {x = c(1:(i-1),seq(ndose,i))}
    
    p.e[i,] = ufit(pe,lmode=i,x=1:ndose,w=n+0.5)[[2]]
    lik[i] = prod(dbinom(yE,n,p.e[i,]))		
  }
  lik = lik/sum(lik)
  pe = t(p.e)%*%lik+0.01*seq(1,ndose)
  
  return(pe)
}
#------------------------------------------------------------------------------#