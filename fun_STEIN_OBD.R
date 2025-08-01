#--------------------------OBD function for STEIN design-----------------------#

# Dependency to perform isotonic regression
# require(Iso)

fun_STEIN_OBD = function(doseDT, pT, w1 = 0.33, w2 = 1.09) {
  
  # Obtain the isotonically transformed values for toxicity probs
  PT = (doseDT$x + 0.05)/(doseDT$n + 0.1)
  ptilde = pava_rcpp(PT, doseDT$n+0.1)+0.001*seq(1,nrow(doseDT))
  ptilde[which(doseDT$n == 0)] = Inf
  ptilde[which(doseDT$keep == 0)] = Inf
  
  # Obtain the unimodal-isotonically transformed values for efficacy probs
  qtilde = peestimate(doseDT$y, doseDT$n)
  qtilde[which(doseDT$keep == 0)] = -Inf
  qtilde[which(doseDT$n == 0)] = -Inf
    
  # Utility computation
  u <- NULL
  for(i in 1:nrow(doseDT))
  {
    u[i] <- qtilde[i] - w1*ptilde[i] - ifelse(ptilde[i] > pT, w2 * ptilde[i], 0)
  }

  if(max(u) == -Inf)
  {
    opt <- 99
  } else {
    opt = doseDT$id[which.max(u)]
  }
  return(opt)
  
}
#------------------------------------------------------------------------------#

#--------------------------Efficacy estimation function------------------------#

peestimate = function(yE,n) {
  ndose = length(yE)
  lik = rep(0,ndose)
  pe = (yE+0.05)/(n+0.1)
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