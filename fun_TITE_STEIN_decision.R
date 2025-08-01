# dT: D: the maximum value of n_T0_tilde to get "D";
#     DU: the maximum value of n_T0_tilde to get "DU";
#     (1) [0, DU] -> DU; (2) (DU, D] -> D; (3) (D, n - n_T1) -> ND;
# dE: TBD: the minimum value of n_E0_tilde to get "TBD";
#     EU: the mimimum value of n_E0_tilde to get "EU";
#     (1) [0, TBD) -> S; (2) [TBD, EU) -> TBD; (3) [EU, n - n_E1] -> EU

fun_TITE_STEIN_decision <- function(n, pT, phi1, phi2, psi1, psi2, cutoff_tox, cutoff_eff, 
                                    alphaT = 1, betaT = 1, alphaE = 1, betaE = 1, onlyDM = TRUE)
{
    # Calculate phi_L and phi_U based on pre-specified phi1, phi2 and pT (phi0)
    phi_L = log((1-phi1)/(1-pT))/log((pT*(1-phi1))/(phi1*(1-pT))) # lower bound
    phi_U = log((1-pT)/(1-phi2))/log((phi2*(1-pT))/(pT*(1-phi2))) # upper bound
    
    # Calculate psi based on pre-specified psi1 and psi2
    psi = log((1-psi1)/(1-psi2))/log((psi2*(1-psi1))/(psi1*(1-psi2)))
    
    dT <- data.frame(n_T1 = 0:n, D = rep(NA, n + 1), DU = rep(-Inf, n+1))  # D: maximum value to get "D"
    dE <- data.frame(n_E1 = 0:n, TBD = rep(NA, n + 1), EU = rep(Inf, n+1))   # TBD: minimum value to get "TBD"

    for(n_T1 in 0:n)
    {
      n_T0_tilde <- seq(from = ifelse(n_T1 == 0, 0.01, 0), to = n - n_T1, by = 0.01)    # n_T0 + m_T: effective sample size without DLT

      p_hat <- n_T1/(n_T1 + n_T0_tilde)
      if(max(p_hat) < phi_U)               # no possible for D
      { 
        dT$D[n_T1 + 1] <- -Inf
      } else if(min(p_hat) >= phi_U)       # must be D 
      {
        dT$D[n_T1 + 1] <- Inf
      } else
      {
        dT$D[n_T1 + 1] <- n_T0_tilde[sum(p_hat >= phi_U)]
      }
      
      DU_list <- (1-pbeta(pT, n_T1+alphaT, n_T0_tilde + betaT)) > cutoff_tox
      if(sum(DU_list) > 0)
      {
        dT$DU[n_T1 + 1] <-  n_T0_tilde[sum(DU_list)]
      }
    }
    
    for(n_E1 in 0:n)
    {
      n_E0_tilde <- seq(from = ifelse(n_E1 == 0, 0.01, 0), to = n - n_E1, by = 0.01)    # n_E0 + m_T: effective sample size without efficacy response

      q_hat <- n_E1/(n_E1 + n_E0_tilde)
      if(max(q_hat) < psi)                      # no possible for S
      {
        dE$TBD[n_E1 + 1] <- -Inf
      } else if(min(q_hat) >= psi)              # no possbile for TBD 
      {
        dE$TBD[n_E1 + 1] <- Inf
      } else
      {
        dE$TBD[n_E1 + 1] <- n_E0_tilde[sum(q_hat >= psi)]
      }
      
      # EU_list <- pbeta(psi1, n_E1+alphaE, n_E0_tilde + betaE) > cutoff_eff
      EU_list <- pbeta(0.25, n_E1+alphaE, n_E0_tilde + betaE) > cutoff_eff
      if(sum(EU_list) > 0)
      {
        dE$EU[n_E1 + 1] <-  n_E0_tilde[length(n_E0_tilde) - sum(EU_list) + 1]
      }
    }
    
    decisionM <- list(dT = dT, dE = dE)
    
    if(onlyDM == TRUE)
    {
      result <- decisionM
    } else
    {
      result <- list("DM" = decisionM, "SPP" = safetyPP, "FPP" = futilPP)
    }
    return(result)
}


# fun_TITE_STEIN_decision(n = 3, pT = 0.3, phi1 = 0.75*0.3, phi2 = 1.25*0.3, psi1 = 0.3, psi2 = 0.8, cutoff_tox = 0.95, cutoff_eff = 0.98)
# fun_TITE_STEIN_decision(n = 6, pT = 0.3, phi1 = 0.75*0.3, phi2 = 1.25*0.3, psi1 = 0.3, psi2 = 0.8, cutoff_tox = 0.95, cutoff_eff = 0.98)
# fun_TITE_STEIN_decision(n = 9, pT = 0.3, phi1 = 0.75*0.3, phi2 = 1.25*0.3, psi1 = 0.3, psi2 = 0.8, cutoff_tox = 0.95, cutoff_eff = 0.98)
