fun_TITE_STEIN_fixsimu <-
  function(dN, pV, qV, pT, pE, qE, pqcorr = 0, csize, cN, design, utility = TRUE,
           current = 1, doselimit = Inf,
           u11 = 100, u00 = 0, cutoff_tox = 0.95, cutoff_eff = 0.9,       
           psi1 = 0.35, psi2 = 0.65, w1 = 0.33, w2 = 1.09,                # STEIN design parameters
           repsize = 10000, n_cores = 10,                                 # replication parameters
           accrual = 10, susp = 0.5, tox_win = 30, eff_win = 60, tox_dist = "UnifCateg",   #TITE parameters
           eff_dist = "UnifCateg", tox_dist_hyper = NULL, eff_dist_hyper = NULL, use_susp = TRUE, 
           accrual_random = FALSE, OBDverify = TRUE, pgrad = 0.1)
  {
    if(utility == FALSE)
    {
      trueOBD <- findOBD(pV, qV, pT, qE)
    } else
    {
      trueOBD <- findOBD_RDS(pV, qV, pT, qE, u11, u00)
    }
    
    phi1 = 0.75*pT
    phi2 = 1.25*pT
    psi = log((1-psi1)/(1-psi2))/log(psi2*(1-psi1)/psi1/(1-psi2))
    phi_L = log((1-phi1)/(1-pT))/log((pT*(1-phi1))/(phi1*(1-pT))) # lower bound
    phi_U = log((1-pT)/(1-phi2))/log((phi2*(1-pT))/(pT*(1-phi2))) # upper bound
    
    # Get decision matrix for each possible sample size n
    decisionM = list()
    for (i in 1:cN) {
      decisionM[[i]] <- fun_TITE_STEIN_decision(n = i*csize, pT, phi1, phi2, psi1, psi2, cutoff_tox, cutoff_eff)
    }
    
    pM <- kronecker(matrix(pV, nrow = 1), rep(1, repsize))
    qM <- kronecker(matrix(qV, nrow = 1), rep(1, repsize))
    
    # replicate the simulation study
    ResultDF <- mclapply(1:repsize, FUN = function(x)
      fun_TITE_STEIN_core_para(pV, qV, trueOBD, decisionM, pT, pE, qE, pqcorr, csize, cN, design,
                               cutoff_tox, cutoff_eff, current, doselimit, psi, w1, w2, 
                               accrual, susp, tox_win, eff_win, tox_dist, eff_dist, tox_dist_hyper, eff_dist_hyper,
                               use_susp, accrual_random, OBDverify, pgrad, seed = x),
      mc.cores = n_cores) %>%
      do.call(rbind, .)
    
    return(data.frame(n = csize * cN,
                      design = design,
                      trueOBD = trueOBD[1],
                      OBD = ResultDF$OBD,
                      rN = ResultDF$rN,
                      select_OBD = ResultDF$select_OBD,
                      num_at_OBD = ResultDF$num_at_OBD,
                      risk_allocate = ResultDF$risk_allocate,
                      num_overdose_OBD = ResultDF$num_overdose_OBD,
                      num_overdose_nOBD = ResultDF$num_overdose_nOBD,
                      duration = ResultDF$duration                            # TITE
    ) %>%
      cbind(., ResultDF[, (ncol(ResultDF) - dN + 1):ncol(ResultDF)])
    )
  }
