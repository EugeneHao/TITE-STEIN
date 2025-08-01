# the core function for one replication study 
fun_TITE_STEIN_core_para <-
  function(plist, qlist, trueOBD, decisionM, pT, pE, qE, pqcorr, csize, cN, design, 
           cutoff_tox, cutoff_eff, current, doselimit, psi, w1, w2,
           accrual, susp, tox_win, eff_win, tox_dist, eff_dist, tox_dist_hyper, eff_dist_hyper,
           use_susp, accrual_random, OBDverify = TRUE, pgrad = 0.1, seed = 1234)
{
  Result <- fun_TITE_STEIN(seed, plist, qlist, trueOBD, pT, psi, w1, w2, cutoff_tox, cutoff_eff, 
                           csize, cN, decisionM, current, doselimit,
                           futil_threshold = qE, alphaE = 1, betaE = 1, 
                           accrual, susp, tox_win, eff_win, tox_dist, eff_dist, tox_dist_hyper, eff_dist_hyper,   # for TITE
                           use_susp, accrual_random, OBDverify, pgrad)
  
  return(Result)
}