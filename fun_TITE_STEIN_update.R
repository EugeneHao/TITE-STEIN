fun_TITE_STEIN_update <- function(cid, patDT, doseDT, current, time_current, tox_win, eff_win, csize,
                                   tox_dist, eff_dist, accrual, susp = 0.5,
                                   tox_dist_hyper = NULL, eff_dist_hyper = NULL, use_susp = TRUE, accrual_random = FALSE)
{
  dN <- nrow(doseDT)
  
  ### Step 1: update the cid'th cohort information in patDT #####
  ###### Accrual (1.1) ######################
  #   1.1 enrollment time for each patient
  if(accrual_random == TRUE)  # can also consider exponential distribution!
  {
    patDT$enroll[patDT$cid == cid] <- time_current + cumsum(runif(csize, 0, 2 * accrual)) + 1    # add 1 here
  } else
  {
    patDT$enroll[patDT$cid == cid] <- time_current + (0:(csize - 1)) * accrual + 1    # add 1 here
  }
  #   1.2 toxicity and efficacy assessment end time
  patDT$toxend[patDT$cid == cid] <- patDT$enroll[patDT$cid == cid] + tox_win
  patDT$effend[patDT$cid == cid] <- patDT$enroll[patDT$cid == cid] + eff_win
  #   1.3 toxicity and efficacy result at the end (0 or 1)
  patDT$toxobs[patDT$cid == cid] <- sample(c(0,1), size = csize, replace = T, prob = c(1-doseDT$tox[current], doseDT$tox[current]))
  patDT$effobs[patDT$cid == cid] <- sample(c(0,1), size = csize, replace = T, prob = c(1-doseDT$eff[current], doseDT$eff[current]))
  ######## DLT and Response Time (1.4)  Very Important #######
  #   1.4 generate the toxicity and efficacy response time
  new_tox = sum(patDT$cid == cid & patDT$toxobs == 1)   # number of DLT in the current cohort
  new_eff = sum(patDT$cid == cid & patDT$effobs == 1)   # number of response in the current cohort
  if(tox_dist == "Uniform") {
    patDT$toxdat[patDT$cid == cid & patDT$toxobs == 1] <- patDT$enroll[patDT$cid == cid & patDT$toxobs == 1] + runif(new_tox, max = tox_win)
  } else if(tox_dist == "UnifCateg") {
    tox_set <- 1:(tox_win/5-1)
    patDT$toxdat[patDT$cid == cid & patDT$toxobs == 1] <- patDT$enroll[patDT$cid == cid & patDT$toxobs == 1] + sample(tox_set, size = new_tox, replace = T) * 5
  }
  
  if(eff_dist == "Uniform") {
    patDT$effdat[patDT$cid == cid & patDT$effobs == 1] <- patDT$enroll[patDT$cid == cid & patDT$effobs == 1] + runif(new_eff, max = eff_win)
  } else if(eff_dist == "UnifCateg") {
    eff_set <- 1:(eff_win/5)
    patDT$effdat[patDT$cid == cid & patDT$effobs == 1] <- patDT$enroll[patDT$cid == cid & patDT$effobs == 1] + sample(eff_set, size = new_eff, replace = T) * 5
  } else if(eff_dist == "Case") {
    patDT$effdat[patDT$cid == cid & patDT$effobs == 1] <- patDT$enroll[patDT$cid == cid & patDT$effobs == 1] + 
      runif(new_eff, max = 28) + 28 * sample(0:2, size = new_eff, replace = T, prob = c(0.6, 0.3, 0.1))
  }
  
  
  #   1.5 update the dose level of the current cohort
  patDT$d[patDT$cid == cid] <- current
  
  #   1.6 update tox_confirm, eff_confirm
  patDT$tox_confirm[patDT$cid == cid] <- pmin(patDT$toxend[patDT$cid == cid], patDT$toxdat[patDT$cid == cid])
  patDT$eff_confirm[patDT$cid == cid] <- pmin(patDT$effend[patDT$cid == cid], patDT$effdat[patDT$cid == cid])
  
  
  ### Step 2: derive time_next ####
  
  tmp <- patDT %>% filter(d == current)
  n_d <- nrow(tmp)                 # number of patients assigned to dose level d (including current cohort)
  
  if(use_susp == TRUE)   # consider suspension rule
  {
    min_n <- min(floor(n_d * susp + 1), n_d)        # susp default value = 0.5 in (0, 1]; if susp = 1, it is the same as non-TITE designs
    tox_next <- sort(tmp$tox_confirm)[min_n]
    eff_next <- sort(tmp$eff_confirm)[min_n]
    if(accrual_random == TRUE)
    {
      time_next <- max(tox_next, eff_next, max(patDT$enroll, na.rm = T) + runif(1, 0, 2 * accrual) + 1)
    } else
    {
      time_next <- max(tox_next, eff_next, time_current + csize * accrual + 1)
    }
    # the earliest time fitting the accrual suspension rule and also the first patient in the next cohort is ready for enrollment
  } else
  {
    if(accrual_random == TRUE)
    {
      time_next <- max(patDT$enroll, na.rm = T) + runif(1, 0, 2 * accrual) + 1
    } else
    {
      time_next <- time_current + csize * accrual + 1
    }
  }
  
  ### Step 3: update doseDT ####
  #  (Note: when do we make the dose decision for the next cohort? )
  #  (We wait until the accrual suspension rule is satisfied and the next patient is available)
  
  # patDT_current initialization (generate at the beginning of each cohort assignment)
  # delta_t: toxicity has been ascertained or pending (0: no DLT; 1: DLT; -1: ascertained);
  # delta_e: efficacy has been ascertained or pending (0: no response; 1: response; -1: ascertained)
  # t_i: from dosing time to time_next
  # w_t: weight adjusting for the DLT has not yet been ascertained
  # w_e: weight adjusting for the response has not yet been ascertained
  # tox_exp: expectation of DLT at the current time (1 = observed, 0 = no DLT, (0, 1): expectation)
  # eff_exp: expectation of response at the current time
  patDT_current <- patDT %>% filter(!is.na(d)) %>%
    mutate(delta_t = -1, delta_e = -1, t_i = time_next - enroll)
  
  # 3.1  delta_t, delta_e (for each patient)
  patDT_current$delta_t[patDT_current$toxdat <= time_next] <- 1     # DLT observed
  patDT_current$delta_t[patDT_current$toxend <= time_next & patDT_current$toxdat > time_next] <- 0 # no DLT (DLT assessment end before time_next, and no DLT observed)
  # if delta_t = -1, the toxicity assessment has not finished at time_next
  
  patDT_current$delta_e[patDT_current$effdat <= time_next] <- 1     # response observed
  patDT_current$delta_e[patDT_current$effend <= time_next & patDT_current$effdat > time_next] <- 0 # no response
  # if delta_e = -1, the efficacy assessment has not finished at time_next
  
  ####  weight adjusting accounts for the partial information (3.2) very important #######
  # 3.2  w_t, w_e (for each patient)
  # remember this part can be changed if we do not assume uniform distribution for the event on [0, A_q]
  patDT_current <- patDT_current %>%
    mutate(w_t = ifelse(t_i < tox_win & delta_t < 0, t_i/tox_win, 0),     # t_i/A_t
           w_e = ifelse(t_i < eff_win & delta_e < 0, t_i/eff_win, 0))     # t_i/A_e
  
  
  
  # 3.3 ESS_t, ESS_e, pi_t_hat, pi_e_hat (for each dose level, should be updated on doseDT)
  # ESS_t: effective sample size for toxicity
  # ESS_e: effective sample size for efficacy
  # pi_t_hat, pi_e_hat: estimate of pi_t and pi_e (dose level)
  # x_d: quasi-event
  for(i in 1:dN)
  {
    if(sum(patDT_current$d == i) > 0)   # at least one cohort has been assigned to dose i
    {
      # update n, x, y
      doseDT$n[i] <- sum(patDT_current$d == i)
      delta_t_i = patDT_current$delta_t[patDT_current$d == i]
      delta_e_i = patDT_current$delta_e[patDT_current$d == i]
      
      doseDT$x[i] <- sum(delta_t_i == 1)    # number of observed DLT
      doseDT$y[i] <- sum(delta_e_i == 1)    # number of observed responses
      
      # update the effective sample size
      w_t_i = patDT_current$w_t[patDT_current$d == i]
      w_e_i = patDT_current$w_e[patDT_current$d == i]
      
      doseDT$ESS_t[i] <- sum(delta_t_i == 1) + sum(delta_t_i == 0) + sum(w_t_i)
      doseDT$pi_t_hat[i] <- doseDT$x[i]/doseDT$ESS_t[i]     # update unless n[i] > 0
      
      doseDT$ESS_e[i] <- sum(delta_e_i == 1) + sum(delta_e_i == 0) + sum(w_e_i)
      doseDT$pi_e_hat[i] <- doseDT$y[i]/doseDT$ESS_e[i]     # update unless n[i] > 0
    }
  }
  
  # 3.4 tox_exp, eff_exp (for each patient)
  patDT_current <- patDT_current %>% mutate(tox_exp = delta_t, eff_exp = delta_e)
  
  patDT_current$tox_exp[patDT_current$delta_t == -1] <-
    (doseDT$pi_t_hat[patDT_current$d] * (1 - patDT_current$w_t) / (1 - doseDT$pi_t_hat[patDT_current$d] * patDT_current$w_t))[patDT_current$delta_t == -1]
  
  patDT_current$eff_exp[patDT_current$delta_e == -1] <-
    (doseDT$pi_e_hat[patDT_current$d] * (1 - patDT_current$w_e) / (1 - doseDT$pi_e_hat[patDT_current$d] * patDT_current$w_e))[patDT_current$delta_e == -1]
  
  # 3.5 x_d (for each dose level)
  for(i in 1:dN) {
    if(sum(patDT_current$d == i) > 0)   # at least one cohort has been assigned to dose i
    {
      tox_exp_i <- patDT_current$tox_exp[patDT_current$d == i]
      eff_exp_i <- patDT_current$eff_exp[patDT_current$d == i]
    }
  }
  
  ### Step 4 return the outputs #####
  return(list(patDT = patDT, doseDT = doseDT, time_next = time_next ))
  
}



