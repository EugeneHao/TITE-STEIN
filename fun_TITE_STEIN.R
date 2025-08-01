#--------------------------Main function for STEIN design----------------------#

fun_TITE_STEIN <- 
  function(seed, plist, qlist, trueOBD, pT, psi, w1, w2, cutoff_tox, cutoff_eff, 
           csize, cN, decisionM, current = 1, doselimit = Inf,
           futil_threshold = 0.25, alphaE = 1, betaE = 1, 
           accrual, susp, tox_win, eff_win, tox_dist, eff_dist, tox_dist_hyper, eff_dist_hyper,   # for TITE
           use_susp = TRUE, accrual_random = FALSE, OBDverify = TRUE, pgrad = 0.1)
{
  dN <- length(plist)
  idlist <- 1:dN
  
  time_current <- 0                # the current time for assigning the next cohort
  time_next <- NA
  set.seed(seed)
  
  phi1 = 0.75*pT
  phi2 = 1.25*pT
  phi_L = log((1-phi1)/(1-pT))/log((pT*(1-phi1))/(phi1*(1-pT))) # lower bound
  phi_U = log((1-pT)/(1-phi2))/log((phi2*(1-pT))/(pT*(1-phi2))) # upper bound
  
  # patDT initialization (update at the beginning of each cohort assignment)
  patDT <- data.frame(cid = rep(1:cN, each = csize), pid = rep(1:csize, cN), id = 1:(csize * cN),
                      enroll = NA, toxend = NA, effend = NA, toxobs = -1, effobs = -1,
                      toxdat = 9999, effdat = 9999, tox_confirm = NA, eff_confirm = NA, d = NA)
  
  # doseDT initialization (update at the beginning of each cohort assignment)
  doseDT <- data.frame(id = 1:dN, tox = plist, eff = qlist, n = rep(0, dN),
                       x = rep(0,dN), y = rep(0, dN), keep = rep(1, dN),
                       ESS_t = 0, ESS_e = 0, pi_t_hat = 0, pi_e_hat = 0
  )
  
  earlystop <- FALSE
  record <- rep(-1, cN)         # record dose selection
  
  for (i in 1:cN) {
    TITE_one <- fun_TITE_STEIN_update(cid = i, patDT, doseDT, current, time_current, tox_win,
                                       eff_win, csize, tox_dist, eff_dist, accrual, susp,
                                       tox_dist_hyper, eff_dist_hyper,
                                       use_susp, accrual_random)
    patDT <- TITE_one$patDT
    doseDT <- TITE_one$doseDT
    time_current <- TITE_one$time_next
    
    
    TITE_STEIN_onestep <- 
      TITE_STEIN_one(doseDT, current, pT, psi, phi_L, phi_U, csize, decisionM, futil_threshold,
                     alphaE, betaE, cutoff_eff = cutoff_eff)
    
    doseDT <- TITE_STEIN_onestep$doseDT
    record[i] <- current
    current <- TITE_STEIN_onestep$newdose   # next cohort dose 
    # check whether the next dose exists or not 
    if(is.na(current)) {
      earlystop = TRUE
      OBD = 99
      break
    } else if(doseDT$n[current] + csize > doselimit) {  # check whether the next dose exceeds the dose limit or not 
      break   # not early stop 
    } 
  }
  
  # after the trial stop, derive the doseDT_final for OBD selection
  doseDT_final <- doseDT %>% select(id, tox, eff, n, keep) %>%
    left_join(., patDT %>% filter(!is.na(d)) %>% group_by(d) %>%
                summarise(x = sum(toxobs == 1), y = sum(effobs == 1)) %>% "colnames<-"(c("id", "x", "y")))
  doseDT_final$x[is.na(doseDT_final$x)] <- 0
  doseDT_final$y[is.na(doseDT_final$y)] <- 0
  
  # OBD selection
  if(earlystop == FALSE) {
    if(OBDverify == FALSE)
    {
      OBD <- fun_STEIN_OBD(doseDT_final, pT, w1, w2)
    } else
    {
      OBD <- fun_STEIN_OBD_new(doseDT_final, pT, psi1 = 0.3, w1, w2, size = 1000, pgrad)
    }

    DTremain = doseDT[doseDT$keep == 1 & doseDT$n > 0, ]  # need to select doses applied to patients
    rN = nrow(DTremain)   # how many number of doses remain for selection
  } else {
    OBD = 99   # 99: early stop
    rN = 0
  }
  
  trueOBDone = trueOBD[1]   # the best one 
  
  # with one OBD
  select_OBD <- ifelse(((trueOBDone %in% idlist) & (OBD %in% idlist) & (OBD == trueOBDone))  | ((!trueOBDone %in% idlist) & (!OBD %in% idlist)), 1, 0)
  num_at_OBD <- ifelse(trueOBDone %in% idlist, doseDT[trueOBDone, 4], NA)
  risk_allocate <- ifelse(trueOBDone %in% idlist, 
                          ifelse(doseDT[trueOBDone, 4] < csize * cN/5, 1, 0), NA)
  
  num_overdose_OBD <- ifelse(trueOBDone %in% idlist, ifelse(max(plist) > (pT + 0.1), sum(doseDT[which(plist > (pT+0.1)), 4]), 0), NA)
  num_overdose_nOBD <- ifelse(!trueOBDone %in% idlist, ifelse(max(plist) > (pT + 0.1), sum(doseDT[which(plist > (pT+0.1)), 4]), 0), NA)
  
  # TITE statistics
  patDT <- patDT %>% filter(!is.na(tox_confirm))
  patDT$tox_confirm[patDT$cid == max(patDT$cid)] <- patDT$toxend[patDT$cid == max(patDT$cid)]  # end the trial until all the patient finish their assessments
  patDT$eff_confirm[patDT$cid == max(patDT$cid)] <- patDT$effend[patDT$cid == max(patDT$cid)]
  
  duration <- max(c(patDT$tox_confirm, patDT$eff_confirm), na.rm = T)   # the end of the trial
  
  # not include doseDT and record here
  return(data.frame(earlystop = earlystop, OBD = OBD, rN = rN, trueOBD = trueOBDone, 
                    select_OBD = select_OBD, num_at_OBD = num_at_OBD, 
                    num_overdose_OBD = num_overdose_OBD, 
                    num_overdose_nOBD = num_overdose_nOBD,
                    risk_allocate = risk_allocate, 
                    duration = duration) %>% 
           cbind(., t(data.frame(doseDT$n)))
  )
}

#------------------------------------------------------------------------------#


TITE_STEIN_one <- function(doseDT, current, pT, psi, phi_L, phi_U, csize, decisionM, 
                           futil_threshold, alphaE=1, betaE=1, cutoff_eff) {
  
  n = doseDT$n[current]
  x = doseDT$x[current]    # observed DLT
  y = doseDT$y[current]    # observed efficacy response
  cn <- round(n/csize)                     # number of cohort assigned to the dose level
  dN <- nrow(doseDT)
  
  onedecisionM <- decisionM[[cn]]
  dT <- onedecisionM$dT      # three columns: n_T1, D, DU
  dE <- onedecisionM$dE      # three columns: n_E1, TBD, EU
  
  above = ifelse(current == dN, NA, current + which(doseDT$keep[(current+1):dN] == 1)[1]) # the next available higher dose, can be NA if no available dose
  below = ifelse(current == 1, NA, current - which(doseDT$keep[(current-1):1] == 1)[1]) # the next available lower dose, can be NA if no available dose
  
  n_T0_tilde <- doseDT$ESS_t[current] - doseDT$x[current]
  n_E0_tilde <- doseDT$ESS_e[current] - doseDT$y[current]
  
  p_hat <- x/(x + n_T0_tilde)   # imputed DLT rate 
  
  DU_limit <- dT$DU[x+1]       # maximum n_T0_tilde to get DU  (-Inf means DU is not possible)
  EU_limit <- dE$EU[y+1]       # minimum n_E0_tilde to get EU  ( Inf means EU is not possible)
  D_limit <- dT$D[x+1]         # maximum n_T0_tilde to get D   (-Inf means D is not possible)
  TBD_limit <- dE$TBD[y+1]     # minimum n_E0_tilde to get E   ( Inf means E is not possible)
  
  if(n_T0_tilde <= DU_limit)                    # DU_T
  {
    doseDT$keep[current:dN] = 0                 # remove current and above doses
    newdose <- below
  } else if(n_T0_tilde <= D_limit & n_E0_tilde >= EU_limit)              # DU_E
  {
    doseDT$keep[current] = 0            # remove current dose
    newdose <- below
  } else if(n_T0_tilde <= D_limit)      # D
  {
    newdose <- ifelse(is.na(below), current, below)  # stay if no available dose below
  } else if(n_E0_tilde >= EU_limit)     # EU
  {
    doseDT$keep[current] = 0                    # remove current dose
    newdose <- ifelse(is.na(above), below, above)
  } else if(n_E0_tilde >= TBD_limit)      # TBD (STEIN special)
  {
    # calculate futility posterior prob at current dose
    futilPP = pbeta(futil_threshold, y+alphaE, n_E0_tilde+betaE)
    
    # posterior prob for current dose
    if(futilPP > cutoff_eff) {
      post_current = -Inf # current dose will not be considered due to futility
      doseDT$keep[current] = 0  # remove current dose due to futility
    } else {
      post_current = 1 - pbeta(psi, y+1, n_E0_tilde + 1) 
    }
    
    # posterior prob for the dose below
    if(is.na(below)) {
      post_below = -Inf
    } else {
      post_below = 1 - pbeta(psi, doseDT$y[below]+1, doseDT$ESS_e[below] - doseDT$y[below] + 1)
    }
    
    # posterior prob for the dose above
    if(is.na(above)) {
      post_above = -Inf
    } else {
      post_above = 1 - pbeta(psi, doseDT$y[above]+1, doseDT$ESS_e[above] - doseDT$y[above] + 1)
    }
    
    # evaluate posterior probs for each dose in the local admissible dose set
    if(p_hat <= phi_L) { # yellow region, we consider three possible doses {j-1, j, j+1}
      post = c(post_below, post_current, post_above) + c(0, 0.001, 0.002)
      select = c(below, current, above)
      newdose = select[which.max(post)]
    } else { # orange region, we only consider two possible doses {j-1, j}
      post = c(post_below, post_current) + c(0, 0.001)
      select = c(below, current)
      newdose = select[which.max(post)]
    }
  } else                                # S
  {
    newdose <- current     # stay
  }
  
  return(list(doseDT = doseDT, newdose = newdose))
}

#------------------------------------------------------------------------------#


