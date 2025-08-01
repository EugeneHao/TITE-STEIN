library(tidyverse)
library(Rcpp)
library(RcppArmadillo)
library(inline)
library(Iso)
library(parallel)

source("fun_STEIN_OBD_new.R")
source("fun_STEIN_OBD.R")
source("fun_TITE_STEIN_core_para.R")
source("fun_TITE_STEIN_decision.R")
source("fun_TITE_STEIN_fixsimu.R")
source("fun_TITE_STEIN_update.R")
source("fun_TITE_STEIN.R")
source("fun_findOBD.R")
source("fun_pava.R")

pT = 0.3; qE = 0.25; pqcorr = 0; 
csize = 3; cN = 15;
repsize = 1000; n_cores = 8;
accrual = 10; tox_win = 30; eff_win = 90; susp = 0.5
tox_dist = "Uniform"; eff_dist = "Uniform"
tox_dist_hyper = NULL; eff_dist_hyper = NULL
u11 = 60; u00 = 40; 

# STEIN design parameters 
psi1 = 0.3; psi2 = 0.8; w1 = 0.33; w2 = 1.09;  
cutoff_tox = 0.95; cutoff_eff = 0.90;

# One simulation setting 
dN = 5; 
pV <- c(0.05, 0.10, 0.15, 0.30, 0.40)    # toxicity probability 
qV <- c(0.30, 0.50, 0.70, 0.75, 0.80)    # efficacy probability 

TITE_STEIN_new_result <- 
  fun_TITE_STEIN_fixsimu(dN, pV, qV, pT, pE, qE, pqcorr, 
                         csize, cN, design = "TITE-STEIN", utility = TRUE,
                         current = 1, doselimit = Inf, u11, u00, cutoff_tox, cutoff_eff,      
                         psi1, psi2, w1, w2,                                      # STEIN parameters
                         repsize, n_cores,        # replication parameters
                         accrual, susp, tox_win, eff_win, tox_dist, eff_dist,     # TITE parameters
                         OBDverify = TRUE)


trueOBD <- findOBD_RDS(pV, qV, pT, qE, u11, u00)

TITE_STEIN_new_result %>% 
  summarise(design = unique(design), 
            OBD99 = mean(OBD %in% c(-1, 0, 99), na.rm = T) * 100,
            OBDother = mean(OBD == trueOBD, na.rm = T) * 100,
            p1 = mean(OBD == 1, na.rm = T) * 100, p2 = mean(OBD == 2, na.rm = T) * 100,
            p3 = mean(OBD == 3, na.rm = T) * 100, p4 = mean(OBD == 4, na.rm = T) * 100, 
            p5 = mean(OBD == 5, na.rm = T) * 100, 
            n1 = mean(`1`, na.rm = T), n2 = mean(`2`, na.rm = T), 
            n3 = mean(`3`, na.rm = T), n4 = mean(`4`, na.rm = T), 
            n5 = mean(`5`, na.rm = T), 
            dur = mean(duration, na.rm = T)/30) %>% 
  mutate(p_OBD = ifelse(trueOBD <= 0, OBD99, OBDother), 
         p_rej = pmax(100 - p1 - p2 - p3 - p4 - p5, 0)) %>% 
  select(-OBD99, -OBDother) %>% "["(,c("design", "dur", "p_OBD", paste0("p", 1:5), "p_rej",
                                       paste0("n", 1:5))) 
