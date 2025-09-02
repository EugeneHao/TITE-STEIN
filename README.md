# TITE-STEIN: A Time-to-Event Model-Assisted Design for Phase I/II Trials

**TITE-STEIN** (Time-to-event Simple Toxicity and Efficacy Interval Design)[^1] is a model-assisted dose-finding design developed to accelerate Phase I/II clinical trials by incorporating time-to-event toxicity and efficacy outcomes. It extends the STEIN framework to handle late-onset outcomes and introduces an **Optimal Biological Dose (OBD) verification** procedure to mitigate the risk of selecting inadmissible doses.

## Design Features
- Incorporates both toxicity and efficacy outcomes for OBD selection.
- Supports real-time dose decision-making using late-onset toxicity and efficacy outcomes.
- Includes an OBD verification step to avoid selecting an unsafe or futile dose as the OBD.
- Straightforward to implement with no complex dose-toxicity and dose-efficacy modeling assumptions.

## User Guide

`run_TITE-STEIN.R`: Main function for conducting a simulation study using TITE-STEIN design with replications.

```R
# Preparation

library(tidyverse)
library(Rcpp)
library(RcppArmadillo)
library(inline)
library(Iso)
library(parallel)

# Load all the other functions 

source("fun_STEIN_OBD_new.R")
source("fun_STEIN_OBD.R")
source("fun_TITE_STEIN_core_para.R")
source("fun_TITE_STEIN_decision.R")
source("fun_TITE_STEIN_fixsimu.R")
source("fun_TITE_STEIN_update.R")
source("fun_TITE_STEIN.R")
source("fun_findOBD.R")
source("fun_pava.R")
```
```R
# Design Settings

pT = 0.3;       # maximum acceptable toxicity probability 
qE = 0.25;      # minimum acceptable efficacy probability 
pqcorr = 0;     # correlation between the toxicity and efficacy probabilities 
csize = 3;      # sample size per cohort 
cN = 15;        # maximum number of cohorts 
repsize = 1000;      # replication size 
n_cores = 8;         # number of cores for parallel computation 
accrual = 10;        # enroll 1 patient per 10 days 
tox_win = 30;        # toxicity assessment window = 30 days 
eff_win = 90;        # efficacy assessment window = 90 days 
susp = 0.5;          # suspend accrual until 50% toxicity and efficacy results observed 

tox_dist = "Uniform";  # time to toxicity outcome follows a uniform distribution 
eff_dist = "Uniform";  # time to efficacy outcome follows a uniform distribution 
tox_dist_hyper = NULL; # time to toxicity outcome does not have other distribution parameters 
eff_dist_hyper = NULL; # time to efficacy outcome does not have other distribution parameters 
u11 = 60;    # utility of having both DLT and efficacy response
u00 = 40;    # utility of no DLT and no efficacy response

# STEIN design parameters 
psi1 = 0.3;      # null hypothesis of efficacy probability 
psi2 = 0.8;      # alternative hypothesis of efficacy probability 
w1 = 0.33; w2 = 1.09;    # weights in the utility function for OBD selection 
cutoff_tox = 0.95;    	 # cutoff probability for safety elimination rule 
cutoff_eff = 0.90;       # cutoff probability for futility elimination rule 

# Probability settings
dN = 5;    # number of doses 
pV = c(0.05, 0.10, 0.15, 0.30, 0.40)    # toxicity probability 
qV = c(0.30, 0.50, 0.70, 0.75, 0.80)    # efficacy probability 
```

```R
TITE_STEIN_result = 
  fun_TITE_STEIN_fixsimu(dN, pV, qV, pT, pE, qE, pqcorr, 
                         csize, cN, design = "TITE-STEIN", utility = TRUE,
                         current = 1, doselimit = Inf, u11, u00, 
                         cutoff_tox, cutoff_eff,  
                         psi1, psi2, w1, w2,     
                         repsize, n_cores,       
                         accrual, susp, tox_win, eff_win, tox_dist, eff_dist, 
                         use_susp = TRUE, 
                         accrual_random = FALSE,
                         OBDverify = TRUE)


```

> Other simulation settings in `fun_TITE_STEIN_fixsimu` :  
>
> + `utility = TRUE`: use utility score to determine the true OBD; otherwise, the OBD is the dose with the highest efficacy probability among those with acceptable toxicity. If multiple doses share the highest utility score or highest efficacy probability, the OBD is the dose with the lowest toxicity probability. 
> + `current = 1`: the start dose is the dose level 1, i.e., the lowest dose level.  
> + `doselimit = Inf`: the trial will end when the maximum sample size is reached or no suitable dose for the next patient cohort during the trial. Users can let the trial stop when a dose has enrolled a specific number of patients by setting `doselimit` to that number. 
> + `use_susp = TRUE`: apply accrual suspension rule 
> + `accrual_random = FALSE`:  enroll patient every `accrual` days; if set as true, then the accrual time follows a uniform distribution U(0, 2`accrual`)
> + `OBDverify = TRUE`: include the OBD verification step. 

## Simulation Results

```R
head(TITE_STEIN_result)

           n     design trueOBD OBD rN select_OBD num_at_OBD risk_allocate num_overdose_OBD
doseDT.n  45 TITE-STEIN       3   3  3          1         36             0                0
doseDT.n1 45 TITE-STEIN       3   2  2          0          0             1                0
doseDT.n2 45 TITE-STEIN       3   3  4          1         24             0                0
doseDT.n3 45 TITE-STEIN       3   3  3          1         18             0                0
doseDT.n4 45 TITE-STEIN       3   3  3          1         18             0                0
doseDT.n5 45 TITE-STEIN       3   2  3          0          9             0                0
          num_overdose_nOBD duration  1  2  3  4 5
doseDT.n                 NA 755.2572  6  3 36  0 0
doseDT.n1                NA 700.0325  9 36  0  0 0
doseDT.n2                NA 836.0901  3  3 24 15 0
doseDT.n3                NA 794.3676 12 15 18  0 0
doseDT.n4                NA 746.6938  9 18 18  0 0
doseDT.n5                NA 798.0410 27  9  9  0 0
```

+ `n`: maximum sample size.
+ `trueOBD`: the true OBD given the true toxicity and efficacy probability values based on utility score or maximum efficacy probability; `trueOBD = 0` if all doses are toxic; `trueOBD = -1` if all safe doses are futile. 
+ `OBD`: the identified OBD through the simulation study.
+ `rN`: number of remaining doses for OBD selection.
+ `select_OBD`: set to 1 if `OBD = trueOBD`; otherwise, set to 0
+ `num_at_OBD`: number of patients assigned to the true OBD; if `trueOBD <= 0`, then `num_at_OBD = NA`.  
+ `risk_allocate`: equal to 1 if we assign less than 20% patients to the true OBD; otherwise, set to 0. `risk_allocate = NA` if the true OBD does not exist. 
+ `num_overdose_OBD`: number of patients assigned to doses with toxicity probability greater than pT + 0.1 when the true OBD exists; set to NA if the true OBD does not exist. 
+ `num_overdose_nOBD`: number of patients assigned to doses with toxicity probability greater than pT + 0.1 when the true OBD does not exist; set to NA if the true OBD exists.
+ `duration`: trial duration in days.
+ `1`, `2`, ..., `d`: number of patients assigned to each dose level.



We can summarize the replicated simulation results through the code below:  

```R
trueOBD = findOBD_RDS(pV, qV, pT, qE, u11, u00)

TITE_STEIN_new_result %>% 
  summarise(design = unique(design), 
            OBD99 = mean(OBD %in% c(-1, 0, 99), na.rm = T) * 100,
            OBDother = mean(OBD == trueOBD, na.rm = T) * 100,
            p1 = mean(OBD == 1, na.rm = T) * 100, 
            p2 = mean(OBD == 2, na.rm = T) * 100,
            p3 = mean(OBD == 3, na.rm = T) * 100, 
            p4 = mean(OBD == 4, na.rm = T) * 100, 
            p5 = mean(OBD == 5, na.rm = T) * 100, 
            n1 = mean(`1`, na.rm = T), 
            n2 = mean(`2`, na.rm = T), 
            n3 = mean(`3`, na.rm = T), 
            n4 = mean(`4`, na.rm = T), 
            n5 = mean(`5`, na.rm = T), 
            dur = mean(duration, na.rm = T)/30) %>% 
  mutate(p_OBD = ifelse(trueOBD <= 0, OBD99, OBDother), 
         p_rej = pmax(100 - p1 - p2 - p3 - p4 - p5, 0)) %>% 
  select(-OBD99, -OBDother) %>% 
  "["(,c("design", "dur", "p_OBD", paste0("p", 1:5), "p_rej", paste0("n", 1:5))) 

```

```R
      design      dur p_OBD  p1   p2 p3  p4  p5 p_rej    n1     n2     n3    n4    n5
1 TITE-STEIN 25.41886    67 1.7 21.9 67 9.1 0.2   0.1 6.429 14.127 20.592 3.711 0.114
```

+ `dur`: average trial duration in months.
+ `p_OBD`: probability of correctly selecting the true OBD.
+ `p1`, ..., `pd`: probability of selecting each dose as the true OBD.
+ `p_rej`: probability of rejecting the trial without selecting any dose level as the OBD.
+ `n1`, ..., `nd`: number of patients assigned to each dose level.

# Case STUDY

This case study uses data from the TRANSCEND NHL 001 trial of lisocabtagene maraleucel (liso-cel) for patients with relapsed or refractory B-cell lymphomas [^2]. Patients were assigned to one of the three traget dose levels: 

+ DL1: $50 \times 10^6$ CAR+ T cells 
+ DL2: $100 \times 10^6$ CAR+ T cells 
+ DL3: $150 \times 10^6$ CAR+ T cells 

> The DL1D (two doses of $50 \times 10^6$ CAR$^+$ T cells given on day 1 and day 15) is not included because of the small sample size. 

```R
# Design Settings for the case study

pT = 0.3;       # maximum acceptable toxicity probability 
qE = 0.25;      # minimum acceptable efficacy probability 
pqcorr = 0;     # correlation between the toxicity and efficacy probabilities 
csize = 3;      # sample size per cohort 
cN = 45;        # maximum number of cohorts 
repsize = 1000;      # replication size 
n_cores = 8;         # number of cores for parallel computation 
accrual = 10;        # enroll 1 patient per 10 days 
tox_win = 28;        # toxicity assessment window = 30 days 
eff_win = 84;        # efficacy assessment window = 90 days 
susp = 0.5;          # suspend accrual until 50% toxicity and efficacy results observed 

tox_dist = "Uniform";  # time to toxicity outcome follows a uniform distribution 
eff_dist = "Uniform";  # time to efficacy outcome follows a uniform distribution 
tox_dist_hyper = NULL; # time to toxicity outcome does not have other distribution parameters 
eff_dist_hyper = NULL; # time to efficacy outcome does not have other distribution parameters 
u11 = 60;    # utility of having both DLT and efficacy response
u00 = 40;    # utility of no DLT and no efficacy response

# STEIN design parameters 
psi1 = 0.3;      # null hypothesis of efficacy probability 
psi2 = 0.8;      # alternative hypothesis of efficacy probability 
w1 = 0.33; w2 = 1.09;    # weights in the utility function for OBD selection 
cutoff_tox = 0.95;    	 # cutoff probability for safety elimination rule 
cutoff_eff = 0.90;       # cutoff probability for futility elimination rule 

# Probability settings
dN = 3;    # number of doses 
pV <- c(0.07, 0.10, 0.12)    # toxicity probability 
qV <- c(0.65, 0.75, 0.75)    # efficacy probability 
```



```R
# Run the case study
fun_TITE_STEIN_fixsimu(dN, pV, qV, pT, pE, qE, pqcorr, 
                       csize, cN, design = "TITE-STEIN", utility = TRUE,
                       current = 1, doselimit = Inf, u11, u00, 
                       cutoff_tox, cutoff_eff,  
                       psi1, psi2, w1, w2,     
                       repsize, n_cores,       
                       accrual, susp, tox_win, eff_win, tox_dist, eff_dist, 
                       use_susp = TRUE, 
                       accrual_random = FALSE,
                       OBDverify = TRUE) %>% 
  summarise(design = unique(design), 
            OBD99 = mean(OBD %in% c(-1, 0, 99), na.rm = T) * 100,
            OBDother = mean(OBD == trueOBD, na.rm = T) * 100,
            p1 = mean(OBD == 1, na.rm = T) * 100, 
            p2 = mean(OBD == 2, na.rm = T) * 100,
            p3 = mean(OBD == 3, na.rm = T) * 100, 
            n1 = mean(`1`, na.rm = T), 
            n2 = mean(`2`, na.rm = T), 
            n3 = mean(`3`, na.rm = T), 
            dur = mean(duration, na.rm = T)/30) %>% 
  mutate(p_OBD = ifelse(trueOBD <= 0, OBD99, OBDother), 
         p_rej = pmax(100 - p1 - p2 - p3, 0)) %>% 
  select(-OBD99, -OBDother) %>% 
  "["(,c("design", "dur", "p_OBD", paste0("p", 1:3), "p_rej", paste0("n", 1:3))) 
```

```R 
      design      dur p_OBD   p1   p2   p3 p_rej     n1     n2     n3
1 TITE-STEIN 30.68859  55.6 28.5 55.6 15.9     0 51.903 65.538 17.559
```



## Reference

[^1]: Sun, H., Tu, J., Ananthakrishnan, R., & Kim, E. (2025). TITE-STEIN: Time-to-event Simple Toxicity and Efficacy Interval Design to Accelerate Phase I/II Trials. *Journal of Biopharmaceutical Statistics* (under review).
[^2]: Abramson, Jeremy S., M. Lia Palomba, Leo I. Gordon, Matthew A. Lunning, Michael Wang, Jon Arnason, Amitkumar Mehta et al. "Lisocabtagene maraleucel for patients with relapsed or refractory large B-cell lymphomas (TRANSCEND NHL 001): a multicentre seamless design study." *The Lancet* 396, no. 10254 (2020): 839-852.
