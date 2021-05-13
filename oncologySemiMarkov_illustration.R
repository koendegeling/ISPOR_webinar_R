##----------------------------------------------------------------------------------------------------------------##
##                                                                                                                ##
##                      ISPOR STUDENT NETWORK WEBINAR: HEALTH ECONOMIC MODELLING IN R                             ##
##                                                                                                                ##
##                               by Petros Pechlivanoglou and Koen Degeling                                       ##
##                                                                                                                ##
##----------------------------------------------------------------------------------------------------------------##
#
# This script has been created by Petros Pechlivanougly and Koen Degeling for the ISPOR Student Network Webinar or
# Health Economic Modelling using R. It demonstrated the implementation on a three-state Markov model using three
# health states: progression-free (F), progression (P) and dead (D). This model structure is ofter referred to as
# the three-state oncology model and is often used to model the cost-effectiveness of advanced cancer treatments.
# The objective of this hypothetical cost-utility study is to estimate the cost per quality-adjusted life year
# gained for an experimental treatment strategy (Exp) compared to the standard of care (SoC). The experimental
# treatment reduces the risk of progression, but also is more expensive.
#
# The script is organized according to the following sections:
# 1. GLOBAL INITIALIZATION      setting R up for the analysis
# 2. SURVIVAL ANALYSIS          survival analysis of individual patient data
# 2. MARKOV MODEL               implementing and evaluating the Markov model
# 3. SENSITIVITY ANALYSIS       defining and performing deterministic sensitivity analyses
# 4. PROBABILISTIC ANALYSIS     defining and performing a probabilistic analysis of the model
#
# The script is heavily based on a script that has been developed by the Decision Analysis in R for Technologies in 
# Health (DARTH) workgroup. This scripts adds the survival analysis and turn the Markov model into a semi-Markov
# model with time-dependent transition probabilities. Please see the README.md file for details on the DARTH 
# workgroup and their work, and make sure to appropriately reference their publications when using their scripts.




### 1. GLOBAL INITIALIZATION ----

# Clear the Global Environment
rm(list = ls()) 

# Check whether the required packages are installed and, if not, install the missing packages
# - the 'pacman' package is used to conveniently install other packages
if (!require("pacman")) install.packages("pacman"); library(pacman)
p_load("flexsurv", "MASS", "dplyr", "devtools", "scales", "ellipse", "ggplot2", "lazyeval", "igraph", "ggraph", "reshape2", "knitr", "stringr", "diagram")   
p_load_gh("DARTH-git/dampack", "DARTH-git/darthtools")

# Load the hypothetical individual patient data for the time-to-progression (TTP)
load(file = "df_TTP.RData")




### 2. SURVIVAL ANALYSIS ----

# To demonstrate the implementation of a semi-Markov model with time-dependent transition probabilities, a 
# parametric survival distribution/model will be fitted to some hypothetical data for the time-to-progression (TTP) 
# for the standard of care. A hazard ratio for the experimental vs. the standard of care will then be applied to 
# obtain transition probabilities for the experimental strategy.

# 2.1 Parametric survival analysis ----

# Using the 'flexsurv' package, several commonly used parametric survival distributions will be fitted
# - note that the data includes a column for the time (in years) and a status indicator whether the time corresponds
#   to an event, i.e. progression (status = 1), or to the last time of follow up, i.e. censoring (status = 0)
head(df_TTP)

l_TTP_SoC_exp      <- flexsurvreg(formula = Surv(time, status) ~ 1, data = df_TTP, dist = "exp")
l_TTP_SoC_gamma    <- flexsurvreg(formula = Surv(time, status) ~ 1, data = df_TTP, dist = "gamma")
l_TTP_SoC_gompertz <- flexsurvreg(formula = Surv(time, status) ~ 1, data = df_TTP, dist = "gompertz")
l_TTP_SoC_llogis   <- flexsurvreg(formula = Surv(time, status) ~ 1, data = df_TTP, dist = "llogis")
l_TTP_SoC_lnorm    <- flexsurvreg(formula = Surv(time, status) ~ 1, data = df_TTP, dist = "lnorm")
l_TTP_SoC_weibull  <- flexsurvreg(formula = Surv(time, status) ~ 1, data = df_TTP, dist = "weibull")


## 2.2 Inspecting the fits ----

# Inspect fit based on visual fit
colors <- rainbow(6)
plot(l_TTP_SoC_exp,       col = colors[1], ci = FALSE, ylab = "Event-free proportion", xlab = "Time in years", las = 1)
lines(l_TTP_SoC_gamma,    col = colors[2], ci = FALSE)
lines(l_TTP_SoC_gompertz, col = colors[3], ci = FALSE)
lines(l_TTP_SoC_llogis,   col = colors[4], ci = FALSE)
lines(l_TTP_SoC_lnorm,    col = colors[5], ci = FALSE)
lines(l_TTP_SoC_weibull,  col = colors[6], ci = FALSE)
legend("right",
       legend = c("exp", "gamma", "gompertz", "llogis", "lnorm", "weibull"),
       col    = colors,
       lty    = 1,
       bty    = "n")

# Compare the fit numerically based on the AIC
(v_AIC <- c(
  exp      = l_TTP_SoC_exp$AIC,
  gamma    = l_TTP_SoC_gamma$AIC,
  gompertz = l_TTP_SoC_gompertz$AIC,
  llogis   = l_TTP_SoC_llogis$AIC,
  lnorm    = l_TTP_SoC_lnorm$AIC,
  weibull  = l_TTP_SoC_weibull$AIC
))

# Weibull has the best visual and numerical fit


## 2.3 Saving the survival parameters ----

# Note that the 'flexsurv' package return the coefficients, which need to be transformed for use in the base R
# functions, but that will be done when the coefficients actually are used
l_TTP_SoC_weibull
l_TTP_SoC_weibull$coefficients

coef_weibull_shape_SoC <- l_TTP_SoC_weibull$coefficients["shape"]
coef_weibull_scale_SoC <- l_TTP_SoC_weibull$coefficients["scale"]





### 3. SEMI-MARKOV MODEL ----

## 3.1 Model parameters ----

# Defining general information regarding the model states and treatment strategies
# - note the abbreviations for the strategies: standard of care (SoC), experimental (Exp)
# - note the abbreviations for the states: progression-free (F), progression (P), dead (D)
v_names_strats <- c("Standard of Care", "Experimental")         # strategy names
v_names_states <- c("ProgressionFree", "Progression", "Dead")   # state names
n_strats       <- length(v_names_strats)                        # number of strategies
n_states       <- length(v_names_states)                        # number of states

# General modelling parameters
t_cycle <- 1/4      # cycle length of 3 months (in years)                                   
n_cycle <- 60       # number of cycles (total runtime 15 years)
d_e     <- 0.03     # discount rate (per year) for effects
d_c     <- 0.03     # discounr rate (per year) for costs

# Time-dependent transition probabilities are obtained in four steps
# 1) Defining the cycle times
# 2) Obtaining the event-free (i.e. survival) probabilities for the cycle times for SoC
# 3) Obtaining the event-free (i.e. survival) probabilities for the cycle times for Exp based on a hazard ratio
# 4) Obtaining the time-dependent transition probabilities from the event-free (i.e. survival) probabilities

# 1) Defining the cycle times
(t <- seq(from = 0, by = t_cycle, length.out = n_cycle + 1))

# 2) Obtaining the event-free (i.e. survival) probabilities for the cycle times for SoC
#    - note that the coefficients need to be transformed to obtain the parameters that the base R function uses
S_FP_SoC <- pweibull(
  q     = t, 
  shape = exp(coef_weibull_shape_SoC), 
  scale = exp(coef_weibull_scale_SoC), 
  lower.tail = FALSE
)

head(cbind(t, S_FP_SoC))

# 3) Obtaining the event-free (i.e. survival) probabilities for the cycle times for Exp based on a hazard ratio
# - note that S(t) = exp(-H(t)) and, hence, H(t) = -ln(S(t))
HR_FP_Exp <- 0.6
H_FP_SoC  <- -log(S_FP_SoC)
H_FP_Exp  <- H_FP_SoC * HR_FP_Exp
S_FP_Exp  <- exp(-H_FP_Exp)

head(cbind(t, S_FP_SoC, H_FP_SoC, H_FP_Exp, S_FP_Exp))

# 4) Obtaining the time-dependent transition probabilities from the event-free (i.e. survival) probabilities
p_FP_SoC <- p_FP_Exp <- rep(NA, n_cycle)

for(i in 1:n_cycle) {
  p_FP_SoC[i] <- 1 - S_FP_SoC[i+1] / S_FP_SoC[i]
  p_FP_Exp[i] <- 1 - S_FP_Exp[i+1] / S_FP_Exp[i]
}

p_FP_SoC
p_FP_Exp

# Time-constant transition probabilities
p_FD      <- 0.02   # probability of dying when progression-free
p_PD      <- 0.1    # probability of dying when in progression

# Costs  
c_F_SoC   <- 400    # cost of one cycle in progression-free state, under the standard of care
c_F_Exp   <- 800    # cost of one cycle in progression-free state, under the experimental treatment
c_P       <- 1000   # cost of one cycle in progression state
c_D       <- 0      # cost of one cycle in dead state

# Health utility values
u_F       <- 0.8    # utility when healthy 
u_P       <- 0.5    # utility when sick
u_D       <- 0      # utility when dead


## 3.2 Markov cohort trace matrix ----

# Initialize matrices to store the Markov cohort traces for each strategy
# - note that the number of rows is n_cycle + 1, because R doesn't use index 0 (i.e. cycle 0)
m_M_SoC <- m_M_Exp <-  matrix(
  data = NA, 
  nrow = n_cycle + 1 ,  
  ncol = n_states, 
  dimnames = list(0:n_cycle, v_names_states)
)

# Specifying the initial state for the cohorts (all patients start in ProgressionFree)
m_M_SoC[1, "ProgressionFree"] <- m_M_Exp[1, "ProgressionFree"] <- 1
m_M_SoC[1, "Progression"]     <- m_M_Exp[1, "Progression"]     <- 0
m_M_SoC[1, "Dead"]            <- m_M_Exp[1, "Dead"]            <- 0

# Inspect whether properly defined
head(m_M_SoC)
head(m_M_Exp)


## 3.3 Transition probability matrix ----

# Initialize matrices or the transition probabilities
# - note that these are now 3-dimensional matrices
# - starting with standard of care
m_P_SoC <- array(
  data = 0,
  dim = c(n_states, n_states, n_cycle),
  dimnames = list(v_names_states, v_names_states, paste0("Cycle", 1:n_cycle))
)

# Setting the transition probabilities from ProgressionFree based on the model parameters
m_P_SoC["ProgressionFree", "ProgressionFree", ] <- (1 - p_FD) * (1 - p_FP_SoC)
m_P_SoC["ProgressionFree", "Progression", ]     <- (1 - p_FD) * p_FP_SoC
m_P_SoC["ProgressionFree", "Dead", ]            <- p_FD

# From Progression
m_P_SoC["Progression", "Progression", ] <- 1 - p_PD
m_P_SoC["Progression", "Dead", ]        <- p_PD

# From Dead
m_P_SoC["Dead", "Dead", ] <- 1

# Using the transition probabilities for standard of care as basis, update the transition probabilities that are
# different for the experimental strategy
m_P_Exp <- m_P_SoC
m_P_Exp["ProgressionFree", "ProgressionFree", ] <- (1 - p_FD) * (1 - p_FP_Exp)
m_P_Exp["ProgressionFree", "Progression", ]     <- (1 - p_FD) * p_FP_Exp

# Check that transition probabilities are in [0, 1] using a function from the 'dampack' package
check_transition_probability(m_P_SoC, verbose = TRUE)
check_transition_probability(m_P_Exp, verbose = TRUE)

# Check that all rows sum to 1 using a function from the 'dampack' package
check_sum_of_transition_array(m_P_SoC, n_states = n_states, n_cycles = n_cycle, verbose = TRUE)
check_sum_of_transition_array(m_P_Exp, n_states = n_states, n_cycles = n_cycle, verbose = TRUE)

# Inspect whether properly defined
# - note that we inspect for the first two cycles
m_P_SoC[ , , 1:2]
m_P_Exp[ , , 1:2]


## 3.4 Model evaluation ----

# Create the Markov cohort trace by looping over all cycles
# - note that the trace can easily be obtained using matrix multiplications
# - note that now the right probabilities for the cycle need to be selected
for(i_cycle in 1:n_cycle) {
  m_M_SoC[i_cycle + 1, ] <- m_M_SoC[i_cycle, ] %*% m_P_SoC[ , , i_cycle]
  m_M_Exp[i_cycle + 1, ] <- m_M_Exp[i_cycle, ] %*% m_P_Exp[ , , i_cycle]
}

# Plotting the Markov cohort traces
matplot(m_M_SoC, 
        type = "l", 
        ylab = "Probability of state occupancy",
        xlab = "Cycle",
        main = "Makrov Cohort Traces",
        lwd  = 3,
        lty  = 1)
matplot(m_M_Exp, 
        type = "l", 
        lwd  = 3,
        lty  = 3,
        add  = TRUE)
legend("right", 
       legend = c(paste(v_names_states, "(SOC)"), paste(v_names_states, "(Exp)")), 
       col    = rep(c("black", "red", "green"), 2), 
       lty    = c(1, 1, 1, 3, 3, 3), 
       lwd    = 3,
       bty    = "n")

# Calculating and plotting overal survival (OS)
v_OS_SoC <- 1 - m_M_SoC[, "Dead"]
v_OS_Exp <- 1 - m_M_Exp[, "Dead"]

plot(v_OS_SoC, 
     type = "l",
     ylim = c(0, 1),
     ylab = "Survival probability",
     xlab = "Cycle",
     main = "Overall Survival",
     lwd  = 3)
lines(v_OS_Exp,
      lty = 3,
      lwd = 3)
legend("right",
       legend = c("SoC", "Exp"),
       lty    = c(1, 3),
       lwd    = 3,
       bty    = "n")


## 3.5 Health economic outcomes ----

# Calculate the costs and QALYs per cycle by multiplying m_M with the cost/utility vectors for the different states
# - note that to obtain QALYs, the utility needs to be mutiplied by the cycle length as well
v_tc_SoC <- m_M_SoC %*% c(c_F_SoC, c_P, c_D)
v_tc_Exp <- m_M_Exp %*% c(c_F_Exp, c_P, c_D)

v_tu_SoC <- m_M_SoC %*% c(u_F, u_P, u_D) * t_cycle
v_tu_Exp <- m_M_Exp %*% c(u_F, u_P, u_D) * t_cycle

# Obtained the discounted costs and QALYs by multiplying the vectors by the discount rate for each cycle
# - note first the discount rate for each cycle needs to be defined accounting for the cycle length
v_dwc <- 1 / ((1 + d_c) ^ ((0:n_cycle) * t_cycle)) 
v_dwe <- 1 / ((1 + d_e) ^ ((0:n_cycle) * t_cycle))

tc_d_SoC <-  t(v_tc_SoC) %*% v_dwc 
tc_d_Exp <-  t(v_tc_Exp) %*% v_dwc

tu_d_SoC <-  t(v_tu_SoC) %*% v_dwe
tu_d_Exp <-  t(v_tu_Exp) %*% v_dwe

# The discounted costs and QALYs can be summarized and visualized using functions from the 'dampack' package
(df_cea <- calculate_icers(cost       = c(tc_d_SoC, tc_d_Exp),
                           effect     = c(tu_d_SoC, tu_d_Exp),
                           strategies = v_names_strats))

plot(df_cea, effect_units = "QALYs")




### 4. DETERMINISTIC SENSITIVITY ANALYSIS ----

# The 'damrack' package provides convenient function to perform a deterministic sensitivity analysis for the model,
# which is demonstrated in this section


## 4.1 Initialization ----

# Load the model as a function that is defined in the supporting script
source(file = "oncologySemiMarkov_function.R")

# Test whether the function works (and generates the same results)
# - to do so, first a list of parameter values needs to be generated
l_params_all <- list(
  coef_weibull_shape_SoC = 0.7700246,
  coef_weibull_scale_SoC = 1.742534,
  HR_FP_Exp = 0.6,
  p_FD      = 0.02,      
  p_PD      = 0.1,      
  c_F_SoC   = 400,  
  c_F_Exp   = 800,
  c_P       = 1000,  
  c_D       = 0,     
  u_F       = 0.8,   
  u_P       = 0.5,   
  u_D       = 0,  
  d_e       = 0.03,  
  d_c       = 0.03,
  n_cycle   = 60,
  t_cycle   = 0.25
)

oncologySemiMarkov(l_params_all = l_params_all, n_wtp = 20000)


## 4.2 Defining and performing a one-way sensitivity analysis ----

# A one-way sensitivity analysis (OWSA) can be defined by specifying the names of the parameters that are to be
# incuded and their minimum and maximum values
# - here, the cost per cycle of the experimental treatment and utility for the progression-free state are included
df_params_OWSA <- data.frame(
  pars = c("c_F_Exp", "u_F"),   # names of the parameters to be changed
  min  = c(400,  0.70),         # min parameter values
  max  = c(1200, 0.90)          # max parameter values
)

# The OWSA is performed using the run_owsa_det function
OWSA_NMB  <- run_owsa_det(
  params_range     = df_params_OWSA,     # dataframe with parameters for OWSA
  params_basecase  = l_params_all,       # list with all parameters
  nsamp            = 100,                # number of parameter values
  FUN              = oncologySemiMarkov, # function to compute outputs
  outcomes         = c("NMB"),           # output to do the OWSA on
  strategies       = c("SoC", "Exp"),    # names of the strategies
  n_wtp            = 20000               # extra argument to pass to FUN to specify the willingness to pay
)

# Plotting the outcomes of the OWSA in a tornado plot
# - note that other plots can also be generated using the plot() and owsa_opt_strat() functions
owsa_tornado(owsa = OWSA_NMB, txtsize = 11)


## 4.2 Defining and performing a two-way sensitivity analysis ----

# To perform a two-way sensitivity analysis (TWSA), a similar data.frame with model parameters is required
df_params_TWSA <- df_params_OWSA

# The TWSA is performed using the run_twsa_det function
TWSA_NMB <- run_twsa_det(
  params_range    = df_params_TWSA,     # dataframe with parameters for TWSA
  params_basecase = l_params_all,       # list with all parameters
  nsamp           = 40,                 # number of parameter values
  FUN             = oncologySemiMarkov, # function to compute outputs
  outcomes        = c("NMB"),           # output to do the TWSA on
  strategies      = c("SoC", "Exp"),    # names of the strategies
  n_wtp           = 20000               # extra argument to pass to FUN to specify the willingness to pay
)

plot(TWSA_NMB)




### 5. PROBABILISTIC ANALYSIS ----

# A probabilistic analysis (PA) of the model is performed to estimate the uncertainty in the model outcomes. This
# involves generating samples for the model parameters from parametric distributions and evaluating the model for
# each set of parameter samples.

## 5.1 Sampling the parameter values ----

# Defining the number of runs and random number seed (for reproducibility)
n_runs <- 1000
set.seed(123)

# Sampling the parameter values from the distributions and storing them in a data.frame
# - note that some of the model parameters have no uncertainty and, hence, are fixed
# - note that correlated sets of coefficients for the survival distribution are generated using the variance-
#   covariance matrix that was obtained from fitting the weibull distribution and a multivariate normal distribution

m_coef_weibull_SoC <- mvrnorm(
  n     = n_runs, 
  mu    = l_TTP_SoC_weibull$coefficients, 
  Sigma = l_TTP_SoC_weibull$cov
)

head(m_coef_weibull_SoC)

df_PA_input <- data.frame(
  coef_weibull_shape_SoC = m_coef_weibull_SoC[ , "shape"],
  coef_weibull_scale_SoC = m_coef_weibull_SoC[ , "scale"],
  HR_FP_Exp = exp(rnorm(n_runs, log(0.6), 0.08)),
  p_FD      = rbeta(n_runs, shape1 = 16, shape2 = 767),
  p_PD      = rbeta(n_runs, shape1 = 22.4, shape2 = 201.6),
  c_F_SoC   = rgamma(n_runs, shape = 16, scale = 25), 
  c_F_Exp   = rgamma(n_runs, shape = 16, scale = 50), 
  c_P       = rgamma(n_runs, shape = 100, scale = 10), 
  c_D       = 0, 
  u_F       = rbeta(n_runs, shape1 =  50.4, shape2 = 12.6), 
  u_P       = rbeta(n_runs, shape1 = 49.5, shape2 = 49.5), 
  u_D       = 0,
  d_c       = 0.03,
  d_e       = 0.03,
  n_cycle   = 60,
  t_cycle   = 0.25
)

# Inspect the data.frame to see what it looks like
head(df_PA_input)


## 5.2 Running the probabilistic analysis ----

# Initialize data.frames to store the output
df_c <- df_e <- data.frame(
  SoC = rep(NA, n_runs),
  Exp = rep(NA, n_runs)
)

head(df_c)

# Evaluate the Markov model for each set (row) of parameter values
# - note this loop can be run in parallel to decrease the runtime
for(i_run in 1:n_runs){
  
  # Evaluate the model and store the outcomes
  l_out_temp    <- oncologySemiMarkov(l_params_all = df_PA_input[i_run, ], n_wtp = 20000)
  df_c[i_run, ] <- l_out_temp$Cost
  df_e[i_run, ] <- l_out_temp$Effect
  
  # Display simulation progress
  if(i_run/(n_runs/10) == round(i_run/(n_runs/10), 0)) { # display progress every 10%
    cat('\r', paste(i_run/n_runs * 100, "% done", sep = " "))
  }
}


## 5.3 Visualizing the probabilistic analysis results ----

# The 'dampack' package contains multiple useful functions that summarize and visualize the results of a
# probabilitic analysis. To use those functions, the data has to be in a particular structure.
l_PA <- make_psa_obj(cost          = df_c, 
                     effectiveness = df_e, 
                     parameters    = df_PA_input, 
                     strategies    = c("SoC", "Exp"))

# Mean outcome estimates
(df_out_ce_PA <- summary(l_PA))

# Summary of health economic outcomes
(df_cea_PA <- calculate_icers(cost       = df_out_ce_PA$meanCost, 
                              effect     = df_out_ce_PA$meanEffect,
                              strategies = df_out_ce_PA$Strategy))

# Incremental cost-effectiveness plane
plot(l_PA)

# Cost-effectiveness acceptability curve (CEAC)
# - note that a vector of willingness to pay values is first generated to define for which values the CEAC is made
v_wtp <- seq(0, 50000, by = 100)
CEAC_obj <- ceac(wtp = v_wtp, psa = l_PA)

# Regions of highest probability of cost-effectiveness for each strategy
summary(CEAC_obj)

# CEAC and cost-effectiveness acceptability frontier (CEAF) plot
plot(CEAC_obj)

# Expected value of perfect information (EVPI)
EVPI_obj <- calc_evpi(wtp = v_wtp, psa = l_PA)
plot(EVPI_obj, effect_units = "QALY")


