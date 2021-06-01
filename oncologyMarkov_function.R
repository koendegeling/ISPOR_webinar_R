##----------------------------------------------------------------------------------------------------------------##
##                                                                                                                ##
##                      ISPOR STUDENT NETWORK WEBINAR: HEALTH ECONOMIC MODELLING IN R                             ##
##                                                                                                                ##
##                               by Petros Pechlivanoglou and Koen Degeling                                       ##
##                                                                                                                ##
##----------------------------------------------------------------------------------------------------------------##
#
# This script is used in the oncologyMarkov_illustration.R script to define the function that evaluated the three- 
# state Markov model for a specific set of parameter values.

oncologyMarkov <- function(l_params_all, n_wtp = 10000) {
  
  with(as.list(l_params_all), {
    
    v_names_strats <- c("Standard of Care", "Experimental")         # strategy names
    v_names_states <- c("ProgressionFree", "Progression", "Dead")   # state names
    n_strats       <- length(v_names_strats)                        # number of strategies
    n_states       <- length(v_names_states)                        # number of states
    
    m_M_SoC <- m_M_Exp <- matrix(
      data = NA, 
      nrow = n_cycle,  
      ncol = n_states, 
      dimnames = list(paste('Cycle', 1:n_cycle), v_names_states)
    )
    
    m_M_SoC[1, "ProgressionFree"] <- m_M_Exp[1, "ProgressionFree"] <- 1
    m_M_SoC[1, "Progression"]     <- m_M_Exp[1, "Progression"]     <- 0
    m_M_SoC[1, "Dead"]            <- m_M_Exp[1, "Dead"]            <- 0
    
    m_P_SoC <- matrix(
      data = 0,
      nrow = n_states, 
      ncol = n_states,
      dimnames = list(v_names_states, v_names_states)
    )
    
    m_P_SoC["ProgressionFree", "ProgressionFree"] <- (1 - p_FD) * (1 - p_FP_SoC)
    m_P_SoC["ProgressionFree", "Progression"]     <- (1 - p_FD) * p_FP_SoC
    m_P_SoC["ProgressionFree", "Dead"]            <- p_FD
    
    m_P_SoC["Progression", "Progression"] <- 1 - p_PD
    m_P_SoC["Progression", "Dead"]        <- p_PD
    
    m_P_SoC["Dead", "Dead"] <- 1
    
    m_P_Exp <- m_P_SoC
    m_P_Exp["ProgressionFree", "ProgressionFree"] <- (1 - p_FD) * (1 - p_FP_Exp)
    m_P_Exp["ProgressionFree", "Progression"]     <- (1 - p_FD) * p_FP_Exp
    
    check_transition_probability(m_P_SoC, verbose = TRUE)
    check_transition_probability(m_P_Exp, verbose = TRUE)
    
    check_sum_of_transition_array(m_P_SoC, n_states = n_states, n_cycles = n_cycle, verbose = TRUE)
    check_sum_of_transition_array(m_P_Exp, n_states = n_states, n_cycles = n_cycle, verbose = TRUE)
    
    for(i_cycle in 1:(n_cycle-1)) {
      m_M_SoC[i_cycle + 1, ] <- m_M_SoC[i_cycle, ] %*% m_P_SoC
      m_M_Exp[i_cycle + 1, ] <- m_M_Exp[i_cycle, ] %*% m_P_Exp
    }
    
    v_tc_SoC <- m_M_SoC %*% c(c_F_SoC, c_P, c_D)
    v_tc_Exp <- m_M_Exp %*% c(c_F_Exp, c_P, c_D)
    
    v_tu_SoC <- m_M_SoC %*% c(u_F, u_P, u_D) * t_cycle
    v_tu_Exp <- m_M_Exp %*% c(u_F, u_P, u_D) * t_cycle
    
    v_dwc <- 1 / ((1 + d_c) ^ ((0:(n_cycle-1)) * t_cycle)) 
    v_dwe <- 1 / ((1 + d_e) ^ ((0:(n_cycle-1)) * t_cycle))
    
    tc_d_SoC <-  t(v_tc_SoC) %*% v_dwc 
    tc_d_Exp <-  t(v_tc_Exp) %*% v_dwc
    
    tu_d_SoC <-  t(v_tu_SoC) %*% v_dwe
    tu_d_Exp <-  t(v_tu_Exp) %*% v_dwe
    
    v_tc_d    <- c(tc_d_SoC, tc_d_Exp)
    v_tu_d    <- c(tu_d_SoC, tu_d_Exp)
    
    v_nmb_d   <- v_tu_d * n_wtp - v_tc_d
    
    df_ce <- data.frame(Strategy = v_names_strats,
                        Cost     = v_tc_d,
                        Effect   = v_tu_d,
                        NMB      = v_nmb_d)
    
    return(df_ce)
    
  })
  
}
