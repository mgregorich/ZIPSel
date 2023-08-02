#' ==============================================================================
#' Author: MG
#' Date: 09/05/2023
#' Info: Functions only relevant for the Monte Carlo simualtion
#' ==============================================================================


# =========================== Run scenario =====================================

simulate_scenario <- function(scen){
  
  # Run replications
  res_all <-  lapply(1:scen$iter, function(x){
    data_iter <- data_generation(n=scen$n, p=scen$p)
    res_iter <-  data_analysis(df=data_iter)
    data_iter$i <- res_iter$i <- x
    return(list("data_gen"=data_iter, "data_ana"=res_iter)) 
  })
  
  data_gen_all <- do.call(rbind, lapply(res_all, function(x) x[[1]]))
  data_ana_all <- do.call(rbind, lapply(res_all, function(x) x[[2]]))
  
  return(list("data_gen_all"=data_gen_all, "data_ana_all"=data_ana_all))
  
}


# ============================ Data generation =================================
data_generation <- function(n, p){
  
  # Step 1: Generate Spike at Zero Variables
  # Remove later (parameter input for function)
  n <- 200     # Number of observations
  p <- 100     # Number of variables
  


  return(out)
}




# ============================ Data analysis =================================
data_analysis <- function(df){
  
  # Prepare data
  
  
  # ridge-garrote
  
  
  # ridge-lasso
  
  
  # lasso-ridge
  
  
  # network-based group penalization
  
  
  # multi univariable
  
  
  # Merge results
  
  return()
}