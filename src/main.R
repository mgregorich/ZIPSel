# ==============================================================================
# Author: MG
# Date: 09/05/2023
# Info: Main execution file for the simulation study 
# ==============================================================================

rm(list=ls())

# Packages
pacman::p_load(ggplot2, parallel, future.apply, stringr, kableExtra,
               openxlsx, dplyr, tidyverse, tableone, concreg, Matrix,
               glmnet, MASS, ranger, simdata)

# Paths
sim.date <- Sys.Date()
sim.file <- paste0("sim_", sim.date,"/")
sim.path <- here::here("output", sim.file)

# Create output folder
if(dir.exists(sim.path)){invisible(do.call(file.remove, list(list.files(sim.path, full.names = T))))
}else{dir.create(sim.path)}

# R files
source(here::here("src","functions_aux.R"))
source(here::here("src","functions_sim.R"))

# Load & save setup
# source("src/setup.R")
setup <- readRDS(here::here("src", "scenario_setup.rds"))
scenarios <- setup[[1]]
sim_design <- setup[[2]]
scenarios$iter <- 2

# ======================= Simulation ===========================================

# --- Run through all scenarios
plan(multisession, workers = detectCores()*.75)
invisible(future_lapply(1:nrow(scenarios), function(k) {
  tryCatch({simulate_scenario(scn=scenarios[k,], dsgn=sim_design[[scenarios[k,]$dsgn]])
  }, error = function(e) {
    # log the error message to a file
    cat(paste0("Error in scenario k=",k,": ", e$message, "\n"), file = paste0(sim.path, "/error_log.txt"), append = TRUE)
  })
}, future.seed = T))
plan(sequential)


# --- Summarize all scenarios
sim.files <- list.files(sim.path, pattern = "sim_")
sim.all <- evaluate_scenarios(sim.files, sim.path)

