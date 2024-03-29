# ============================================================================ #
# Author: MG
# Date: 09/05/2023
# Info: Main execution file for the simulation study 
# ============================================================================ #

rm(list=ls())
set.seed(666)

# Packages
# remotes::install_github("matherealize/simdata")
pacman::p_load(parallel, future.apply, stringr, dplyr, tidyr, Matrix, 
               glmnet, MASS, simdata, here, simdata)

# Paths
sim.date <- "final"
sim.file <- paste0("sim_", sim.date,"/")
sim.path <- here::here("output", sim.file)


# Create output folder
if(dir.exists(sim.path)){
 invisible(do.call(file.remove, list(list.files(sim.path, full.names = T))))
}else{dir.create(sim.path)}

# R files
source(here::here("src","functions_aux.R"))
source(here::here("src","functions_sim.R"))

# Load & save setup
# source("src/setup.R")
setup <- readRDS(here::here("src", "scenario_setup.rds"))
scenarios <- setup[[1]]
sim_design <- setup[[2]]
scenarios$iter <- 500
nlambdas = 50
nrep = 10


# ======================= Simulation ===========================================

# --- Run through all scenarios
plan(multisession, workers = 6)
# Generate seeds for each simulation replicate
if(dir.exists(here::here(sim.path, "random_seeds.rds"))){
  seeds <- readRDS(here::here(sim.path, "random_seeds.rds"))
}else{
  seeds <- future_lapply(1:scenarios$iter[1], FUN = function(x) .Random.seed, future.chunk.size = Inf, future.seed = 666L)
  saveRDS(seeds, file=here::here(sim.path, "random_seeds.rds"))
}

# Run scenarios
for(k in 1:nrow(scenarios)){
  print(paste0(k, "/", nrow(scenarios), " scenarios"))
  tryCatch({
    simulate_scenario(scn=scenarios[k,], dsgn=sim_design[[scenarios[k,]$dsgn]])
  }, error = function(e) {
    # log the error message to a file
    cat(paste0("Error in scenario k=",k,": ", e$message, "\n"), file = paste0(sim.path, "/error_log.txt"), append = TRUE)
  })
}
plan(sequential)


# --- Summarize all scenarios
sim.files <- list.files(sim.path, pattern = "sim_")
sim.all <- evaluate_scenarios(sim.files, sim.path)

