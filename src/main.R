# ==============================================================================
# Author: MG
# Date: 09/05/2023
# Info: Main execution file for the simulation study 
# ==============================================================================

rm(list=ls())

# Load & save setup
source("src/setup.R")

# Create sim folder
if(dir.exists(sim.path)){invisible(do.call(file.remove, list(list.files(sim.path, full.names = T))))
}else{dir.create(sim.path)}

setup <- readLines("src/setup.R")
write.table(setup, file = here::here(sim.path, "info_setup.txt"))

# ======================= Simulation ===========================================

# --- Run through all scenarios
# plan(multisession, workers = detectCores()*.5)
# invisible(future_lapply(1:nrow(scenarios), function(k) {
#   tryCatch({simulate_scenario(scn=scenarios[k,])
#   }, error = function(e) {
#     # log the error message to a file
#     cat(paste0("Error in scenario k=",k,": ", e$message, "\n"), file = paste0(sim.path, "/error_log.txt"), append = TRUE)
#   })
# }, future.seed = T))
# plan(sequential)
