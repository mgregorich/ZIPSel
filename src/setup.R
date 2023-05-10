#' ==============================================================================
#' Author: MG
#' Date: 09/05/2023
#' Info: Setup of parameters, packages, paths, etc.
#' ==============================================================================


# Packages
pacman::p_load(ggplot2, parallel, future.apply, stringr, kableExtra,
               openxlsx, dplyr, tidyverse,tableone,  concreg,
               glmnet, gglasso, oem, WGCNA)

# Paths
sim.date <- Sys.Date()
sim.file <- paste0("sim_", sim.date,"/")
sim.path <- here::here("output", sim.file)

# R files
source(here::here("src","functions.R"))

# Parameter
iter = 1000
n = c(100,200,400)         # sample size
p = c(200,400)             # number of candidate predictors
m = 4                      # number of modules
setting = c("one-module-effect", "two-module-effect", "random-selection")

# Scenario matrix
scenarios <- expand.grid(
  iter = iter,
  n = n,
  p = p,
  m = m,
  setting = setting)

