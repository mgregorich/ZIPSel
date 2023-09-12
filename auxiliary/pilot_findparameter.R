#' ==============================================================================
#' Author: MG
#' Date: 09/05/2023
#' Info: Pilot study to find a and residual variance
#' ==============================================================================

rm(list=ls())

# Packages
pacman::p_load(ggplot2, parallel, future.apply, stringr, kableExtra,
               openxlsx, dplyr, tidyverse, tableone, concreg, Matrix,
               glmnet, MASS, ranger, simdata)
source(here::here("src","functions_aux.R"))
source(here::here("src","functions_sim.R"))

# Load & save setup
# source("src/setup.R")
setup <- readRDS(here::here("src", "scenario_setup.rds"))
scenarios <- setup[[1]]
sim_design <- setup[[2]]


# --- Determine optimal residual variance to achieve R2: 0.2 and 0.5
scenarios <- scenarios %>% 
  mutate(epsstd = ifelse(scenario == "A" & p == 200 & a == .2165, epsstd * 5.5, epsstd),
         epsstd = ifelse(scenario == "A" & p == 200 & a == .28, epsstd * 6.5, epsstd),
         epsstd = ifelse(scenario == "A" & p == 200 & a == 1, epsstd * 20.25, epsstd),
         epsstd = ifelse(scenario == "A" & p == 400 & a == .2165, epsstd * 10, epsstd),
         epsstd = ifelse(scenario == "A" & p == 400 & a == .28, epsstd * 12, epsstd),
         epsstd = ifelse(scenario == "A" & p == 400 & a == 1, epsstd * 40, epsstd)) %>% 
  arrange(p, a, epsstd) %>%
  filter(scenario %in% "A")
scenarios$R2 <- NA
for(i in 1: nrow(scenarios)){
  print(i)
  scn = scenarios[i,]
  dsgn = sim_design[[scenarios[i,]$dsgn]]
  data.val <- data_generation(dsgn = dsgn, scenario = scn$scenario, n = 100000, p=scn$p, beta_max = scn$beta_max, a = scn$a, epsstd = scn$epsstd, 
                              prop.nonzero = scn$prop.nonzero, sampthresh = scn$sampthresh)
  
  # R2
  y_noisy <- data.val$data_gen$y
  y_true <- y_noisy - data.val$data_gen$eps
  SST <- sum((y_noisy - mean(y_noisy))^2)
  SSE <- sum((y_true-y_noisy)^2)
  R2 <- 1- (SSE/SST)
  scenarios[i,]$R2 <- R2
}

