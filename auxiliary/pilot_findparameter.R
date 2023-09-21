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
#source("src/setup.R")
setup <- readRDS(here::here("src", "scenario_setup.rds"))
scenarios <- setup[[1]]
sim_design <- setup[[2]]

# --- Determine parameter a ---
scenarios <- scenarios %>% 
  filter(epslvl %in% "none" & a!=1) %>%
  mutate(n = 1000) %>%
  filter(!duplicated(.)) %>%
  arrange(a,  revzi, struczero, scenario, p, propzi)
scenarios$var_Xb_Db <- scenarios$varXb <-scenarios$varDb <- scenarios$newa <- NA
for(i in 1:nrow(scenarios)){
  print(i)
  scn = scenarios[i,]
  dsgn = sim_design[[scenarios[i,]$dsgn]]
  data.val <- data_generation(dsgn = dsgn, scenario = scn$scenario, n = scn$n, p=scn$p, beta_max = scn$beta_max, a = scn$a, epsstd = scn$epsstd, 
                              propzi = scn$propzi, revzi = scn$revzi, struczero = scn$struczero)
  
  Xb <- data.val$data_gen$X_true %*% data.val$true_coef$beta_X
  Db <- data.val$data_gen$D_struc %*% data.val$true_coef$beta_D
  scenarios[i,]$varXb <- var(scn$a*Xb)
  scenarios[i,]$varDb <- var((1-scn$a)*Db)
  b <- data.val$true_coef$beta_D
  
  a_U2D <- (4 - sqrt(8*(var(Xb) / var(Db)))) / (2*(2 - var(Xb) / var(Db)))
  a_UD <- (2 - sqrt(4*(var(Xb) / var(Db)))) / (2*(1 - var(Xb) / var(Db)))
  scenarios[i,]$newa <- ifelse(scn$UDdepen %in% "U=2D", a_U2D, a_UD)
  scenarios[i,]$var_Xb_Db <- var(scenarios[i,]$newa*Xb)/var((1-scenarios[i,]$newa)*Db)

}




# --- Determine optimal residual variance to achieve R2: 0.2 and 0.5
scenarios <- scenarios %>% 
  filter(epslvl != "none" & revzi == TRUE) %>%
  arrange(p, a, epsstd) 
scenarios$sd_eps <- scenarios$checkR2 <- NA
for(i in 1:10){
  print(i)
  scn = scenarios[i,]
  dsgn = sim_design[[scenarios[i,]$dsgn]]
  data.val <- data_generation(dsgn = dsgn, scenario = scn$scenario, n = 100000, p=scn$p, beta_max = scn$beta_max, 
                              a = scn$a, epsstd = scn$epsstd, 
                              propzi = scn$propzi, revzi = scn$revzi, struczero = scn$struczero)
  
  # R2
  y_noisy <- data.val$data_gen$y
  y_true <- y_noisy - data.val$data_gen$eps
  # SST <- sum((y_noisy - mean(y_noisy))^2)
  # SSE <- sum((y_true-y_noisy)^2)
  # R2 <- 1- (SSE/SST)
  R2 <- 1-(var(data.val$data_gen$eps)/var(y_noisy))
  
  R2 <- ifelse(scn$epslvl == "moderate", 0.5, 0.2)
  sd_eps <- sqrt((1-R2)* var(y_noisy))
  scenarios[i,]$epsstd <- sqrt(sd_eps)
  
  # Check r2
  data.val <- data_generation(dsgn = dsgn, scenario = scn$scenario, n = 100000, p=scn$p, beta_max = scn$beta_max, 
                              a = scn$a, epsstd = sd_eps, 
                            propzi = scn$propzi, revzi = scn$revzi, struczero = scn$struczero)
  checkR2 <- 1-(var(data.val$data_gen$eps)/var(data.val$data_gen$y))
  
  scenarios[i,]$checkR2 <- cor(y_noisy, y_true)^2
}



