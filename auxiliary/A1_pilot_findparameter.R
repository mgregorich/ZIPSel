#' ==============================================================================
#' Author: MG
#' Date: 09/05/2023
#' Info: Pilot study to find a and residual variance
#' ==============================================================================

rm(list=ls())

# Packages
# devtools::install_github("matherealize/simdata", ref="fix_missing_colapply")
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


# --- (1) Determine regression coefficients ----
scenarios <- scenarios %>% 
  filter(epslvl %in% "none" & UDdepen == "U=D" & struczero==0.66 & propzi==0.75 & revzi==TRUE) %>%
  mutate(n = 10000) %>%
  filter(!duplicated(.)) %>%
  arrange(a,  revzi, struczero, scenario, p, propzi) 

 scenarios$varXb <- NA
for(i in 1:nrow(scenarios)){
  print(i)
  scn = scenarios[i,]
  dsgn = sim_design[[scenarios[i,]$dsgn]]
  data.val <- data_generation(dsgn = dsgn, scenario = scn$scenario, n = scn$n, p=scn$p, beta_max = scn$beta_max, a = scn$a, epsstd = scn$epsstd, 
                              propzi = scn$propzi, revzi = scn$revzi, struczero = scn$struczero)
  
  Xb <- data.val$data_gen$X_true %*% data.val$true_coef$beta_X
  scenarios[i,]$varXb <- as.numeric(var(Xb))
}


# --- (2) Determine parameter a ----
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
  scenarios[i,]$varXb <- var(Xb)
  scenarios[i,]$varDb <- var(Db)
  b <- data.val$true_coef$beta_D
  
  a_U2D <- (4 - sqrt(8*(var(Xb) / var(Db)))) / (2*(2 - var(Xb) / var(Db)))
  a_UD <- (2 - sqrt(4*(var(Xb) / var(Db)))) / (2*(1 - var(Xb) / var(Db)))
  scenarios[i,]$newa <- ifelse(scn$UDdepen %in% "U=2D", a_U2D, a_UD)
  scenarios[i,]$var_Xb_Db <- var(scenarios[i,]$newa*Xb)/var((1-scenarios[i,]$newa)*Db)
}




# --- (3) Determine optimal residual variance to achieve R2: 0.2 and 0.5 ----
scenarios <- scenarios %>% 
  filter(epslvl != "none"  & n == 100) %>%
  arrange(p, a, epsstd) 
scenarios$new_epsstd <- scenarios$checkR2 <- NA
for(i in 1:10){
  print(i)
  scn = scenarios[i,]
  dsgn = sim_design[[scenarios[i,]$dsgn]]
  data.val <- data_generation(dsgn = dsgn, scenario = scn$scenario, n = 100000, p=scn$p, beta_max = scn$beta_max, 
                              a = scn$a, epsstd = scn$epsstd, 
                              propzi = scn$propzi, revzi = scn$revzi, struczero = scn$struczero)
  
  y_true <- data.val$data_gen$y-data.val$data_gen$eps
  r2 <- ifelse(scn$epslvl == "moderate", 0.5, 0.2)
  var.eps <- var(y_true) *((1 - r2)/r2)
  sd.eps <- sqrt(var.eps)
  scenarios[i,]$new_epsstd <- sd.eps
  
  # Check R2
  data.val <- data_generation(dsgn = dsgn, scenario = scn$scenario, n = 10000, p=scn$p, beta_max = scn$beta_max, 
                              a = scn$a, epsstd = sd.eps, 
                            propzi = scn$propzi, revzi = scn$revzi, struczero = scn$struczero)
  checkR2 <- 1-(var(data.val$data_gen$eps)/var(data.val$data_gen$y))
  print(paste0("True R2 = ", r2, " and est R2 = ", round(checkR2,3)))
  scenarios[i,]$checkR2 <- checkR2
}



