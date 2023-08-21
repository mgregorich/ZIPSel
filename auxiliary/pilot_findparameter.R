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

scenarios$iter <- 5
scenarios$epsstd <- ifelse(scenarios$epsstd==3, 5, 10)

# --- Data generation
nr=3; scn=scenarios[nr,]; dsgn=sim_design[[scenarios[nr,]$dsgn]]

data.val <- data_generation(dsgn = dsgn, n = 100000, p=scn$p, beta_max = scn$beta_max, a = scn$a, epsstd = scn$epsstd, 
                            prop.nonzero = scn$prop.nonzero, sampthresh = scn$sampthresh)
data.val <- generate_dataobj(y = data.val$data_ana$y, x = data.val$data_ana$x, clinical = NULL)


# --- Oracle
X <- cbind(data.obj$u, data.obj$d)
beta_true <- c(df$true_coef$beta_X, df$true_coef$beta_D)
true_pred <- X[,which(beta_true!=0)]
fit.oracle <- lm(df$data_ana$y~true_pred)
beta_oracle <- rep(0, ncol(X))
beta_oracle[which(beta_true!=0)] <- coef(fit.oracle)[-1]
beta_oracle[is.na(beta_oracle)] <- 0 
data.val$pred.oracle <-  fit.oracle$coefficients[1] + c(cbind(data.val$u, data.val$d) %*% beta_oracle) # intercept!!!!
tbl_perf[1, 3:7] <- eval_performance(pred = data.val$pred.oracle, obs = data.val$y)
tbl_coef$beta_oracle_u <- beta_oracle[1:(p)]
tbl_coef$beta_oracle_d <- beta_oracle[(p+1):(2*p)]
