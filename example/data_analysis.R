# ============================================================================ #
# Author: Mariella Gregorich
# Date: 17/08/2023
# Info: Analysis of the Mosaiques Diagnostics data
# ============================================================================ #


rm(list = ls())

# Packages
pacman::p_load(ggplot2, parallel, future.apply, stringr, 
               dplyr, tidyverse, tableone, concreg, gt, tictoc,
               glmnet, patchwork, ranger, tuneRanger, RColorBrewer, ggvenn)
source(here::here("src","functions_aux.R"))


# ---------- Initialization & Data --------------
ncv = 10
nR  = 2
propnz = 0.9
data_mos <- readRDS(here::here("data", "mosaique", "Large_Mariella.rds"))
data_peptides <- data_mos[,7:ncol(data_mos)]
colnames(data_peptides) <- paste0("P",1:ncol(data_peptides))
rownames(data_peptides) <- paste0("N",1:nrow(data_peptides))
data_peptides <- apply(as.matrix(data_peptides), 2, as.numeric)
clin.vars <- c("Age", "Sex", "eGFR")

# Peptide subset selection due to computational feasibility based on max non-zero entries
ind <- which(apply(data_peptides, 2, assess_nonzeros, prop_nonzero = propnz))
data_peptides <- as.matrix(data_peptides[,ind])
data_mos <- data.frame(cbind(data_mos[,1:6],data_peptides)) %>%
  dplyr::rename("Sex" = Gender)

# Parameter
n = nrow(data_mos)
p = ncol(data_peptides)
folds <- sample(rep(1:ncv, ceiling(n/ncv)))[1:n]
pflist <- list(c(1,2), c(2,1), c(1,3), c(3,1))


# Data preparation
data_mos$Sex <- to_factor(ifelse(data_mos$Sex == "female", 1, 0))
data_mos$Age <- as.numeric(data_mos$Age)
data_mos$eGFR <- as.numeric(data_mos$eGFR)
data_mos$log2eGFR <- log(data_mos$eGFR, 2)
marker.vars <- colnames(data_mos[,!colnames(data_mos) %in% c("fidAuswertung",clin.vars)])

# Data transformation
x <- data_peptides
colnames(x) <- paste0("x.", colnames(x))
d <- (data_peptides != 0)*1
colnames(d) <- paste0("d.", colnames(d))
utmp <- apply(data_peptides, 2, function(x) ifelse(x == 0, 0, log2(x)))
u <- apply(utmp, 2, function(x) ifelse(x == 0, mean(x[x > 0]), x))
colnames(u) <- paste0("u.", colnames(u))
global_min <- min(data_peptides[data_peptides > 0])
ximp <- apply(data_peptides, 2, function(x) log2(ifelse(x == 0, global_min*(1/2), x)))
c <- data.frame("age" = data_mos$Age, "sex" = data_mos$Sex)

# Data object for models
data.obj <- list("y" = data_mos$log2eGFR, "clinical" = c, x = ximp, u = u, d = d,
                 pred.ridge.cb = NA, pred.ridge.cp = NA,
                 pred.lasso = NA,
                 pred.lridge.cb = NA, pred.lridge.cp = NA,
                 pred.rgarrote.cb = NA, pred.rgarrote.cp = NA,
                 pred.rf.notune = NA, pred.rf.tune = NA)


# Data frame for random forest
dataset <- data.frame(y = log2(data_mos$eGFR), age = data_mos$Age, sex = data_mos$Sex, x, pred.rf.notune = NA, pred.rf.tune = NA)

# ----------- CV Predictions --------------------

names_models <- rep(c("lasso", "ridge", "lasso-ridge", "ridge-lasso", "ridge-garrote", "random forest"), each=2)[-2]
names_penalty <- c(rep(c("combined", "component"), 5)[-2], "without tuning", "tuning")

perform_CVmethods <- function(i, data.obj){
  print(paste0("NCV iteration = ", i,"/", ncv))
  
  tbl_out <- cbind.data.frame(names_models, names_penalty, matrix(NA, nrow=length(names_models), ncol=4))
  colnames(tbl_out) <- c("model", "penalty", "R2", "RMSPE", "MAE", "CS")
  
  obj.train <- list("y" = data.obj$y[(1:n)[folds != i]],
                    "clinical" = data.obj$clinical[(1:n)[folds != i], ], 
                    "x" = data.obj$x[(1:n)[folds != i], ], 
                    "u" = data.obj$u[(1:n)[folds != i], ],
                    "d" = data.obj$d[(1:n)[folds != i], ])
  obj.test <- list("y" = data.obj$y[(1:n)[folds == i]],
                   "clinical" = data.obj$clinical[(1:n)[folds == i],], 
                   "x" = data.obj$x[(1:n)[folds == i], ],
                   "u" = data.obj$u[(1:n)[folds == i], ], 
                   "d" = data.obj$d[(1:n)[folds == i], ])
  train <- dataset[(1:n)[folds != i], ]
  test <- dataset[(1:n)[folds == i], ]
  
  # Models
  fit.lasso.cb <- perform_penreg(obj.train, family = "gaussian", alpha1 = 1, cv = ncv, R = nR, penalty = "combined")
  fit.ridge.cb <- perform_penreg(obj.train, family="gaussian", alpha1 = 0, cv=ncv, R=nR, penalty = "combined")
  fit.ridge.cp <- perform_penreg(obj.train, family = "gaussian", alpha1 = 0, cv = ncv, R = nR, penalty  =  "component", pflist = pflist)
  fit.lridge.cb <- perform_lridge(obj.train, family = "gaussian", cv = ncv, R = nR, penalty = "combined")
  fit.lridge.cp <- perform_lridge(obj.train, family = "gaussian", cv = ncv, R = nR, penalty = "component", pflist = pflist)
  fit.rlasso.cb <- perform_rlasso(obj.train, family = "gaussian", cv = ncv, R = nR, penalty = "combined")
  fit.rlasso.cp <- perform_rlasso(obj.train, family = "gaussian", cv = ncv, R = nR, penalty = "component", pflist = pflist)
  fit.rgarrote.cb <- perform_rgarrote(obj.train, family = "gaussian", cv = ncv, R = nR, penalty = "combined")
  fit.rgarrote.cp <- perform_rgarrote(obj.train, family = "gaussian", cv = ncv, R = nR, penalty = "component", pflist = pflist)
  fit.rf.notune <- ranger(y~., data = train %>% dplyr::select(-c(pred.rf.notune, pred.rf.tune)), num.trees = 1000)
  fit.rf.tune <- tuneMtryFast(y~., data = train %>% dplyr::select(-c(pred.rf.notune, pred.rf.tune)), num.trees = 1000, stepFactor = 0.5,
                              doBest = TRUE, plot = FALSE, trace = FALSE)
  
  
  # Predictions
  data.obj$pred.lasso.cb[(1:n)[folds == i]] <- predict_penreg(obj = fit.lasso.cb, newdata = obj.test, model = "lasso")
  data.obj$pred.ridge.cb[(1:n)[folds == i]] <- predict_penreg(obj = fit.ridge.cb,newdata = obj.test,model = "ridge")
  data.obj$pred.ridge.cp[(1:n)[folds == i]] <- predict_penreg(obj = fit.ridge.cp, newdata = obj.test, model = "ridge")
  data.obj$pred.lridge.cb[(1:n)[folds == i]] <- predict_lridge(obj = fit.lridge.cb, newdata = obj.test)
  data.obj$pred.lridge.cp[(1:n)[folds == i]] <- predict_lridge(obj = fit.lridge.cp, newdata = obj.test)
  data.obj$pred.rlasso.cb[(1:n)[folds == i]] <- predict_rlasso(obj = fit.rlasso.cb, newdata = obj.test)
  data.obj$pred.rlasso.cp[(1:n)[folds == i]] <- predict_rlasso(obj = fit.rlasso.cp, newdata = obj.test)
  data.obj$pred.rgarrote.cb[(1:n)[folds == i]] <- predict_rgarrote(obj = fit.rgarrote.cb, newdata = obj.test)
  data.obj$pred.rgarrote.cp[(1:n)[folds == i]] <- predict_rgarrote(obj = fit.rgarrote.cp, newdata = obj.test)
  dataset$pred.rf.notune[(1:n)[folds == i]] <- predict(fit.rf.notune,  data = test %>% dplyr::select(-c(y, pred.rf.notune, pred.rf.tune)))$predictions
  dataset$pred.rf.tune[(1:n)[folds == i]] <- predict(fit.rf.tune,  data = test %>% dplyr::select(-c(y, pred.rf.tune, pred.rf.notune)))$predictions
  
  # Performance
  tbl_out[1, 3:6] <- eval_performance(data.obj$pred.lasso.cb[(1:n)[folds == i]], obj.test$y)
  tbl_out[2, 3:6] <- eval_performance(data.obj$pred.ridge.cb[(1:n)[folds == i]], obj.test$y)
  tbl_out[3, 3:6] <- eval_performance(data.obj$pred.ridge.cp[(1:n)[folds == i]], obj.test$y)
  tbl_out[4, 3:6] <- eval_performance(data.obj$pred.ridge.cb[(1:n)[folds == i]], obj.test$y)
  tbl_out[5, 3:6] <- eval_performance(data.obj$pred.lridge.cp[(1:n)[folds == i]], obj.test$y)
  tbl_out[6, 3:6] <- eval_performance(data.obj$pred.rlasso.cb[(1:n)[folds == i]], obj.test$y)
  tbl_out[7, 3:6] <- eval_performance(data.obj$pred.rlasso.cp[(1:n)[folds == i]], obj.test$y)
  tbl_out[8, 3:6] <- eval_performance(data.obj$pred.rgarrote.cb[(1:n)[folds == i]], obj.test$y)
  tbl_out[9, 3:6] <- eval_performance(data.obj$pred.rgarrote.cp[(1:n)[folds == i]], obj.test$y)
  tbl_out[10, 3:6] <- eval_performance(pred = dataset$pred.rf.notune[(1:n)[folds == i]], obs = test$y)
  tbl_out[11, 3:6] <- eval_performance(pred = dataset$pred.rf.tune[(1:n)[folds == i]], obs = test$y)
  return(tbl_out)
}

plan(multisession, workers=10)
tbl_performance <- future_lapply(1:ncv, function(x) perform_CVmethods(x, data.obj = data.obj), future.seed=TRUE)
plan(sequential)
tbl_performance <- do.call(rbind, tbl_performance)

# -------------------  Final model ----------------------
list_model <- list()
list_model$fit.lasso.cb <- perform_penreg(data.obj = data.obj, family = "gaussian", cv = ncv, R = nR, alpha1 = 1, penalty = "combined")
list_model$fit.ridge.cb <- perform_penreg(data.obj = data.obj, family = "gaussian", cv = ncv, R = nR, penalty = "combined")
list_model$fit.ridge.cp <- perform_penreg(data.obj = data.obj, family = "gaussian", cv = ncv, R = nR, penalty = "component", pflist = pflist)
list_model$fit.lridge.cb <- perform_lridge(data.obj = data.obj, family = "gaussian", cv = ncv, R = nR, penalty = "combined")
list_model$fit.lridge.cp <- perform_lridge(data.obj = data.obj, family = "gaussian", cv = ncv, R = nR, penalty = "component", pflist = pflist)
list_model$fit.rlasso.cb <- perform_rlasso(data.obj = data.obj, family = "gaussian", cv = ncv, R = nR, penalty = "combined")
list_model$fit.rlasso.cp <- perform_rlasso(data.obj = data.obj, family = "gaussian", cv = ncv, R = nR, penalty = "component", pflist = pflist)
list_model$fit.rgarrote.cb <- perform_rgarrote(data.obj = data.obj, family = "gaussian", cv = ncv, R = nR, penalty = "combined")
list_model$fit.rgarrote.cp <- perform_rgarrote(data.obj = data.obj, family = "gaussian", cv = ncv, R = nR, penalty = "component", pflist = pflist)


# Coefficients
list_coef <- list()
list_coef$lasso.cb <- list_model$fit.lasso.cb$fit$coefficients
list_coef$ridge.cb <- list_model$fit.ridge.cb$fit$coefficients
list_coef$ridge.cp <- list_model$fit.ridge.cp$fit$coefficients
list_coef$lridge.cb <- list_model$fit.lridge.cb$fit$coefficients
list_coef$lridge.cp <- list_model$fit.lridge.cp$fit$coefficients
list_coef$rlasso.cb <- list_model$fit.rlasso.cb$fit$coefficients
list_coef$rlasso.cp <- list_model$fit.rlasso.cp$fit$coefficients
list_coef$rgarrote.cb <- list_model$fit.rgarrote.cb$fit$coefficients
list_coef$rgarrote.cp <- list_model$fit.rgarrote.cp$fit$coefficients

# Variable selection
extract_varsel <- function(modelcoef){
  model.varsel <- rownames(modelcoef)[which(modelcoef != 0)][-1]
  model.varsel <- str_remove_all(model.varsel, "u.|d.|x.")
  return(model.varsel[!duplicated(model.varsel)])
}

# List of variable selection for models 
list_varsel <- list()
list_varsel$all <- colnames(data_peptides)
list_varsel$lasso.cb <- extract_varsel(list_coef$lasso.cb)
list_varsel$ridge.cb <- extract_varsel(list_coef$ridge.cb)
list_varsel$ridge.cp <- extract_varsel(list_coef$ridge.cp)
list_varsel$lridge.cb <- extract_varsel(list_coef$lridge.cb)
list_varsel$lridge.cp <- extract_varsel(list_coef$lridge.cp)
list_varsel$rlasso.cb <- extract_varsel(list_coef$rlasso.cb)
list_varsel$rlasso.cp <- extract_varsel(list_coef$rlasso.cp)
list_varsel$rgarrote.cb <- extract_varsel(list_coef$rgarrote.cb)
list_varsel$rgarrote.cp <- extract_varsel(list_coef$rgarrote.cp)


tbl_perf <- tbl_performance %>%
  data.frame() %>%
  mutate(model=fct_relevel(model, names_models[!duplicated(names_models)])) %>%
  group_by(model, penalty) %>%
  summarise("R2"=mean(R2, na.rm=TRUE), "RMSPE"=mean(RMSPE, na.rm=TRUE), 
            "MAE"=mean(MAE, na.rm=TRUE), "CS"=mean(CS, na.rm=TRUE)) %>%
  ungroup() %>%
  mutate("df"=c(unlist(lapply(list_model, function(x) x$dfvars)), rep(ncol(data_peptides),2)),
         "propnz" = propnz) %>%
  relocate(propnz, .after=penalty)

# Runtime
library(microbenchmark)
mbm <- microbenchmark("lasso" = { b <- perform_penreg(data.obj = data.obj, family = "gaussian", cv = ncv, R = nR, alpha1 = 1, penalty = "combined")},
                      "ridge_cb" = { b <- perform_penreg(data.obj = data.obj, family = "gaussian", cv = ncv, R = nR, penalty = "combined")},
                      "ridge_cp" = { b <- perform_penreg(data.obj = data.obj, family = "gaussian", cv = ncv, R = nR, penalty = "component", pflist = pflist)},
                      "lridge_cb" = { b <- perform_lridge(data.obj = data.obj, family = "gaussian", cv = ncv, R = nR, penalty = "combined")},
                      "lridge_cp" = { b <- perform_lridge(data.obj = data.obj, family = "gaussian", cv = ncv, R = nR, penalty = "component", pflist = pflist)},
                      "rlasso_cb" = { b <- perform_rlasso(data.obj = data.obj, family = "gaussian", cv = ncv, R = nR, penalty = "combined")},
                      "rlasso_cp" = { b <- perform_rlasso(data.obj = data.obj, family = "gaussian", cv = ncv, R = nR, penalty = "component", pflist = pflist)},
                      "rgarrote_cb" = { b <- perform_rgarrote(data.obj = data.obj, family = "gaussian", cv = ncv, R = nR, penalty = "combined")},
                      "rgarrote_cp" = { b <- perform_rgarrote(data.obj = data.obj, family = "gaussian", cv = ncv, R = nR, penalty = "component", pflist = pflist)},
                      "rf_notune" = { b <- ranger(y~., data = dataset %>% dplyr::select(-c(pred.rf.notune, pred.rf.tune)), num.trees = 1000)},
                      "rf_tune" =  { b <- tuneMtryFast(y~., data = dataset %>% dplyr::select(-c(pred.rf.notune, pred.rf.tune)), num.trees = 1000, 
                                                       stepFactor = 0.5, doBest = TRUE, plot = FALSE, trace = FALSE)}, times = 1
                      )


# ------------ Save results ---------------
list_results <- list("performance" = tbl_perf, 
                     "performance_ncv" = tbl_performance, 
                     "coef" = list_coef,
                     "models" = list_model,
                     "varsel" = list_varsel,
                     "mbm" = mbm)
filename <- paste0("results_mosaique_pnz", propnz, "_nR", nR,"_ncv", ncv,".rds")
saveRDS(list_results, here::here("output", "example", filename))
