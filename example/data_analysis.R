# ============================================================================ #
# Author: Mariella Gregorich
# Date: 17/08/2023
# Info: Analysis of the Mosaiques Diagnostics data
# ============================================================================ #


rm(list = ls())

# Packages
pacman::p_load(ggplot2, parallel, future.apply, stringr, 
               dplyr, tidyverse, tableone, concreg, gt, microbenchmark,
               glmnet, patchwork, ranger, tuneRanger, RColorBrewer)
source(here::here("src","functions_aux.R"))
options(dplyr.summarise.inform = FALSE)


# --------- Functions --------------------

perform_CVmethods <- function(i, data.obj, dataset){
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
  fit.lasso <- perform_penreg(obj.train, family = "gaussian", alpha1 = 1, cv = ncv, R = nR, penalty = "combined", split_vars = FALSE)
  fit.ridge.ud <- perform_penreg(obj.train, family="gaussian", alpha1 = 0, cv=ncv, R=nR, penalty = "combined", split_vars = TRUE)
  fit.ridge.x <- perform_penreg(obj.train, family = "gaussian", alpha1 = 0, cv = ncv, R = nR, penalty = "combined", split_vars = FALSE)
  fit.lridge.ud <- perform_lridge(obj.train, family = "gaussian", cv = ncv, R = nR, penalty = "combined", split_vars = TRUE)
  fit.lridge.x <- perform_lridge(obj.train, family = "gaussian", cv = ncv, R = nR, penalty = "combined", split_vars = FALSE)
  fit.rlasso.ud <- perform_rlasso(obj.train, family = "gaussian", cv = ncv, R = nR, penalty = "combined", split_vars = TRUE)
  fit.rlasso.x <- perform_rlasso(obj.train, family = "gaussian", cv = ncv, R = nR, penalty = "combined", split_vars = FALSE)
  fit.rgarrote.ud <- perform_rgarrote(obj.train, family = "gaussian", cv = ncv, R = nR, penalty = "combined", split_vars = TRUE)
  fit.rgarrote.x <- perform_rgarrote(obj.train, family = "gaussian", cv = ncv, R = nR, penalty = "combined", split_vars = FALSE)
  fit.rf.notune <- ranger(y~., data = train %>% dplyr::select(-c(pred.rf.notune, pred.rf.tune)), num.trees = 1000)
  fit.rf.tune <- tuneMtryFast(y~., data = train %>% dplyr::select(-c(pred.rf.notune, pred.rf.tune)), num.trees = 1000, stepFactor = 0.5,
                              doBest = TRUE, plot = FALSE, trace = FALSE)
  
  
  # Predictions
  data.obj$pred.lasso[(1:n)[folds == i]] <- predict_penreg(obj = fit.lasso, newdata = obj.test, model = "lasso")
  data.obj$pred.ridge.ud[(1:n)[folds == i]] <- predict_penreg(obj = fit.ridge.ud,newdata = obj.test,model = "ridge")
  data.obj$pred.ridge.x[(1:n)[folds == i]] <- predict_penreg(obj = fit.ridge.x, newdata = obj.test, model = "ridge")
  data.obj$pred.lridge.ud[(1:n)[folds == i]] <- predict_lridge(obj = fit.lridge.ud, newdata = obj.test)
  data.obj$pred.lridge.x[(1:n)[folds == i]] <- predict_lridge(obj = fit.lridge.x, newdata = obj.test)
  data.obj$pred.rlasso.ud[(1:n)[folds == i]] <- predict_rlasso(obj = fit.rlasso.ud, newdata = obj.test, split_vars = TRUE)
  data.obj$pred.rlasso.x[(1:n)[folds == i]] <- predict_rlasso(obj = fit.rlasso.x, newdata = obj.test, split_vars = FALSE)
  data.obj$pred.rgarrote.ud[(1:n)[folds == i]] <- predict_rgarrote(obj = fit.rgarrote.ud, newdata = obj.test)
  data.obj$pred.rgarrote.x[(1:n)[folds == i]] <- predict_rgarrote(obj = fit.rgarrote.x, newdata = obj.test)
  dataset$pred.rf.notune[(1:n)[folds == i]] <- predict(fit.rf.notune,  data = test %>% dplyr::select(-c(y, pred.rf.notune, pred.rf.tune)))$predictions
  dataset$pred.rf.tune[(1:n)[folds == i]] <- predict(fit.rf.tune,  data = test %>% dplyr::select(-c(y, pred.rf.tune, pred.rf.notune)))$predictions
  
  # Performance
  tbl_out[1, 3:6] <- eval_performance(data.obj$pred.lasso[(1:n)[folds == i]], obj.test$y)
  tbl_out[2, 3:6] <- eval_performance(data.obj$pred.ridge.ud[(1:n)[folds == i]], obj.test$y)
  tbl_out[3, 3:6] <- eval_performance(data.obj$pred.ridge.x[(1:n)[folds == i]], obj.test$y)
  tbl_out[4, 3:6] <- eval_performance(data.obj$pred.ridge.ud[(1:n)[folds == i]], obj.test$y)
  tbl_out[5, 3:6] <- eval_performance(data.obj$pred.lridge.x[(1:n)[folds == i]], obj.test$y)
  tbl_out[6, 3:6] <- eval_performance(data.obj$pred.rlasso.ud[(1:n)[folds == i]], obj.test$y)
  tbl_out[7, 3:6] <- eval_performance(data.obj$pred.rlasso.x[(1:n)[folds == i]], obj.test$y)
  tbl_out[8, 3:6] <- eval_performance(data.obj$pred.rgarrote.ud[(1:n)[folds == i]], obj.test$y)
  tbl_out[9, 3:6] <- eval_performance(data.obj$pred.rgarrote.x[(1:n)[folds == i]], obj.test$y)
  tbl_out[10, 3:6] <- eval_performance(pred = dataset$pred.rf.notune[(1:n)[folds == i]], obs = test$y)
  tbl_out[11, 3:6] <- eval_performance(pred = dataset$pred.rf.tune[(1:n)[folds == i]], obs = test$y)
  return(tbl_out)
}

extract_varsel <- function(modelcoef){
  model.varsel <- names(modelcoef)[which(modelcoef != 0)][-1]
  model.varsel <- str_remove_all(model.varsel, "u.|d.|x.")
  return(model.varsel[!duplicated(model.varsel)])
}


# ---------- Initialization & Data --------------
data_mos <- readRDS(here::here("data", "mosaique", "Large_Mariella.rds"))
data_allpeptides <- data_mos[,7:ncol(data_mos)]
colnames(data_allpeptides) <- paste0("P",1:ncol(data_allpeptides))
rownames(data_allpeptides) <- paste0("N",1:nrow(data_allpeptides))
data_allpeptides <- apply(as.matrix(data_allpeptides), 2, as.numeric)
data_allpeptides <- scale(data_allpeptides)
clin.vars <- c("Age", "Sex", "eGFR")

# Data preparation
data_mos <- data_mos %>%
  dplyr::rename("Sex" = Gender) %>%
  mutate(Sex = to_factor(ifelse(Sex %in% "female", 1, 0)),
         Age = to_numeric(Age),
         eGFR = to_numeric(eGFR),
         log2eGFR = log(eGFR, 2))


ncv = 10
nR  = 10
vec_propz <-  rev(seq(0.1, 0.9, 0.1))


# Peptide subset selection due to computational feasibility based on max non-zero entries
list_results <- list()
for(l in 1:length(vec_propz)){
  propz <- vec_propz[l]
  ind <- which(apply(data_allpeptides, 2, max_proportionZI, prop_zero = propz))
  data_peptides <- as.matrix(data_allpeptides[,ind])
  data_mos <- data.frame(cbind(data_mos[,1:6], data_peptides)) 
  
  # Parameter
  n = nrow(data_peptides)
  p = ncol(data_peptides)
  folds <- sample(rep(1:ncv, ceiling(n/ncv)))[1:n]
  pflist <- list(c(1,2), c(2,1), c(1,3), c(3,1))
  
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
  data.obj <- list("y" = log2(data_mos$eGFR), "clinical" = c, x = ximp, u = u, d = d,
                   pred.ridge.ud = NA, pred.ridge.x = NA,
                   pred.lasso = NA,
                   pred.lridge.ud = NA, pred.lridge.x = NA,
                   pred.rgarrote.ud = NA, pred.rgarrote.x = NA,
                   pred.rf.notune = NA, pred.rf.tune = NA)
  
  
  # Data frame for random forest
  dataset <- data.frame(y = log2(data_mos$eGFR), age = data_mos$Age, sex = data_mos$Sex, x, pred.rf.notune = NA, pred.rf.tune = NA)
  
  # ----------- CV Predictions --------------------
  
  names_models <- rep(c("lasso", "ridge", "lasso-ridge", "ridge-lasso", "ridge-garrote", "random forest"), each=2)[-2]
  names_penalty <- c(rep(c("ud", "x"), 5)[-1], "without tuning", "with tuning")
  
  plan(multisession, workers=10)
  tbl_performance <- future_lapply(1:ncv, function(x) perform_CVmethods(x, data.obj = data.obj, dataset = dataset), future.seed=TRUE)
  plan(sequential)
  tbl_performance <- do.call(rbind, tbl_performance)
  
  # -------------------  Final model ----------------------
  list_model <- list()
  list_model$fit.lasso <- perform_penreg(data.obj = data.obj, family = "gaussian", cv = ncv, R = nR, alpha1 = 1, penalty = "combined", split_vars = FALSE)
  list_model$fit.ridge.ud <- perform_penreg(data.obj = data.obj, family = "gaussian", cv = ncv, R = nR, penalty = "combined", split_vars = TRUE)
  list_model$fit.ridge.x <- perform_penreg(data.obj = data.obj, family = "gaussian", cv = ncv, R = nR, penalty = "combined", split_vars = FALSE)
  list_model$fit.lridge.ud <- perform_lridge(data.obj = data.obj, family = "gaussian", cv = ncv, R = nR, penalty = "combined", split_vars = TRUE)
  list_model$fit.lridge.x <- perform_lridge(data.obj = data.obj, family = "gaussian", cv = ncv, R = nR, penalty = "combined", split_vars = FALSE)
  list_model$fit.rlasso.ud <- perform_rlasso(data.obj = data.obj, family = "gaussian", cv = ncv, R = nR, penalty = "combined", split_vars = TRUE)
  list_model$fit.rlasso.x <- perform_rlasso(data.obj = data.obj, family = "gaussian", cv = ncv, R = nR, penalty = "combined", split_vars = FALSE)
  list_model$fit.rgarrote.ud <- perform_rgarrote(data.obj = data.obj, family = "gaussian", cv = ncv, R = nR, penalty = "combined", split_vars = TRUE)
  list_model$fit.rgarrote.x <- perform_rgarrote(data.obj = data.obj, family = "gaussian", cv = ncv, R = nR, penalty = "combined", split_vars = FALSE)
  
  
  # Coefficients
  list_coef <- list()
  list_coef$lasso <- list_model$fit.lasso$fit$coefficients
  list_coef$ridge.ud <- list_model$fit.ridge.ud$fit$coefficients
  list_coef$ridge.x <- list_model$fit.ridge.x$fit$coefficients
  list_coef$lridge.ud <- list_model$fit.lridge.ud$fit$coefficients
  list_coef$lridge.x <- list_model$fit.lridge.x$fit$coefficients
  list_coef$rlasso.ud <- list_model$fit.rlasso.ud$fit$coefficients
  list_coef$rlasso.x <- list_model$fit.rlasso.x$fit$coefficients
  list_coef$rgarrote.ud <- list_model$fit.rgarrote.ud$fit$coefficients
  list_coef$rgarrote.x <- list_model$fit.rgarrote.x$fit$coefficients
  
  # Variable selection
  # List of variable selection for models 
  tbl_varsel <- data.frame(matrix(NA, ncol=14, nrow=ncol(data_allpeptides)))
  colnames(tbl_varsel) <- c("propz", "all", "included", "lasso", "ridge_ud", "ridge_x", "lridge_ud", "lridge_x",
                            "rlasso_ud", "rlasso_x", "rgarrote_ud", "rgarrote_x", "rf_notune", "rf_tune")
  tbl_varsel$propz <- propz
  tbl_varsel$all <- colnames(data_allpeptides)
  tbl_varsel$included <- ifelse(tbl_varsel$all %in% colnames(data_peptides), 1, 0)
  tbl_varsel$lasso <- ifelse(tbl_varsel$all %in% extract_varsel(list_coef$lasso), 1, 0)
  tbl_varsel$ridge_ud <- ifelse(tbl_varsel$all %in% extract_varsel(list_coef$ridge.ud), 1, 0)
  tbl_varsel$ridge_x <- ifelse(tbl_varsel$all %in% extract_varsel(list_coef$ridge.x), 1, 0)
  tbl_varsel$lridge_ud <- ifelse(tbl_varsel$all %in% extract_varsel(list_coef$lridge.ud), 1, 0)
  tbl_varsel$lridge_x <- ifelse(tbl_varsel$all %in% extract_varsel(list_coef$lridge.x), 1, 0) 
  tbl_varsel$rlasso_ud <- ifelse(tbl_varsel$all %in% extract_varsel(list_coef$rlasso.ud), 1, 0)
  tbl_varsel$rlasso_x <- ifelse(tbl_varsel$all %in% extract_varsel(list_coef$rlasso.x), 1, 0) 
  tbl_varsel$rgarrote_ud <- ifelse(tbl_varsel$all %in% extract_varsel(list_coef$rgarrote.ud), 1, 0)
  tbl_varsel$rgarrote_x <- ifelse(tbl_varsel$all %in% extract_varsel(list_coef$rgarrote.x), 1, 0) 
  tbl_varsel$rf_notune <- ifelse(tbl_varsel$all %in% colnames(data_peptides), 1, 0)
  tbl_varsel$rf_tune <- ifelse(tbl_varsel$all %in% colnames(data_peptides), 1, 0)
  
  
  tbl_perf <- tbl_performance %>%
    data.frame() %>%
    mutate(model=fct_relevel(model, names_models[!duplicated(names_models)])) %>%
    group_by(model, penalty) %>%
    summarise("R2"=mean(R2, na.rm=TRUE), "RMSPE"=mean(RMSPE, na.rm=TRUE), 
              "MAE"=mean(MAE, na.rm=TRUE), "CS"=mean(CS, na.rm=TRUE)) %>%
    ungroup() %>%
    mutate("df"=c(unlist(lapply(list_model, function(x) x$dfvars)), rep(ncol(data_peptides),2)),
           "propz" = propz) %>%
    relocate(propz, .after=penalty)
  
  tbl_performance <- tbl_performance %>%
    mutate("propz"=propz) %>%
    relocate(propz, .after=penalty)
  
  # Runtime
  mbm <- microbenchmark("lasso" = { b <- perform_penreg(data.obj = data.obj, family = "gaussian", cv = ncv, R = nR, alpha1 = 1, penalty = "combined", split_vars = FALSE)},
                        "ridge_ud" = { b <- perform_penreg(data.obj = data.obj, family = "gaussian", cv = ncv, R = nR, penalty = "combined", split_vars = TRUE)},
                        "ridge_x" = { b <- perform_penreg(data.obj = data.obj, family = "gaussian", cv = ncv, R = nR, penalty = "combined", split_vars = FALSE)},
                        "lasso-ridge_ud" = { b <- perform_lridge(data.obj = data.obj, family = "gaussian", cv = ncv, R = nR, penalty = "combined", split_vars = TRUE)},
                        "lasso-ridge_x" = { b <- perform_lridge(data.obj = data.obj, family = "gaussian", cv = ncv, R = nR, penalty = "combined", split_vars = FALSE)},
                        "ridge-lasso_ud" = { b <- perform_rlasso(data.obj = data.obj, family = "gaussian", cv = ncv, R = nR, penalty = "combined", split_vars = TRUE)},
                        "ridge-lasso_x" = { b <- perform_rlasso(data.obj = data.obj, family = "gaussian", cv = ncv, R = nR, penalty = "combined", split_vars = FALSE)},
                        "ridge-garrote_ud" = { b <- perform_rgarrote(data.obj = data.obj, family = "gaussian", cv = ncv, R = nR, penalty = "combined", split_vars = TRUE)},
                        "ridge-garrote_x" = { b <- perform_rgarrote(data.obj = data.obj, family = "gaussian", cv = ncv, R = nR, penalty = "combined", split_vars = FALSE)},
                        "rf_notune" = { b <- ranger(y~., data = dataset %>% dplyr::select(-c(pred.rf.notune, pred.rf.tune)), num.trees = 1000)},
                        "rf_tune" =  { b <- tuneMtryFast(y~., data = dataset %>% dplyr::select(-c(pred.rf.notune, pred.rf.tune)), num.trees = 1000, 
                                                         stepFactor = 0.5, doBest = TRUE, plot = FALSE, trace = FALSE)}, times = 1
                        )
  mbm <- mbm %>%
    data.frame() %>%
    mutate(propz = propz,
           model = str_split_fixed(expr, "_", 2)[,1],
           penalty = str_split_fixed(expr, "_", 2)[,2],
           model = ifelse(model %in% "rf", "random forest", model),
           penalty = ifelse(penalty %in% "notune", "without tuning", penalty),
           penalty = ifelse(penalty %in% "tune", "with tuning", penalty),
           penalty = ifelse(penalty %in% "", "x", penalty))
  
  
# ------------ Save results ---------------
  list_results[[l]] <- list("performance" = tbl_perf, 
                       "performance_ncv" = tbl_performance, 
                       "varsel" = tbl_varsel,
                       "coef" = list_coef,
                       "models" = list_model,
                       "mbm" = mbm)
  filename <- paste0("results_mosaique_pnz", propz, "_nR", nR,"_ncv", ncv,".rds")
  saveRDS(list_results[[l]], here::here("output", "example", filename))
}

# --- Concatenate results across proportions of non-zero values ----
res.files <- list.files(here::here("output", "example"), pattern = "results_")
list_results <- lapply(res.files, function(x) readRDS(here::here("output", "example", x)))
tbl_performance <- do.call(rbind, lapply(list_results, function(x) x$performance))
tbl_performance_ncv <- do.call(rbind, lapply(list_results, function(x) x$performance_ncv))
tbl_mbm <- do.call(rbind, lapply(list_results, function(x) x$mbm))
tbl_varsel <- do.call(rbind, lapply(list_results, function(x) x$varsel))
list_coef <- lapply(list_results, function(x) x$coef)

results <- list("performance"=tbl_performance, "performance_ncv"=tbl_performance_ncv, "mbm"=tbl_mbm, "varsel"=tbl_varsel, "coef"=list_coef)
filename <-  paste0("results_mosaique_all_nR", nR,"_ncv", ncv,".rds")
saveRDS(results, here::here("output", "example", filename))


