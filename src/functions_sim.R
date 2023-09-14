#' ==============================================================================
#' Author: MG
#' Date: 09/05/2023
#' Info: Functions only relevant for the Monte Carlo simulation
#' ==============================================================================


# =========================== Run scenario =====================================

simulate_scenario <- function(scn, dsgn){
  
  ## Remove later
  # nr=1; scn=scenarios[nr,]; dsgn=sim_design[[scenarios[nr,]$dsgn]]
  
  filename <- paste0("sim_i", scn$iter, "_scn", scn$scenario, "_n", scn$n, "_p", scn$p, "_beta", scn$beta_max, "_a", scn$a, 
                     "_eps", scn$epslvl, "_propzi", scn$propzi, "_revzi", tolower(scn$revzi), "_struczero", scn$struczero, ".rds")

  # Generate large validation dataset
  data.val <- data_generation(dsgn = dsgn, scenario = scn$scenario, n = 100000, p=scn$p, beta_max = scn$beta_max, a = scn$a, epsstd = scn$epsstd, 
                              propzi = scn$propzi, revzi = scn$revzi, struczero = scn$struczero)
  data.val <- generate_dataobj(y = data.val$data_ana$y, x = data.val$data_ana$x, clinical = NULL)
  
  # Run replications
  scn_res <-  lapply(1:scn$iter, function(x){
    data_iter <- data_generation(dsgn = dsgn, scenario = scn$scenario, n = scn$n, p = scn$p, beta_max = scn$beta_max, a = scn$a, epsstd = scn$epsstd, 
                                 propzi=scn$propzi, revzi = scn$revzi, struczero=scn$struczero)
    res_iter <-  data_analysis(df = data_iter, data.val = data.val, n = scn$n, p = scn$p, 
                               ncv = 10, nR = 2, nlams = 10)
    data_iter <- mapply(cbind, data_iter, "i" = x, SIMPLIFY = F)
    res_iter <- mapply(cbind, res_iter, "i" = x, SIMPLIFY = F)
    
    return(c(data_iter, res_iter)) 
  })
  
  # Summarize
  summarize_scenario(filename = filename, scn = scn, scn_res = scn_res)
  
  return(NULL)
  
}


# ============================ Data generation =================================
data_generation <- function(dsgn, scenario, n, p, beta_max, a, epsstd, propzi, revzi, struczero){
  
  # # Remove later (parameter input for function)
  # dsgn=dsgn; scenario = scn$scenario; n=100000; p=scn$p; ngroups=scn$ngroups; a=scn$a; epsstd=scn$epsstd;
  # propzi=scn$propzi; struczero=scn$struczero; beta_max=scn$beta_max
  
  # Parameter
  ngroups = 4
  
  # Groupwise Hub correlation design filled with toeplitz
  data_logvars <- simulate_data(dsgn, n)

  # Structural zeros for outcome generation
  groupzi <- seq(0,propzi, length.out = p/ngroups)
  if(revzi){groupzi <- rev(groupzi)}
  seq_propzi <- rep(groupzi, ngroups)
  seq_strucz <- seq_propzi * struczero

  D_struc <- as.matrix(sapply(1-seq_strucz, function(x) rbinom(n, size = 1, prob = x))) # per group linear increase in zero-inflation 0-80%
  X <- as.matrix(data_logvars * D_struc)
  
  # Extract hubs and indices of true predictors (scenario D)
  groupsize <- rep(p/ngroups, ngroups)
  hubindex <- cumsum(groupsize) - groupsize + 1 # identify index of hub
  groupindex <- cbind(hubindex, hubindex + groupsize - 1)
  
  
  if(scenario %in% "A"){
    # Scenario A
    nelem <- c(groupsize[1], 0, 0, 0)
    truepredindex <- 1:groupsize[1]
    ptrue <- length(truepredindex)    
    
    # Generate coeffs and error 
    beta_X <- beta_D <- rep(0, p)
    beta_X[truepredindex] <- rev(seq(0, beta_max, length.out = ptrue+1)[-1])
    beta_D[truepredindex] <-  sample(seq(0, beta_max, length.out = ptrue+1)[-1])  
  }else if(scenario %in% "B"){
    # Scenario B
    nelem <- c(5, 5, 5, 5)  # number of true predictors in each group
    truepredindex <- c(sapply(1:ngroups, function(x) seq(groupindex[x,1], groupindex[x,2], 5)[1:nelem[x]]))
    ptrue <- length(truepredindex)
    
    # Generate coeffs and error 
    beta_X <- beta_D <- rep(0, p)
    beta_X[truepredindex] <- c(sapply(nelem, function(x) rev(seq(0, beta_max, length.out = x+1)[-1])))
    beta_D[truepredindex] <-  c(sapply(nelem, function(x) sample(seq(0, beta_max, length.out = x+1)[-1], x)))    
  }else{
    stop("Scenario must be 'A' or 'B'")
  }

  # Generate outcome: a controls influence of X and D components
  eps <- rnorm(n, mean = 0, sd = epsstd)
  y <- c(a * X %*% beta_X + (1-a) * D_struc %*% beta_D  + eps)
  
  # Structural and sampling zeros for data analyst
  sampzero <- 1 - struczero
  seq_sampz <- seq_propzi * sampzero
  sampz_thresh <- sapply(1:length(seq_sampz), function(i) quantile(data_logvars[,i], seq_sampz[i]))
  D_samp <- sapply(1:ncol(data_logvars), function(j) (data_logvars[,j] > sampz_thresh[j])*1)
  Xs <- (data_logvars * D_struc) * D_samp
  # plot(x=1:ncol(Xs), sort(apply(Xs,2, function(x) sum(x==0)/length(x))))
  
  out <- list("data_ana" = list("y" = y, "x" = Xs),
              "data_gen" = data.frame("y" = y, "X_true" = X, "D_struc" = D_struc, "D_samp" = D_samp,"eps" = eps),
              "true_coef" = data.frame("beta_X" = beta_X, "beta_D" = beta_D),
              "info_zi" = data.frame("vind" = 1:length(seq_propzi), "propzi" = seq_propzi, "struczi" = seq_strucz, "samplzi" = seq_sampz))
  return(out)
}




# ============================ Data analysis =================================
data_analysis <- function(df, data.val, n, p, ncv=10, nR=10, nlams=100, pflist=list(c(1,2), c(2,1), c(1,3))){
  
  # Remove later
  # df = data_iter
  # data.val = data.val
  # n=scn$n
  # p=scn$p
  # ngroups=scn$ngroups
  # ncv=10; nR=2; nlams=10; pflist=list(c(1,2), c(2,1), c(1,3))

  # Parameter
  ngroups = 4
  
  # Prepare data
  data.obj <- generate_dataobj(y = df$data_ana$y, x = df$data_ana$x, clinical = NULL)
  
  # CV 
  methnames <- c("oracle", "lasso", "ridge", "lasso-ridge", "ridge-lasso", "ridge-garrote", "random forest")
  tbl_perf <- data.frame("model" = methnames,
                         "penalty" = c("-",rep(c("combined"), times = 5), "-"),
                         R2 = NA, RMSPE = NA, MAE = NA, C = NA, CS = NA)
  tbl_coef <- cbind.data.frame("var"=paste0("V", 1:nrow(df$true_coef)), df$true_coef)
  
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
  
  # --- Lasso
  fit.lasso.cv <- perform_penreg(data.obj, family = "gaussian", alpha1 = 1, nl1 = nlams, cv = ncv, R = nR, penalty = "combined")
  data.val$pred.lasso <- predict_penreg(obj = fit.lasso.cv, newdata = data.val, model = "lasso")
  tbl_perf[2, 3:7] <- eval_performance(pred = data.val$pred.lasso, obs = data.val$y)
  tbl_coef$beta_lasso_cb <- fit.lasso.cv$coefficients[-1]
  
  # --- Ridge 
  fit.ridge.cb <- perform_penreg(data.obj, family = "gaussian",  alpha = 0, nl1 = nlams, cv = ncv, R = nR, 
                                 penalty = "combined", split_vars = TRUE)
  data.val$pred.ridge.cb <- predict_penreg(obj = fit.ridge.cb, newdata = data.val, model = "ridge")
  tbl_perf[3, 3:7] <- eval_performance(pred = data.val$pred.ridge.cb, obs = data.val$y)
  tbl_coef$beta_ridge_cb_u <- fit.ridge.cb$coefficients[2:(p+1)]
  tbl_coef$beta_ridge_cb_d <- fit.ridge.cb$coefficients[(p+2):(2*p+1)]
  
  
  # --- Lasso-ridge 
  fit.lridge.cb <- perform_lridge(data.obj, family = "gaussian", cv = ncv, R = nR, nlambda = rep(nlams, 2), 
                                  penalty = "combined", split_vars = TRUE)
  data.val$pred.lridge.cb <- predict_lridge(obj = fit.lridge.cb, newdata = data.val)
  tbl_perf[4, 3:7] <- eval_performance(pred = data.val$pred.lridge.cb, obs = data.val$y)
  tbl_coef$beta_lridge_cb_u <- fit.lridge.cb$coefficients[2:(p+1)]
  tbl_coef$beta_lridge_cb_d <- fit.lridge.cb$coefficients[(p+2):(2*p+1)]

  
  # --- Ridge-lasso 
  fit.rlasso.cb <- perform_rlasso(data.obj, family = "gaussian", cv = ncv, R = nR, nlambda = rep(nlams, 2), 
                                  penalty = "combined", split_vars = TRUE)
  data.val$pred.rlasso.cb <- predict_rlasso(obj = fit.rlasso.cb, newdata = data.val)
  tbl_perf[5, 3:7] <- eval_performance(pred = data.val$pred.rlasso.cb, obs = data.val$y)
  tbl_coef$beta_rlasso_cb_u <- fit.rlasso.cb$coefficients[2:(p+1)]
  tbl_coef$beta_rlasso_cb_d <- fit.rlasso.cb$coefficients[(p+2):(2*p+1)]
  
  # --- Ridge-garrote
  fit.rgarrote.cb <- perform_rgarrote(data.obj, family = "gaussian", cv = ncv, R = nR, nlambda = rep(nlams, 2),
                                      penalty = "combined", alpha1 = 0, split_vars = TRUE)
  data.val$pred.rgarrote.cb <- predict_rgarrote(obj = fit.rgarrote.cb, newdata = data.val)
  tbl_perf[6, 3:7] <- eval_performance(pred = data.val$pred.rgarrote.cb, obs = data.val$y)
  tbl_coef$beta_rgarrote_cb_u <- fit.rgarrote.cb$coefficients[2:(p+1)]
  tbl_coef$beta_rgarrote_cb_d <- fit.rgarrote.cb$coefficients[(p+2):(2*p+1)]
  
  # --- Random forest
  train <- data.frame(y = data.obj$y, data.obj$x)
  fit.rf <- ranger(y~., data = train, num.trees = 1000)
  data.val$pred.rf <- predict(fit.rf, data = data.val$x)$predictions
  tbl_perf[7, 3:7] <- eval_performance(pred = data.val$pred.rf , obs = data.val$y)
  tbl_perf$relRMSPE <- tbl_perf$RMSPE/min(tbl_perf$RMSPE)
  
  # Merge results
  # Variable selection
  groupsize <- p/ngroups
  list_models <- list(fit.lasso.cv, 
                      fit.ridge.cb, 
                      fit.rlasso.cb,
                      fit.lridge.cb,
                      fit.rgarrote.cb)
  list_sel <- list()
  for(l in 1:length(list_models)){
    fit.model <-  list_models[[l]]
    list_sel[[l]] <- eval_selection(model = attr(fit.model, "class"), penalty = attr(fit.model, "penalty"), varnames = tbl_coef$var,
                                    true_coef = df$true_coef, pred_coef = fit.model$coefficients, ngroups = ngroups, p = p)
  }
  tbl_varsel <- do.call(rbind, lapply(list_sel, function(x) x[[1]]))
  tbl_groupsel <- do.call(rbind, lapply(list_sel, function(x) x[[2]])) 
  tbl_allsel <- do.call(rbind, lapply(list_sel, function(x) x[[3]])) 
  
  out <- list("est_perf" = tbl_perf, 
              "est_coef" = tbl_coef,
              "est_varsel" = tbl_varsel, 
              "est_groupsel" = tbl_groupsel,
              "est_allsel" = tbl_allsel)
  return(out)
}

# ============================ Summarize scenario =================================

summarize_scenario <- function(filename, scn, scn_res){
 
  options(dplyr.summarise.inform = FALSE)
  # Concatenate results in corresponding oject
  tbl_iters_dataana <- do.call(rbind, lapply(scn_res, function(x) data.frame(cbind("y" = x$data_ana[1,1][[1]], x$data_ana[2,1][[1]], "i" = x$data_ana[1,2][[1]]))))
  tbl_iters_datagen <- do.call(rbind, lapply(scn_res, function(x) x$data_gen))
  tbl_iters_performance <- do.call(rbind, lapply(scn_res, function(x) x$est_perf))
  tbl_iters_varsel <- do.call(rbind, lapply(scn_res, function(x) x$est_varsel))
  tbl_iters_groupsel <- do.call(rbind, lapply(scn_res, function(x) x$est_groupsel))
  tbl_iters_allsel <- do.call(rbind, lapply(scn_res, function(x) x$est_allsel))
  tbl_iters_truecoef <- do.call(rbind, lapply(scn_res, function(x) x$true_coef)) %>% merge(scn, . )
  tbl_iters_estcoef <- do.call(rbind, lapply(scn_res, function(x) x$est_coef)) %>% merge(scn, . )
  niter <- max(tbl_iters_varsel$i)
  
  # -- Oracle 
  tbl_performance <- tbl_iters_performance %>%
    data.frame() %>%
    group_by(model, penalty) %>%
    summarise("RMSPE.est" = mean(RMSPE, na.rm = T), "RMSPE.med" = median(RMSPE, na.rm = T), 
              "RMSPE.lo" = quantile(RMSPE, 0.05, na.rm = T), "RMSPE.up" = quantile(RMSPE, 0.95, na.rm = T),
              "relRMSPE.est" = mean(relRMSPE, na.rm = T), "relRMSPE.med" = median(relRMSPE, na.rm = T), 
              "relRMSPE.lo" = quantile(relRMSPE, 0.05, na.rm = T), "relRMSPE.up" = quantile(relRMSPE, 0.95, na.rm = T),
              "R2.est" = mean(R2, na.rm = T), "R2.med" = median(R2, na.rm = T),
              "R2.lo" = quantile(R2, 0.05, na.rm = T), "R2.up" = quantile(R2, 0.95, na.rm = T),
              "CS.est" = mean(CS, na.rm = T), "CS.med" = median(CS, na.rm = T),
              "CS.lo" = quantile(CS, 0.05, na.rm = T), "CS.up" = quantile(CS, 0.95, na.rm = T),
              "MAE.est" = mean(MAE, na.rm = T), "MAE.med" = median(MAE, na.rm = T),
              "MAE.lo" = quantile(MAE, 0.05, na.rm = T), "MAE.up" = quantile(MAE, 0.95, na.rm = T),
              "C.est" = mean(C, na.rm = T), "C.med" = median(C, na.rm = T),
              "C.lo" = quantile(C, 0.05, na.rm = T), "C.up" = quantile(C, 0.95, na.rm = T)) %>%
    data.frame()  %>%
    merge(scn, . )
  
  
  tbl_varsel <- tbl_iters_varsel %>% 
    group_by(model, penalty, varname) %>% 
    summarise("truepos_var" = sum(truepos_var)/niter, 
              "trueneg_var" = sum(trueneg_var)/niter,
              "falsepos_var" = sum(falsepos_var)/niter, 
              "falseneg_var" = sum(falseneg_var)/niter) %>%
    data.frame() %>%
    arrange(varname) %>%
    merge(scn, . )
  
  tbl_groupsel <- tbl_iters_groupsel %>% 
    group_by(model, penalty, group) %>% 
    summarise("truepos_any" = sum(truepos_any)/niter, 
              "truepos_group" = sum(truepos_group)/niter, 
              "trueneg_group" = sum(trueneg_group)/niter,
              "falsepos_group" = sum(falsepos_group)/niter, 
              "falseneg_group" = sum(falseneg_group)/niter) %>%
    data.frame()  %>%
    merge(scn, . )
  
  tbl_allsel <- tbl_iters_allsel %>% 
    group_by(model, penalty) %>% 
    summarise("FPDR" = sum(FPDR)/niter, 
              "FNDR" = sum(FNDR)/niter,
              "TPDR" = sum(TPDR)/niter, 
              "TNDR" = sum(TNDR)/niter) %>%
    data.frame()  %>%
    merge(scn, . )
  
  tbl_iters_performance <- tbl_iters_performance %>% merge(scn, . )
  tbl_iters_varsel <- tbl_iters_varsel %>% merge(scn, . )
  tbl_iters_groupsel <- tbl_iters_groupsel %>% merge(scn, . )
  tbl_iters_allsel <- tbl_iters_allsel %>% merge(scn, . )
  
  # Save results
  list_results <- list("results" = list("performance" = tbl_performance, "varsel" = tbl_varsel, "groupsel" = tbl_groupsel, "allsel" = tbl_allsel),
                       "iters" = list("performance" = tbl_iters_performance, "varsel" = tbl_iters_varsel, "groupsel" = tbl_iters_groupsel, "allsel" = tbl_iters_allsel),
                       "data" = list("gen" = tbl_iters_datagen, "ana" = tbl_iters_dataana),
                       "coef" = list("truecoef" = tbl_iters_truecoef, "estcoef" = tbl_iters_estcoef))
  saveRDS(list_results, here::here(sim.path, paste0(filename , ".rds")))  
}

# ============================ Save all scenarios =================================

evaluate_scenarios <- function(sim.files, sim.path){
  #' Preparation of simulation results for analysis
  
  res <- lapply(sim.files, function(x) readRDS(here::here(sim.path, x)))
  tbl_perf <- do.call(rbind, lapply(res, function(x) x$results$performance))
  tbl_varsel <- do.call(rbind, lapply(res, function(x) x$results$varsel))
  tbl_groupsel <- do.call(rbind, lapply(res, function(x) x$results$groupsel))
  tbl_allsel <- do.call(rbind, lapply(res, function(x) x$results$allsel))
  
  tbl_iters_perf <- do.call(rbind, lapply(res, function(x) x$iters$performance))
  tbl_iters_varsel <- do.call(rbind, lapply(res, function(x) x$iters$varsel))
  tbl_iters_groupsel <- do.call(rbind, lapply(res, function(x) x$iters$groupsel))
  tbl_iters_allsel <- do.call(rbind, lapply(res, function(x) x$iters$allsel))
  
  tbl_estcoef <- do.call(rbind, lapply(res, function(x) x$coef$estcoef))

  # Results per iteration
  tbl_res <- list("performance"=list("tbl_perf"=tbl_perf, "tbl_varsel"=tbl_varsel, "tbl_groupsel"=tbl_groupsel, "tbl_allsel"=tbl_allsel), 
                  "iters" = list("tbl_iters_perf"=tbl_iters_perf, "tbl_iters_varsel"=tbl_iters_varsel, 
                                 "tbl_iters_groupsel"=tbl_iters_groupsel, "tbl_iters_allsel"=tbl_iters_allsel),
                  "coef" = tbl_estcoef)
  saveRDS(tbl_res, here::here(sim.path, "tbl_scenario_results.rds"))

  return(NULL)
}
