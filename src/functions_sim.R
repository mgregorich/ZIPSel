#' ==============================================================================
#' Author: MG
#' Date: 09/05/2023
#' Info: Functions only relevant for the Monte Carlo simulation
#' ==============================================================================


# =========================== Run scenario =====================================

simulate_scenario <- function(scn, dsgn){
  
  ## Remove later
  # nr=100; scn=scenarios[nr,]; dsgn=sim_design[[scenarios[nr,]$dsgn]]

  filename <- paste0("sim_i", scn$iter, "_scn", scn$scenario, "_n", scn$n, "_p", scn$p,"_beta", scn$beta_max, "_UDdepen", scn$UDdepen, 
                     "_eps", scn$epslvl, "_propzi", scn$propzi, "_revzi", tolower(scn$revzi), "_struczero", scn$struczero, ".rds")

  # Generate large validation dataset
  data.val <- data_generation(dsgn = dsgn, scenario = scn$scenario, n = 100000, p=scn$p, beta_max = scn$beta_max, 
                              a = scn$a, epsstd = scn$epsstd, xmean = scn$xmean, xstd = scn$xstd,
                              propzi = scn$propzi, revzi = scn$revzi, struczero = scn$struczero)
  data.val <- generate_dataobj(y = data.val$data_ana$y, x = data.val$data_ana$X, 
                               clinical = NULL, logtransform = FALSE)
  
  # Run replications in parallel
  scn_res <- future_lapply(1:scn$iter, function(x){
    data_iter <- data_generation(dsgn = dsgn, scenario = scn$scenario, n = scn$n, p = scn$p, beta_max = scn$beta_max, 
                                 a = scn$a, epsstd = scn$epsstd, xmean = scn$xmean, xstd = scn$xstd,
                                 propzi=scn$propzi, revzi = scn$revzi, struczero=scn$struczero)
    res_iter <-  data_analysis(df = data_iter, data.val = data.val, n = scn$n, p = scn$p, 
                               ncv = 10, nR = 1, nlams = 50)
    data_iter <- mapply(cbind, data_iter, "i" = x, SIMPLIFY = F)
    res_iter <- mapply(cbind, res_iter, "i" = x, SIMPLIFY = F)
    
    return(c(data_iter, res_iter)) 
  }, future.seed = TRUE)
  
  # Summarize
  summarize_scenario(filename = filename, scn = scn, scn_res = scn_res)
}


# ============================ Data generation =================================
data_generation <- function(dsgn, scenario, n, p, beta_max, a, epsstd, xmean, xstd, propzi, revzi, struczero){
  
  # # Remove later (parameter input for function)
  # dsgn=dsgn; scenario = scn$scenario; n=scn$n; p=scn$p; ngroups=scn$ngroups; a=scn$a; epsstd=scn$epsstd;
  # xmean=scn$xmean; xstd=scn$xstd;
  # propzi=scn$propzi; struczero=scn$struczero; beta_max=scn$beta_max; revzi=scn$revzi

  # Parameter
  ngroups = 4
  
  # Draw data from simdesign
  data_logvars <- simulate_data(dsgn, n)

  # Structural zeros for outcome generation
  groupzi <- seq(0, propzi, length.out = p/ngroups)
  if(revzi){groupzi <- rev(groupzi)}
  seq_propzi <- rep(groupzi, ngroups)
  seq_strucz <- seq_propzi * struczero

  D_struc <- sapply(1-seq_strucz, function(x) rbinom(n, size = 1, prob = x), simplify = "matrix") # per group linear increase in zero-inflation 0-80%
  X_org <- as.matrix(data_logvars * D_struc)
  struczero_X <- apply(X_org, 2, function(x) sum(x==0)/nrow(X_org))
  
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
    beta_X[truepredindex] <- rev(seq(0.1, beta_max, length.out = ptrue+1)[-1])
    beta_D[truepredindex] <-  sample(seq(0.1, beta_max, length.out = ptrue+1)[-1])  
  }else if(scenario %in% "B"){
    # Scenario B
    nelem <- c(5, 5, 5, 5)  # number of true predictors in each group
    truepredindex <- floor(c(sapply(1:ngroups, function(x) seq(groupindex[x,1], groupindex[x,2], length.out=5)[1:nelem[x]])))
    ptrue <- length(truepredindex)
    
    # Generate coeffs and error 
    beta_X <- beta_D <- rep(0, p)
    beta_X[truepredindex] <- c(sapply(nelem, function(x) rev(seq(1, beta_max, length.out = x+1)[-1])))
    beta_D[truepredindex] <-  c(sapply(nelem, function(x) sample(seq(1, beta_max, length.out = x+1)[-1], x)))    
  }else{
    stop("Scenario must be 'A' or 'B'")
  }

  # Generate outcome: a controls influence of X and D components
  X_struc <- X_org
  X_struc[X_struc>0] <- log(X_struc[X_struc>0])
  X_struc[X_struc==0] <- xmean
  eps <- rnorm(n, mean = 0, sd = epsstd)
  y <- a * c(X_struc %*% beta_X) + (1-a) * c(D_struc %*% beta_D) + eps
  
  
  # Structural and sampling zeros for data analyst
  sampzero <- 1 - struczero
  seq_sampz <- seq_propzi * sampzero
  sampz_thresh <- sapply(1:length(seq_sampz), function(i) qlnorm(seq_sampz[i], meanlog = xmean, sdlog = xstd))
  # sampz_thresh <- sapply(1:length(seq_sampz), function(i) quantile(data_logvars[,i], seq_sampz[i]))
  D_samp <- sapply(1:ncol(data_logvars), function(j) (data_logvars[,j] > sampz_thresh[j])*1)
  X_samp <- as.matrix((X_struc * D_struc) *D_samp)
  samzero_X <- apply(D_samp, 2, function(x) sum(x==0)/nrow(D_samp))
  zero_X <- apply(X_samp, 2, function(x) sum(x==0)/nrow(X_samp))
  
  out <- list("data_ana" = list("y" = y, "X" = X_samp),
              "data_gen" = list("y" = y, "X_org" = X_org, "X_struc" = X_struc,
                                "D_struc" = D_struc, "D_samp" = D_samp, "eps"=eps),
              "data_coef" = data.frame("var" = paste0("V",1:length(seq_propzi)), "beta_X" = beta_X, "beta_D" = beta_D, 
                                       "propzi" = zero_X, "struczi" = struczero_X, "samplzi" = samzero_X))
  return(out)
}




# ============================ Data analysis =================================
data_analysis <- function(df, data.val, n, p, ncv, nR, nlams, pflist=list(c(1,2), c(2,1), c(1,3))){
  
  # Remove later
  # df = data_iter
  # data.val = data.val
  # n=scn$n
  # p=scn$p
  # ngroups=scn$ngroups
  # ncv=10; nR=1; nlams=10; pflist=NULL

  # Parameter
  ngroups = 4
  
  # Prepare data
  data.obj <- generate_dataobj(y = df$data_ana$y, x = df$data_ana$X, 
                               clinical = NULL, logtransform = FALSE)
  true_beta <- df$data_coef
  
  # CV 
  methnames <- c("oracle", "lasso", "ridge", "lasso-ridge", "ridge-lasso", "ridge-garrote")
  tbl_perf <- data.frame("model" = methnames,
                         "penalty" = c("-",rep(c("combined"), times = 5)),
                         R2 = NA, RMSPE = NA, MAE = NA, CS = NA, df = NA, extime=NA)
  tbl_coef <- cbind.data.frame("var"=paste0("V", 1:nrow(true_beta)), true_beta)
  tbl_pred <- data.frame("y" = data.val$y)
  
  # --- Oracle
  start <- Sys.time()
  X <- cbind(data.obj$u, data.obj$d)
  beta_true <- c(true_beta$beta_X, true_beta$beta_D)
  true_pred <- X[,which(beta_true!=0)]
  fit.oracle <- lm(data.obj$y~true_pred)
  beta_oracle <- rep(0, ncol(X)+1)
  beta_oracle[c(1,which(beta_true!=0)+1)] <- coef(fit.oracle)
  beta_oracle[is.na(beta_oracle)] <- 0 
  tbl_pred$pred.oracle <-  beta_oracle[1] + c(cbind(data.val$u, data.val$d) %*% beta_oracle[-1])
  tbl_perf[1, 3:6] <- eval_performance(pred = tbl_pred$pred.oracle, obs = tbl_pred$y)
  tbl_perf[1, "df"] <- floor(sum(beta_oracle[2:(2*(p+1))]!=0, na.rm=TRUE)/2)
  tbl_coef$beta_oracle_u <- beta_oracle[2:(p+1)]
  tbl_coef$beta_oracle_d <- beta_oracle[(p+2):(2*p+1)]
  fit.oracle$coefficients <- beta_oracle
  attr(fit.oracle, "class") <- "oracle"
  attr(fit.oracle, "penalty") <- "-"
  end <- Sys.time()
  tbl_perf[1, "extime"] <- as.numeric(end-start, units="secs")
  
  
  # --- Lasso
  fit.lasso <- perform_penreg(data.obj, family = "gaussian", alpha1 = 1, nl1 = nlams, cv = ncv, R = nR, penalty = "combined")
  tbl_pred$pred.lasso <- predict_penreg(obj = fit.lasso, newdata = data.val, model = "lasso")
  tbl_perf[2, 3:6] <- eval_performance(pred = tbl_pred$pred.lasso, obs = tbl_pred$y)
  tbl_perf[2, "df"] <- fit.lasso$dfvars
  tbl_perf[2, "extime"] <- fit.lasso$extime
  tbl_coef$beta_lasso <- fit.lasso$coefficients[-1]
  
  # --- Ridge 
  fit.ridge.cb <- perform_penreg(data.obj, family = "gaussian",  alpha = 0, nl1 = nlams, cv = ncv, R = nR, 
                                 penalty = "combined", split_vars = TRUE)
  tbl_pred$pred.ridge.cb <- predict_penreg(obj = fit.ridge.cb, newdata = data.val, model = "ridge")
  tbl_perf[3, 3:6] <- eval_performance(pred = tbl_pred$pred.ridge.cb, obs = tbl_pred$y)
  tbl_perf[3, "df"] <- fit.ridge.cb$dfvars
  tbl_perf[3, "extime"] <- fit.ridge.cb$extime
  tbl_coef$beta_ridge_cb_u <- fit.ridge.cb$coefficients[2:(p+1)]
  tbl_coef$beta_ridge_cb_d <- fit.ridge.cb$coefficients[(p+2):(2*p+1)]
  
  
  # --- Lasso-ridge 
  fit.lridge.cb <- perform_lridge(data.obj, family = "gaussian", cv = ncv, R = nR, nlambda = rep(nlams, 2), 
                                  penalty = "combined", split_vars = TRUE)
  tbl_pred$pred.lridge.cb <- predict_lridge(obj = fit.lridge.cb, newdata = data.val)
  tbl_perf[4, 3:6] <- eval_performance(pred = tbl_pred$pred.lridge.cb, obs = tbl_pred$y)
  tbl_perf[4, "df"] <- fit.lridge.cb$dfvars
  tbl_perf[4, "extime"] <- fit.lridge.cb$extime
  tbl_coef$beta_lridge_cb_u <- fit.lridge.cb$coefficients[2:(p+1)]
  tbl_coef$beta_lridge_cb_d <- fit.lridge.cb$coefficients[(p+2):(2*p+1)]

  
  # --- Ridge-lasso 
  fit.rlasso.cb <- perform_rlasso(data.obj, family = "gaussian", cv = ncv, R = nR, nlambda = rep(nlams, 2), 
                                  penalty = "combined", split_vars = TRUE)
  tbl_pred$pred.rlasso.cb <- predict_rlasso(obj = fit.rlasso.cb, newdata = data.val)
  tbl_perf[5, 3:6] <- eval_performance(pred = tbl_pred$pred.rlasso.cb, obs = tbl_pred$y)
  tbl_perf[5, "df"] <- fit.rlasso.cb$dfvars
  tbl_perf[5, "extime"] <- fit.rlasso.cb$extime
  tbl_coef$beta_rlasso_cb_u <- fit.rlasso.cb$coefficients[2:(p+1)]
  tbl_coef$beta_rlasso_cb_d <- fit.rlasso.cb$coefficients[(p+2):(2*p+1)]
  
  # --- Ridge-garrote
  fit.rgarrote.cb <- perform_rgarrote(data.obj, family = "gaussian", cv = ncv, R = nR, nlambda = rep(nlams, 2),
                                      penalty = "combined", alpha1 = 0, split_vars = TRUE)
  tbl_pred$pred.rgarrote.cb <- predict_rgarrote(obj = fit.rgarrote.cb, newdata = data.val)
  tbl_perf[6, 3:6] <- eval_performance(pred = tbl_pred$pred.rgarrote.cb, obs = tbl_pred$y)
  tbl_perf[6, "df"] <- fit.rgarrote.cb$dfvars
  tbl_perf[6, "extime"] <- fit.rgarrote.cb$extime
  tbl_coef$beta_rgarrote_cb_u <- fit.rgarrote.cb$coefficients[2:(p+1)]
  tbl_coef$beta_rgarrote_cb_d <- fit.rgarrote.cb$coefficients[(p+2):(2*p+1)]
  
  # --- Random forest
  # train <- data.frame(y = data.obj$y, data.obj$x)
  # fit.rf <- ranger(y~., data = train, num.trees = 500, num.threads = 1)
  # data.val$pred.rf <- predict(fit.rf, data = data.val$x, num.threads = 1)$predictions
  # tbl_perf[7, 3:6] <- eval_performance(pred = data.val$pred.rf , obs = data.val$y)
  # fit.rf$coefficients <- rep(1, p+1)
  # attr(fit.rf, "class") <- "random forest"
  # attr(fit.rf, "penalty") <- "-"
  
  # Merge results
  # Variable selection
  groupsize <- p/ngroups
  list_models <- list(fit.oracle,
                      fit.lasso, 
                      fit.ridge.cb, 
                      fit.rlasso.cb,
                      fit.lridge.cb,
                      fit.rgarrote.cb)
  list_sel <- list()
  for(l in 1:length(list_models)){
    fit.model <-  list_models[[l]]
    list_sel[[l]] <- eval_selection(model = attr(fit.model, "class"), penalty = attr(fit.model, "penalty"), varnames = tbl_coef$var,
                                    true_coef = true_beta, pred_coef = fit.model$coefficients, ngroups = ngroups, p = p)
  }
  tbl_varsel <- do.call(rbind, lapply(list_sel, function(x) x[[1]]))
  tbl_groupsel <- do.call(rbind, lapply(list_sel, function(x) x[[2]])) 
  tbl_allsel <- do.call(rbind, lapply(list_sel, function(x) x[[3]])) 
  tbl_perf[,"relRMSPE"] <- tbl_perf$RMSPE/min(tbl_perf$RMSPE)
  
  out <- list("est_perf" = tbl_perf, 
              "est_coef" = tbl_coef,
              "est_pred" = tbl_pred,
              "est_varsel" = tbl_varsel, 
              "est_groupsel" = tbl_groupsel,
              "est_allsel" = tbl_allsel)
  return(out)
}

# ============================ Summarize scenario =================================

summarize_scenario <- function(filename, scn, scn_res){
 
  options(dplyr.summarise.inform = FALSE)
  # Concatenate results in corresponding object
  tbl_iters_performance <- do.call(rbind, lapply(scn_res, function(x) x$est_perf)) 
  tbl_iters_varsel <- do.call(rbind, lapply(scn_res, function(x) x$est_varsel))
  tbl_iters_groupsel <- do.call(rbind, lapply(scn_res, function(x) x$est_groupsel))
  tbl_iters_allsel <- do.call(rbind, lapply(scn_res, function(x) x$est_allsel))
  tbl_iters_truecoef <- do.call(rbind, lapply(scn_res, function(x) x$true_coef)) %>% merge(scn, . )
  tbl_iters_estcoef <- do.call(rbind, lapply(scn_res, function(x) x$est_coef)) %>% merge(scn, . )
  tbl_iters_estpred <- do.call(rbind, lapply(scn_res, function(x) x$est_pred)) %>% merge(scn, . )
  niter <- max(tbl_iters_varsel$i)
  
  # -- Oracle 
  tbl_performance <- tbl_iters_performance %>%
    data.frame() %>%
    group_by(model, penalty) %>%
    summarise("RMSPE.est" = mean(RMSPE, na.rm = T), "RMSPE.med" = median(RMSPE, na.rm = T), "RMSPE.sd" = sd(RMSPE, na.rm = T),
              "RMSPE.lo" = quantile(RMSPE, 0.05, na.rm = T), "RMSPE.up" = quantile(RMSPE, 0.95, na.rm = T),
              "relRMSPE.est" = mean(relRMSPE, na.rm = T), "relRMSPE.med" = median(relRMSPE, na.rm = T), "relRMSPE.sd" = sd(relRMSPE, na.rm = T), 
              "relRMSPE.lo" = quantile(relRMSPE, 0.05, na.rm = T), "relRMSPE.up" = quantile(relRMSPE, 0.95, na.rm = T),
              "R2.est" = mean(R2, na.rm = T), "R2.med" = median(R2, na.rm = T), "R2.sd" = sd(R2, na.rm = T),
              "R2.lo" = quantile(R2, 0.05, na.rm = T), "R2.up" = quantile(R2, 0.95, na.rm = T),
              "CS.est" = mean(CS, na.rm = T), "CS.med" = median(CS, na.rm = T), "CS.sd" = sd(CS, na.rm = T),
              "CS.lo" = quantile(CS, 0.05, na.rm = T), "CS.up" = quantile(CS, 0.95, na.rm = T),
              "MAE.est" = mean(MAE, na.rm = T), "MAE.med" = median(MAE, na.rm = T), "MAE.sd" = sd(MAE, na.rm = T),
              "MAE.lo" = quantile(MAE, 0.05, na.rm = T), "MAE.up" = quantile(MAE, 0.95, na.rm = T),
              "df.est" = mean(df, na.rm = T), "df.med" = median(df, na.rm = T), "df.sd" = sd(df, na.rm = T),
              "df.lo" = quantile(df, 0.05, na.rm = T), "df.up" = quantile(df, 0.95, na.rm = T),
              "extime.est" = mean(extime, na.rm = T), "extime.med" = median(extime, na.rm = T), "extime.sd" = sd(extime, na.rm = T),
              "extime.lo" = quantile(extime, 0.05, na.rm = T), "extime.up" = quantile(extime, 0.95, na.rm = T)) %>%
    data.frame()  %>%
    merge(scn, . ) %>%
    mutate(RMSPE.std = RMSPE.est / epsstd)
  
  
  tbl_varsel <- tbl_iters_varsel %>% 
    group_by(model, penalty, varname) %>% 
    summarise("truepos_var" = sum(truepos_var)/niter, 
              "trueneg_var" = sum(trueneg_var)/niter,
              "falsepos_var" = sum(falsepos_var)/niter, 
              "falseneg_var" = sum(falseneg_var)/niter,
              "vif" = (sum(truepos_var) + sum(falsepos_var))/niter,
              "vef" =  (sum(trueneg_var) + sum(falseneg_var))/niter) %>%
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
  list_results <- list("results" = list("performance" = tbl_performance, "varsel" = tbl_varsel,  
                                        "groupsel" = tbl_groupsel, "allsel" = tbl_allsel),
                       "iters" = list("performance" = tbl_iters_performance, "varsel" = tbl_iters_varsel, 
                                      "groupsel" = tbl_iters_groupsel, "allsel" = tbl_iters_allsel),
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
  tbl_vif <- do.call(rbind, lapply(res, function(x) x$results$vif))
  
  tbl_iters_perf <- do.call(rbind, lapply(res, function(x) x$iters$performance))
  tbl_iters_varsel <- do.call(rbind, lapply(res, function(x) x$iters$varsel))
  tbl_iters_groupsel <- do.call(rbind, lapply(res, function(x) x$iters$groupsel))
  tbl_iters_allsel <- do.call(rbind, lapply(res, function(x) x$iters$allsel))
  
  tbl_estcoef <- do.call(rbind, lapply(res, function(x) x$coef$estcoef))
  tbl_estpred <- do.call(rbind, lapply(res, function(x) x$pred$estpred))
  
  # Additional analysis

  
  # Results per iteration
  tbl_res_perf <- list("tbl_perf"=tbl_perf,  "tbl_varsel"=tbl_varsel, 
                       "tbl_groupsel"=tbl_groupsel, "tbl_allsel"=tbl_allsel, 
                       "tbl_vif"=tbl_vif)
  tbl_res_iters <- list("tbl_iters_perf"=tbl_iters_perf, "tbl_iters_varsel"=tbl_iters_varsel, 
                        "tbl_iters_groupsel"=tbl_iters_groupsel, "tbl_iters_allsel"=tbl_iters_allsel)
  tbl_res_coef <- list("coef" = tbl_estcoef, "pred"=tbl_estpred)
  saveRDS(tbl_res_perf, here::here(sim.path, "tbl_results_perf.rds"))
  saveRDS(tbl_res_iters, here::here(sim.path, "tbl_results_iters.rds"))
  saveRDS(tbl_res_coef, here::here(sim.path, "tbl_results_coef.rds"))
  
  return(NULL)
}
