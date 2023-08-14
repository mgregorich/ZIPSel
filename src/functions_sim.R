#' ==============================================================================
#' Author: MG
#' Date: 09/05/2023
#' Info: Functions only relevant for the Monte Carlo simualtion
#' ==============================================================================


# =========================== Run scenario =====================================

simulate_scenario <- function(scn, dsgn){
  
  # Remove later
  # scn=scenarios[1,]
  # dsgn=sim_design[[scenarios[1,]$dsgn]]
  
  filename <- paste0("sim_n", scn$n, "_p", scn$p, "_beta", scn$beta_max, "_a", scn$a, 
                     "_epsstd", scn$epsstd, "_propnz", scn$prop.nonzero, "_sampthresh", scn$sampthresh, ".rds")

  # Generate large validation dataset
  data.val <- data_generation(dsgn = dsgn, n = 10000, p=scn$p, beta_max = scn$beta_max, a = scn$a, epsstd = scn$epsstd, 
                              prop.nonzero = scn$prop.nonzero, sampthresh = scn$sampthresh)
  data.val <- generate_dataobj(y = data.val$data_ana$y, x = data.val$data_ana$x, clinical = NULL)
  
  # Run replications
  res_all <-  lapply(1:scn$iter, function(x){
    data_iter <- data_generation(dsgn = dsgn, n = scn$n, p = scn$p, beta_max = scn$beta_max, a = scn$a, epsstd = scn$epsstd, 
                                 prop.nonzero=scn$prop.nonzero, sampthresh=scn$sampthresh)
    res_iter <-  data_analysis(df = data_iter, data.val = data.val, n = scn$n, p = scn$p, 
                               ncv = 10, nR = 2, nlams = 10, pflist = list(c(1,2), c(2,1), c(1,3)))
    data_iter <- mapply(cbind, data_iter, "i" = x, SIMPLIFY = F)
    res_iter <- mapply(cbind, res_iter, "i" = x, SIMPLIFY = F)
    
    return(c(data_iter, res_iter)) 
  })
  
  # Summarize
  res_dataana <- do.call(rbind, lapply(res_all, function(x) data.frame(cbind("y"=x$data_ana[1,1][[1]], x$data_ana[2,1][[1]], "i"=x$data_ana[1,2][[1]]))))
  res_datagen <- do.call(rbind, lapply(res_all, function(x) x$data_gen))
  res_truecoef <- do.call(rbind, lapply(res_all, function(x) x$true_coef))
  res_performance <- do.call(rbind, lapply(res_all, function(x) x$est_perf))
  res_estcoef <- do.call(rbind, lapply(res_all, function(x) x$est_coef))
  res_varsel <- do.call(rbind, lapply(res_all, function(x) x$est_varsel))
  res_groupsel <- do.call(rbind, lapply(res_all, function(x) x$est_groupsel))
  
  # Save results
  list_results <- list(res_datagen, res_dataana, res_truecoef, res_estcoef, res_performance, res_varsel, res_groupsel)
  saveRDS(list_results, here::here(sim.path, paste0(filename , ".rds")))  
  return(list_results)
  
}


# ============================ Data generation =================================
data_generation <- function(dsgn, n, p, beta_max, a, epsstd, prop.nonzero, sampthresh){
  
  # # Remove later (parameter input for function)
  # dsgn=dsgn; n=scn$n; p=scn$p; ngroups=scn$ngroups; 
  # prop.nonzero=scn$prop.nonzero; sampthresh=scn$sampthresh
  
  # Parameter
  xmean = 0
  xstd = 0.5
  ngroups = 4
  
  # Groupwise Hub correlation design filled with toeplitz
  data_logvars <- simulate_data(dsgn, n, seed = 2)

  # Structural zeros
  D <- as.matrix(sapply(rep(rev(seq(0.3, 1, length.out = p/4)), 4), function(x) rbinom(n, size = 1, prob = x))) # per group linear increase in zero-inflation 0-80%
  X <- as.matrix(data_logvars * D)
  
  # Extract hubs and indices of true predictors (scenario D)
  groupsize <- p/ngroups
  hubindex <- cumsum(groupsize) - groupsize + 1 # identify index of hub
  nelem <- c(5,5,5,5)  # number of true predictors in each group
  truepredindex <- c(apply(cbind(hubindex, hubindex + nelem - 1),1, function(x) seq(x[1], x[2], 1)))
  ptrue <- length(truepredindex)
  
  # Generate coeffs and error 
  beta_X <- beta_D <- rep(0, p)
  beta_X[truepredindex] <- c(sapply(nelem, function(x) rev(seq(0, beta_max, length.out = x+1)[-1])))
  beta_D[truepredindex] <-  c(sapply(nelem, function(x) sample(seq(0, beta_max, length.out = x+1)[-1], x)))
  
  # Generate outcome: a controls influence of X and D components
  eps <- rnorm(n, mean = 0, sd = epsstd)
  y <- c(a * X %*% beta_X + (1-a) * D %*% beta_D  + eps)
  
  # Sampling zeros for data analyst
  Xs <- X
  Xs[Xs < sampthresh] <- 0
  #plot(x=1:ncol(Xs), sort(apply(Xs,2, function(x) sum(x==0)/length(x))))
  
  # R2
  data <- data.frame(y, X[, truepredindex], D[, truepredindex])
  fit.lm <- lm(y~., data = data)
  sum.lm <- summary(fit.lm)
  adjR2 <- sum.lm$adj.r.squared
  
  out <- list("data_ana" = list("y" = y, "x" = Xs),
              "data_gen" = data.frame("y" = y, "X_true" = X, "D_true" = D, "eps" = eps),
              "true_coef" = data.frame("beta_X" = beta_X, "beta_D" = beta_D),
              "add" = list("adjR2" = adjR2))
  return(out)
}




# ============================ Data analysis =================================
data_analysis <- function(df, data.val, n, p,  ncv=10, nR=2, nlams=10, pflist=list(c(1,2), c(2,1), c(1,3))){
  
  # Remove later
  # df = data_iter
  # n=scn$n
  # p=scn$p
  # ngroups=scn$ngroups
  # ncv=10; nR=2; nlams=10; pflist=list(c(1,2), c(2,1), c(1,3))
  
  # Parameter
  ngroups = 4
  
  # Prepare data
  data.obj <- generate_dataobj(y = df$data_ana$y, x = df$data_ana$x, clinical = NULL)
  
  # CV 
  methnames <- c("lasso", "ridge", "lridge", "rlasso", "rgarrote", "random forest")
  tbl_perf <- data.frame("methods" = rep(methnames, each = 2)[-c(1,12)],
                         "penalty" = c(rep(c("combined", "component"), times = 5)[-2], "-"),
                         R2 = NA, RMSE = NA, MAE = NA, C = NA, CS = NA)
  tbl_coef <- df$coef
  
  # --- Lasso
  fit.lasso.cv <- perform_penreg(data.obj, family = "gaussian", alpha1 = 1, nl1 = nlams, cv = ncv, R = nR, penalty = "combined")
  data.val$pred.lasso <- predict_penreg(obj = fit.lasso.cv, newdata = data.val, model = "lasso")
  tbl_perf[1, 3:7] <- eval_performance(pred = data.val$pred.lasso, obs = data.val$y)
  tbl_coef$beta_lasso_cb <- fit.lasso.cv$coefficients[-1]
  
  # --- Ridge (penalty: combined and component)
  fit.ridge.cb <- perform_penreg(data.obj, family = "gaussian",  alpha = 0, nl1 = nlams, cv = ncv, R = nR, penalty = "combined")
  data.val$pred.ridge.cb <- predict_penreg(obj = fit.ridge.cb, newdata = data.val, model = "ridge")
  tbl_perf[2, 3:7] <- eval_performance(pred = data.val$pred.ridge.cb, obs = data.val$y)
  tbl_coef$beta_ridge_cb_u <- fit.ridge.cb$coefficients[2:(p+1)]
  tbl_coef$beta_ridge_cb_d <- fit.ridge.cb$coefficients[(p+2):(2*p+1)]
  
  fit.ridge.cp <- perform_penreg(data.obj, family = "gaussian", alpha = 0, nl1 = nlams, cv = ncv, R = nR, penalty = "component", pflist = pflist)
  data.val$pred.ridge.cp <- predict_penreg(obj=fit.ridge.cp, newdata = data.val, model = "ridge")
  tbl_perf[3, 3:7] <- eval_performance(pred = data.val$pred.ridge.cp, obs = data.val$y)
  tbl_coef$beta_ridge_cp_u <- fit.ridge.cp$coefficients[2:(p+1)]
  tbl_coef$beta_ridge_cp_d <- fit.ridge.cp$coefficients[(p+2):(2*p+1)]
  
  
  # --- Lasso-ridge (penalty: combined and component)
  fit.lridge.cb <- perform_lridge(data.obj, family = "gaussian", cv = ncv, R = nR, nlambda = rep(nlams, 2), penalty = "combined")
  data.val$pred.lridge.cb <- predict_lridge(obj = fit.lridge.cb, newdata = data.val)
  tbl_perf[4, 3:7] <- eval_performance(pred = data.val$pred.lridge.cb, obs = data.val$y)
  tbl_coef$beta_lridge_cb_u <- fit.lridge.cb$coefficients[2:(p+1)]
  tbl_coef$beta_lridge_cb_d <- fit.lridge.cb$coefficients[(p+2):(2*p+1)]
  
  fit.lridge.cp <- perform_lridge(data.obj, family = "gaussian", cv = ncv, R = nR, nlambda = rep(nlams, 2), penalty = "component", pflist = pflist)
  data.val$pred.lridge.cp <- predict_lridge(obj = fit.lridge.cp, newdata = data.val)
  tbl_perf[5, 3:7] <- eval_performance(pred = data.val$pred.lridge.cp, obs = data.val$y)
  tbl_coef$beta_lridge_cp_u <- fit.lridge.cp$coefficients[2:(p+1)]
  tbl_coef$beta_lridge_cp_d <- fit.lridge.cp$coefficients[(p+2):(2*p+1)]
  
  # --- Ridge-lasso (penalty: combined and component)
  fit.rlasso.cb <- perform_rlasso(data.obj, family = "gaussian", cv = ncv, R = nR, nlambda = rep(nlams, 2), penalty = "combined")
  data.val$pred.rlasso.cb <- predict_rlasso(obj = fit.rlasso.cb, newdata = data.val)
  tbl_perf[6, 3:7] <- eval_performance(pred = data.val$pred.rlasso.cb, obs = data.val$y)
  tbl_coef$beta_rlasso_cb_x <- fit.rlasso.cb$coefficients[-1]

  fit.rlasso.cp <- perform_rlasso(data.obj, family = "gaussian", cv = ncv, R = nR, nlambda = rep(nlams, 2), penalty = "component", pflist = pflist)
  data.val$pred.rlasso.cp <- predict_rlasso(obj = fit.rlasso.cp, newdata = data.val)
  tbl_perf[7, 3:7] <- eval_performance(pred = data.val$pred.rlasso.cp, obs = data.val$y)
  tbl_coef$beta_rlasso_cp_x <- fit.rlasso.cp$coefficients[-1]

  # --- Ridge-garrote
  fit.rgarrote.cb <- perform_rgarrote(data.obj, family = "gaussian", cv = ncv, R = nR, nlambda = rep(nlams, 2), penalty = "combined")
  data.val$pred.rgarrote.cb <- predict_rgarrote(obj = fit.rgarrote.cb, newdata = data.val)
  tbl_perf[8, 3:7] <- eval_performance(pred = data.val$pred.rgarrote.cb, obs = data.val$y)
  tbl_coef$beta_rgarrote_cb_u <- fit.rgarrote.cb$coefficients[2:(p+1)]
  tbl_coef$beta_rgarrote_cb_d <- fit.rgarrote.cb$coefficients[(p+2):(2*p+1)]
  
  fit.rgarrote.cp <- perform_rgarrote(data.obj, family = "gaussian", cv = ncv, R = nR, nlambda = rep(nlams, 2), penalty = "component", pflist = pflist)
  data.val$pred.rgarrote.cp <- predict_rgarrote(obj = fit.rgarrote.cp, newdata = data.val)
  tbl_perf[9, 3:7] <- eval_performance(pred = data.val$pred.rgarrote.cp, obs = data.val$y)
  tbl_coef$beta_rgarrote_cp_u <- fit.rgarrote.cp$coefficients[2:(p+1)]
  tbl_coef$beta_rgarrote_cp_d <- fit.rgarrote.cp$coefficients[(p+2):(2*p+1)]
  
  # --- Random forest
  train <- data.frame(y = data.obj$y, data.obj$x)
  fit.rf <- ranger(y~., data = train, num.trees = 500)
  data.val$pred.rf <- predict(fit.rf, data = data.val$x)$predictions
  tbl_perf[10, 3:7] <- eval_performance(pred = data.val$pred.rf , obs = data.val$y)
  
  # Merge results
  # Variable selection
  groupsize <- p/ngroups
  list_models <- list(fit.lasso.cv, fit.ridge.cb, fit.ridge.cp, fit.rlasso.cb, fit.rlasso.cp, fit.lridge.cb, fit.lridge.cp,
                      fit.rgarrote.cb, fit.rgarrote.cp)
  list_sel <- list()
  for(l in 1:length(list_models)){
    fit.model <-  list_models[[l]]
    list_sel[[l]] <- eval_selection(model = attr(fit.model, "class"), penalty = attr(fit.model, "penalty"), 
                                    true_coef = df$true_coef, pred_coef = fit.model$coefficients, groupsize = groupsize, p = p)
  }
  tbl_varsel <- do.call(rbind, lapply(list_sel, function(x) x[[1]]))
  tbl_groupsel <- do.call(rbind, lapply(list_sel, function(x) x[[2]])) 
  
  out <- list("est_perf" = tbl_perf, 
              "est_coef" = tbl_coef,
              "est_varsel" = tbl_varsel, 
              "est_groupsel" = tbl_groupsel)
  return(out)
}

# ============================ Evaluate all scenarios =================================

evaluate_scenarios <- function(sim.files, sim.path){
  #' Preparation of simulation results for analysis
  
  res <- lapply(sim.files, function(x) readRDS(here::here(sim.path, x)))
  tbl_perf <- do.call(rbind, lapply(res, function(x) x[[3]]))
  tbl_truecoef <- do.call(rbind, lapply(res, function(x) x[[4]]))
  tbl_estcoef <- do.call(rbind, lapply(res, function(x) x[[5]]))
  tbl_varsel <- do.call(rbind, lapply(res, function(x) x[[6]]))
  tbl_groupsel <- do.call(rbind, lapply(res, function(x) x[[7]]))
  
  out <- list("tbl_perf"=tbl_perf, "tbl_truecoef"=tbl_truecoef, "tbl_estcoef"=tbl_estcoef, 
              "tbl_varsel"=tbl_varsel, "tbl_groupsel"=tbl_groupsel)
  return(out)
}
