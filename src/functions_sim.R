#' ==============================================================================
#' Author: MG
#' Date: 09/05/2023
#' Info: Functions only relevant for the Monte Carlo simulation
#' ==============================================================================


# =========================== Run scenario =====================================

simulate_scenario <- function(scn, dsgn){
  
  filename <- paste0("sim_i", scn$iter, "_scn", scn$OGM, "_n", scn$n, "_p", scn$p,"_beta", scn$beta_max, "_UDdepen", scn$UDdepen, 
                     "_eps", scn$epslvl, "_propzi", scn$propzi, "_revzi", tolower(scn$revzi), "_struczero", scn$struczero, ".rds")

  # Generate large validation dataset
  data.val <- data_generation(dsgn = dsgn, OGM = scn$OGM, n = 100000, p=scn$p, beta_max = scn$beta_max, 
                              a = scn$a, epsstd = scn$epsstd, xmean = scn$xmean, xstd = scn$xstd,
                              propzi = scn$propzi, revzi = scn$revzi, struczero = scn$struczero)
  data.val <- generate_dataobj(y = data.val$data_ana$y, z = data.val$data_ana$X, 
                               clinical = NULL, logtransform = FALSE)
  
  # Run replications in parallel
  scn_res <- future_lapply(1:scn$iter, function(x){
    data_iter <- data_generation(dsgn = dsgn, OGM = scn$OGM, n = scn$n, p = scn$p, beta_max = scn$beta_max, 
                                 a = scn$a, epsstd = scn$epsstd, xmean = scn$xmean, xstd = scn$xstd,
                                 propzi=scn$propzi, revzi = scn$revzi, struczero=scn$struczero)
    res_iter <-  data_analysis(df = data_iter, data.val = data.val, n = scn$n, p = scn$p, 
                               ncv = 10, nR = nrep, nlams = nlambdas)
    data_iter <- mapply(cbind, data_iter, "i" = x, SIMPLIFY = F)
    res_iter <- mapply(cbind, res_iter, "i" = x, SIMPLIFY = F)
    
    return(c(data_iter, res_iter)) 
  }, future.seed = seeds)
  
  # Summarize
  summarize_scenario(filename = filename, scn = scn, scn_res = scn_res)
}


# ============================ Data generation =================================
data_generation <- function(dsgn, OGM, n, p, beta_max, a, epsstd, xmean, xstd, propzi, revzi, struczero){
  
  # Parameter
  ngroups = 4
  
  # Draw data from simdesign
  data_logvars <- simulate_data(dsgn, n)

  # Structural zeros for outcome generation
  groupzi <- seq(0, propzi, length.out = p/ngroups)
  if(revzi){groupzi <- rev(groupzi)}
  seq_propzi <- rep(groupzi, ngroups)
  seq_strucz <- seq_propzi * struczero
  
  # Extract hubs and indices of true predictors
  groupsize <- rep(p/ngroups, ngroups)
  hubindex <- cumsum(groupsize) - groupsize + 1 # identify index of hub
  groupindex <- cbind(hubindex, hubindex + groupsize - 1)
  

  # Scenario specification
  if(OGM %in% "A"){
    # OGM A
    nelem <- c(groupsize[1], 0, 0, 0)
    truepredindex <- c(1:groupsize[1])
    ptrue <- length(truepredindex)    
    
    # Generate coeffs and error 
    beta_X <- beta_D <- rep(0, p)
    perm_X <- c(38, 39, 28, 31, 9, 20, 21, 6, 17, 44, 42, 48, 19, 25, 2, 49, 50, 8, 41, 15, 36, 
                18, 24, 40, 26, 10, 46, 16, 12, 22, 45, 34, 7, 47, 4, 14, 23, 37, 43, 29, 32, 30, 
                3, 1, 5, 33, 27, 11, 35, 13)
    scn_specs <- data.frame("rank_D"=50:1, "rank_X"=perm_X, "sum_rank_XD"=NA,"rank_XD"=NA)
    scn_specs$sum_rank_XD <- scn_specs$rank_D + scn_specs$rank_X
    scn_specs$rank_XD[with(scn_specs, order(scn_specs$sum_rank_XD, scn_specs$rank_D, decreasing=TRUE))] <- 1:50
    scn_specs$predind[order(scn_specs$rank_XD)] <- truepredindex[1:50]
    scn_specs$predind_grp3[order(scn_specs$rank_XD)] <- scn_specs[order(scn_specs$predind),]$predind+cumsum(groupsize)[2]
    beta_values <- rev(seq(0.1, beta_max, length.out = ptrue))
    scn_specs$beta_D[order(scn_specs$rank_D, decreasing = TRUE)] <- beta_values
    scn_specs$beta_X[order(scn_specs$rank_X, decreasing = TRUE)] <- beta_values
    beta_D[scn_specs$predind] <- scn_specs$beta_D
    beta_X[scn_specs$predind] <-  scn_specs$beta_X
    
    beta_D[rev(scn_specs$predind_grp3)] <- beta_D[scn_specs$predind]
    beta_X[rev(scn_specs$predind_grp3)] <- beta_X[scn_specs$predind]
    
    # Zero-inflation
    zi_assignment <- 1:p
    
 }else if(OGM %in% "B"){
    # OGM B
    nelem <- c(5, 5, 5, 5)  # number of true predictors in each group
    truepredindex <- floor(c(sapply(1:ngroups, function(x) seq(groupindex[x,1], groupindex[x,2], length.out=5)[1:nelem[x]])))
    ptrue <- length(truepredindex)
    
    # Generate coeffs and error 
    beta_X <- beta_D <- rep(0, p/ngroups)
    scn_specs <- data.frame( "rank_ZI"=c(3,2,4,5,1), "rank_D"=5:1, "rank_X"=c(3,1,4,2,5),
                             "sum_rank_XD"=NA,"rank_XD"=NA)
    scn_specs$sum_rank_XD <- scn_specs$rank_D + scn_specs$rank_X
    scn_specs$rank_XD[order(scn_specs$sum_rank_XD, decreasing=TRUE)] <- 1:5
    scn_specs$predind[order(scn_specs$rank_XD)] <- truepredindex[1:5]
    beta_values <- rev(seq(1, beta_max, length.out = ptrue/ngroups))
    scn_specs$beta_D[order(scn_specs$rank_D, decreasing = TRUE)] <- beta_values
    scn_specs$beta_X[order(scn_specs$rank_X, decreasing = TRUE)] <- beta_values
    scn_specs <- scn_specs[order(scn_specs$rank_ZI),]
    beta_D[scn_specs$predind] <- scn_specs$beta_D
    beta_X[scn_specs$predind] <-  scn_specs$beta_X
    beta_D <- rep(beta_D, ngroups)
    beta_X <- rep(beta_X, ngroups)
    
    # Zero-inflation specification
    zi_specs <- data.frame( "rank_ZI"=c(3,2,4,5,1), "rank_XD"=c(5,2,4,1,3))
    zi_specs$predind[order(zi_specs$rank_XD, decreasing = TRUE)]  <-   floor(seq(groupindex[1,1], groupindex[1,2], length.out=5))
    zi_specs <- zi_specs[order(zi_specs$rank_ZI),]
    vec_with_next_largest <- cbind(zi_specs$predind, sapply(zi_specs$predind, function(x) find_next_largest(x, zi_specs$predind)))
    one_group_zi_assignment <- unlist(c(sapply(1:4, function(i) zi_specs$predind[i] + 0:(apply(vec_with_next_largest, 1, diff)[i]-1)), zi_specs$predind[5]))
    zi_assignment <- c(rep(one_group_zi_assignment,2), rep(rev(one_group_zi_assignment),2)) + rep(c(0, cumsum(groupsize)[-ngroups]), times=groupsize)

  }else if(OGM %in% "C"){
    
    beta_X <- beta_D <- rep(0, p)
    beta_X[c(1:cumsum(groupsize)[1],(cumsum(groupsize)[2]+1):cumsum(groupsize)[3])] <- 0.1
    beta_D[c(1:cumsum(groupsize)[1],(cumsum(groupsize)[2]+1):cumsum(groupsize)[3])] <- 0.1
    
    strongpredindex <- floor(c(seq(groupindex[1,1], groupindex[1,2], length.out=5), seq(groupindex[3,1], groupindex[3,2], length.out=5)))
    beta_X[strongpredindex] <- beta_X[strongpredindex]*2
    beta_D[strongpredindex] <- beta_D[strongpredindex]*2
    
    negativepredindex <- floor(c(seq(groupindex[1,1], groupindex[1,2], length.out=15), seq(groupindex[3,1], groupindex[3,2], length.out=15)))
    beta_X[negativepredindex] <- beta_X[negativepredindex]*2
    beta_D[negativepredindex] <- beta_D[negativepredindex]*2
    
    # Zero-inflation
    zi_assignment <- c(1:(p/2), p:(p/2+1))
    
  }else{
    stop("OGM must be 'A','B' or 'C'")}
  
  
  # Data generation with structural zeros
  D_struc <- sapply(1-seq_strucz, function(x) rbinom(n, size = 1, prob = x), simplify = "matrix")[, zi_assignment] # per group linear increase in zero-inflation 0-80%
  X_org <- as.matrix(data_logvars * D_struc)
  struczero_X <- apply(X_org, 2, function(x) sum(x==0)/nrow(X_org))
  
  # Outcome generation: a controls influence of X and D components
  X_struc <- X_org
  X_struc[X_struc>0] <- log(X_struc[X_struc>0])
  X_struc[X_struc==0] <- xmean
  eps <- rnorm(n, mean = 0, sd = epsstd)
  y <- a * c(X_struc %*% beta_X) + (1-a) * c(D_struc %*% beta_D) + eps
  
  # Outcome modification: Inclusion of sampling zeros for data analyst
  sampzero <- 1 - struczero
  seq_sampz <- seq_propzi * sampzero
  sampz_thresh <- sapply(1:length(seq_sampz), function(i) qlnorm(seq_sampz[i], meanlog = xmean, sdlog = xstd))
  D_samp <- sapply(1:ncol(data_logvars), function(j) (data_logvars[,j] > sampz_thresh[zi_assignment[j]])*1)
  X_samp <- as.matrix((X_struc * D_struc) *D_samp)
  samzero_X <- apply(D_samp, 2, function(x) sum(x==0)/nrow(D_samp))
  zero_X <- apply(X_samp, 2, function(x) sum(x==0)/nrow(X_samp))

  
  out <- list("data_ana" = list("y" = y, "X" = X_samp),
              "data_gen" = list("y" = y, "X_org" = X_org, "X_struc" = X_struc, "data_logn" = data_logvars,
                                "D_struc" = D_struc, "D_samp" = D_samp, "eps"=eps),
              "data_coef" = data.frame("var" = paste0("V",1:length(seq_propzi)), "beta_X" = beta_X, "beta_D" = beta_D, 
                                       "propzi" = zero_X, "struczi" = struczero_X, "samplzi" = samzero_X))
  return(out)
}




# ============================ Data analysis =================================
data_analysis <- function(df, data.val, n, p, ncv, nR, nlams){

  # Parameter
  ngroups = 4
  
  # Prepare data
  data.obj <- generate_dataobj(y = df$data_ana$y, z = df$data_ana$X, 
                               clinical = NULL, logtransform = FALSE)
  true_beta <- df$data_coef
  
  # CV 
  methnames <- c("oracle-OLS", "oracle-ridge","lasso", "ridge", "lasso-ridge", "ridge-lasso", "ridge-garrote")
  tbl_perf <- data.frame("model" = methnames,
                         R2 = NA, RMSPE = NA, MAE = NA, CS = NA, npeps = NA, extime=NA)
  tbl_coef <- cbind.data.frame("var"=paste0("V", 1:nrow(true_beta)), true_beta)
  tbl_pred <- data.frame("y" = data.val$y)
  beta_true <- c(true_beta$beta_X, true_beta$beta_D)
  
  # --- Oracle
  start <- Sys.time()
  X <- cbind(data.obj$u, data.obj$d)
  true_pred <- X[,which(beta_true!=0)]
  fit.oracle <- lm(data.obj$y~true_pred)
  beta_oracle <- rep(0, ncol(X)+1)
  beta_oracle[c(1,which(beta_true!=0)+1)] <- coef(fit.oracle)
  beta_oracle[is.na(beta_oracle)] <- 0 
  tbl_pred$pred.oracle <-  beta_oracle[1] + c(cbind(data.val$u, data.val$d) %*% beta_oracle[-1])
  tbl_perf[1, 2:5] <- eval_performance(pred = tbl_pred$pred.oracle, obs = tbl_pred$y)
  tbl_perf[1, "npeps"] <- floor(sum(beta_oracle[2:(2*(p+1))]!=0, na.rm=TRUE)/2)
  tbl_coef$beta_oracle_u <- beta_oracle[2:(p+1)]
  tbl_coef$beta_oracle_d <- beta_oracle[(p+2):(2*p+1)]
  fit.oracle$coefficients <- NULL
  fit.oracle$coefficients <- beta_oracle
  attr(fit.oracle, "class") <- "oracle-OLS"
  end <- Sys.time()
  tbl_perf[1, "extime"] <- as.numeric(end-start, units="secs")
  
  # --- Oracle-ridge
  X <- cbind(data.obj$u, data.obj$d)
  true_pred_ind <- which(beta_true!=0)
  data.obj.true <- list("y" = df$data_ana$y, "x" = data.obj$x[,which(true_beta$beta_X!=0)], "u" = data.obj$u[,which(true_beta$beta_X!=0)], 
                        "d" = data.obj$d[,which(true_beta$beta_D!=0)], "clinical" = data.obj$clinical)    
  fit.roracle <- perform_penreg(data.obj.true, family = "gaussian",  alpha = 0, nl1 = nlams, cv = ncv, R = nR, split_vars = TRUE)
  tbl_pred$pred.roracle <- predict_penreg(obj = fit.roracle, newdata = data.val, model = "roracle")
  tbl_perf[2, 2:5] <- eval_performance(pred = tbl_pred$pred.roracle, obs = tbl_pred$y)
  tbl_perf[2, "npeps"] <- fit.roracle$dfvars
  tbl_perf[2, "extime"] <- fit.roracle$extime
  beta_roracle <- rep(0,ncol(X)+1)
  beta_roracle[c(1, true_pred_ind+1)] <- fit.roracle$coefficients
  beta_roracle[is.na(beta_roracle)] <- 0 
  fit.roracle$coefficients <- NULL
  fit.roracle$coefficients <- beta_roracle
  tbl_coef$beta_roracle_u <- beta_roracle[2:(p+1)]
  tbl_coef$beta_roracle_d <- beta_roracle[(p+2):(2*p+1)]
  attr(fit.roracle, "class") <- "oracle-ridge"

  # --- Lasso
  fit.lasso <- perform_penreg(data.obj, family = "gaussian", alpha1 = 1, nl1 = nlams, cv = ncv, R = nR)
  tbl_pred$pred.lasso <- predict_penreg(obj = fit.lasso, newdata = data.val, model = "lasso")
  tbl_perf[3, 2:5] <- eval_performance(pred = tbl_pred$pred.lasso, obs = tbl_pred$y)
  tbl_perf[3, "npeps"] <- fit.lasso$dfvars
  tbl_perf[3, "extime"] <- fit.lasso$extime
  tbl_coef$beta_lasso <- fit.lasso$coefficients[-1]
  
  # --- Ridge 
  fit.ridge <- perform_penreg(data.obj, family = "gaussian",  alpha = 0, nl1 = nlams, cv = ncv, R = nR, split_vars = TRUE)
  tbl_pred$pred.ridge <- predict_penreg(obj = fit.ridge, newdata = data.val, model = "ridge")
  tbl_perf[4, 2:5] <- eval_performance(pred = tbl_pred$pred.ridge, obs = tbl_pred$y)
  tbl_perf[4, "npeps"] <- fit.ridge$dfvars
  tbl_perf[4, "extime"] <- fit.ridge$extime
  tbl_coef$beta_ridge_u <- fit.ridge$coefficients[2:(p+1)]
  tbl_coef$beta_ridge_d <- fit.ridge$coefficients[(p+2):(2*p+1)]
  
  # --- Lasso-ridge 
  fit.lridge <- perform_lridge(data.obj, family = "gaussian", cv = ncv, R = nR, nlambda = rep(nlams, 2), split_vars = TRUE)
  tbl_pred$pred.lridge <- predict_lridge(obj = fit.lridge, newdata = data.val)
  tbl_perf[5, 2:5] <- eval_performance(pred = tbl_pred$pred.lridge, obs = tbl_pred$y)
  tbl_perf[5, "npeps"] <- fit.lridge$dfvars
  tbl_perf[5, "extime"] <- fit.lridge$extime
  tbl_coef$beta_lridge_u <- fit.lridge$coefficients[2:(p+1)]
  tbl_coef$beta_lridge_d <- fit.lridge$coefficients[(p+2):(2*p+1)]

  # --- Ridge-lasso 
  fit.rlasso <- perform_rlasso(data.obj, family = "gaussian", cv = ncv, R = nR, nlambda = rep(nlams, 2), split_vars = TRUE)
  tbl_pred$pred.rlasso <- predict_rlasso(obj = fit.rlasso, newdata = data.val)
  tbl_perf[6, 2:5] <- eval_performance(pred = tbl_pred$pred.rlasso, obs = tbl_pred$y)
  tbl_perf[6, "npeps"] <- fit.rlasso$dfvars
  tbl_perf[6, "extime"] <- fit.rlasso$extime
  tbl_coef$beta_rlasso_u <- fit.rlasso$coefficients[2:(p+1)]
  tbl_coef$beta_rlasso_d <- fit.rlasso$coefficients[(p+2):(2*p+1)]
  
  # --- Ridge-garrote
  fit.rgarrote <- perform_rgarrote(data.obj, family = "gaussian", cv = ncv, R = nR, nlambda = rep(nlams, 2), split_vars = TRUE)
  tbl_pred$pred.rgarrote <- predict_rgarrote(obj = fit.rgarrote, newdata = data.val)
  tbl_perf[7, 2:5] <- eval_performance(pred = tbl_pred$pred.rgarrote, obs = tbl_pred$y)
  tbl_perf[7, "npeps"] <- fit.rgarrote$dfvars
  tbl_perf[7, "extime"] <- fit.rgarrote$extime
  tbl_coef$beta_rgarrote_u <- fit.rgarrote$coefficients[2:(p+1)]
  tbl_coef$beta_rgarrote_d <- fit.rgarrote$coefficients[(p+2):(2*p+1)]
  
  # Merge results
  # Variable selection
  groupsize <- p/ngroups
  list_models <- list(fit.oracle, fit.roracle, fit.lasso, fit.ridge, 
                      fit.rlasso, fit.lridge, fit.rgarrote)
  list_sel <- list()
  for(l in 1:length(list_models)){
    fit.model <-  list_models[[l]]
    list_sel[[l]] <- eval_selection(model = attr(fit.model, "class"), varnames = tbl_coef$var,
                                    true_coef = true_beta, pred_coef = fit.model$coefficients, ngroups = ngroups, p = p)
  }
  tbl_varsel <- do.call(rbind, lapply(list_sel, function(x) x[[1]]))
  tbl_groupsel <- do.call(rbind, lapply(list_sel, function(x) x[[2]])) 
  tbl_allsel <- do.call(rbind, lapply(list_sel, function(x) x[[3]])) 
  tbl_perf[,"relRMSPE.or"] <- tbl_perf$RMSPE/tbl_perf[tbl_perf$model == "oracle-ridge",]$RMSPE
  tbl_perf[,"relRMSPE.best"] <- tbl_perf$RMSPE/min(tbl_perf$RMSPE)
  
  out <- list("est_perf" = tbl_perf, "est_coef" = tbl_coef, "est_pred" = tbl_pred,
              "est_varsel" = tbl_varsel, "est_groupsel" = tbl_groupsel, "est_allsel" = tbl_allsel)
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
  
  # -- Performance 
  tbl_performance <- tbl_iters_performance %>%
    data.frame() %>%
    group_by(model) %>%
    summarise("RMSPE.est" = mean(RMSPE, na.rm = T), "RMSPE.med" = median(RMSPE, na.rm = T), "RMSPE.sd" = sd(RMSPE, na.rm = T),
              "RMSPE.lo" = quantile(RMSPE, 0.05, na.rm = T), "RMSPE.up" = quantile(RMSPE, 0.95, na.rm = T),
              "relRMSPE.est" = mean(relRMSPE.or, na.rm = T), "relRMSPE.med" = median(relRMSPE.or, na.rm = T), "relRMSPE.sd" = sd(relRMSPE.or, na.rm = T), 
              "relRMSPE.lo" = quantile(relRMSPE.or, 0.05, na.rm = T), "relRMSPE.up" = quantile(relRMSPE.or, 0.95, na.rm = T),
              "relRMSPEbest.est" = mean(relRMSPE.best, na.rm = T), "relRMSPEbest.med" = median(relRMSPE.best, na.rm = T), "relRMSPEbest.sd" = sd(relRMSPE.best, na.rm = T), 
              "relRMSPEbest.lo" = quantile(relRMSPE.best, 0.05, na.rm = T), "relRMSPEbest.up" = quantile(relRMSPE.best, 0.95, na.rm = T),
              "R2.est" = mean(R2, na.rm = T), "R2.med" = median(R2, na.rm = T), "R2.sd" = sd(R2, na.rm = T),
              "R2.lo" = quantile(R2, 0.05, na.rm = T), "R2.up" = quantile(R2, 0.95, na.rm = T),
              "CS.est" = mean(CS, na.rm = T), "CS.med" = median(CS, na.rm = T), "CS.sd" = sd(CS, na.rm = T),
              "CS.lo" = quantile(CS, 0.05, na.rm = T), "CS.up" = quantile(CS, 0.95, na.rm = T),
              "MAE.est" = mean(MAE, na.rm = T), "MAE.med" = median(MAE, na.rm = T), "MAE.sd" = sd(MAE, na.rm = T),
              "MAE.lo" = quantile(MAE, 0.05, na.rm = T), "MAE.up" = quantile(MAE, 0.95, na.rm = T),
              "npeps.est" = mean(npeps, na.rm = T), "npeps.med" = median(npeps, na.rm = T), "npeps.sd" = sd(npeps, na.rm = T),
              "npeps.lo" = quantile(npeps, 0.05, na.rm = T), "npeps.up" = quantile(npeps, 0.95, na.rm = T),
              "extime.est" = mean(extime, na.rm = T), "extime.med" = median(extime, na.rm = T), "extime.sd" = sd(extime, na.rm = T),
              "extime.lo" = quantile(extime, 0.05, na.rm = T), "extime.up" = quantile(extime, 0.95, na.rm = T)) %>%
    data.frame()  %>%
    merge(scn, . ) %>%
    mutate(RMSPE.std = RMSPE.est / epsstd,
           RMSPE.nstd = (RMSPE.est / epsstd)*sqrt(scn$n[1]),
           MCSE.rmspe = sqrt(RMSPE.sd^2/niter))
  
  # -- Variable selection
  tbl_varsel <- tbl_iters_varsel %>% 
    group_by(model, varname) %>% 
    summarise("vif" = (sum(truepos_var) + sum(falsepos_var))/niter,
              "vef" =  (sum(trueneg_var) + sum(falseneg_var))/niter,
              "truepos_var" = sum(truepos_var)/niter, 
              "trueneg_var" = sum(trueneg_var)/niter,
              "falsepos_var" = sum(falsepos_var)/niter, 
              "falseneg_var" = sum(falseneg_var)/niter) %>%
    data.frame() %>%
    arrange(varname) %>%
    merge(scn, . )
  

  tbl_groupsel <- tbl_iters_groupsel %>% 
    group_by(model, group) %>% 
    summarise("truepos_any" = sum(truepos_any)/niter, 
              "truepos_group" = sum(truepos_group)/niter, 
              "trueneg_group" = sum(trueneg_group)/niter,
              "falsepos_group" = sum(falsepos_group)/niter, 
              "falseneg_group" = sum(falseneg_group)/niter) %>%
    data.frame()  %>%
    merge(scn, . )
  
  tbl_allsel <- tbl_iters_allsel %>% 
    group_by(model) %>% 
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
  # Preparation of simulation results for analysis
  
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

