#' ==============================================================================
#' Author: MG
#' Date: 09/05/2023
#' Info: Functions only relevant for the Monte Carlo simualtion
#' ==============================================================================


# =========================== Run scenario =====================================

simulate_scenario <- function(scen){
  
  # Generate large validation dataset
  data.val <- data_generation(n=10000)
  data.val <- generate_dataobj(y = data.val$data_ana$y, x = data.val$data_ana$x, clinical = NULL)
  
  # Run replications
  res_all <-  lapply(1:scen$iter, function(x){
    data_iter <- data_generation(n = 200)
    res_iter <-  data_analysis(df = data_iter, data.val = data.val)
    data_iter$i <- res_iter$i <- x
    return(list("data_gen" = data_iter, "data_ana" = res_iter)) 
  })
  
  data_gen_all <- do.call(rbind, lapply(res_all, function(x) x[[1]]))
  data_ana_all <- do.call(rbind, lapply(res_all, function(x) x[[2]]))
  
  return(list("data_gen_all"=data_gen_all, "data_ana_all"=data_ana_all))
  
}


# ============================ Data generation =================================
data_generation <- function(n){
  
  # Remove later (parameter input for function)
  # n <- 200     # Number of observations
  p <- 100     # Number of variables
  groupsize <- c(25,25,25,25)
  ngroups <- 4
  rhomat <- rbind(c(.8,.2), c(.8,.2), c(.8,.2), c(.4,.2))
  power  <- 1
  x_mean <- -2.5
  coefmax <- 5
  a <- 0.5
  epsstd <- 2.5
  coreps <- .075
  prop.nonzero <- 0.5
  sampthresh <- 0.05
  
  # Groupwise Hub correlation design filled with toeplitz
  hub_cormat <- simcor.H(k = 4, size = groupsize, rho = rhomat, power = power,
                         epsilon = coreps, eidim = 2)
  hub_cormat <- nearPD(hub_cormat, base.matrix = TRUE, keepDiag = TRUE)$mat
  
  # Normally distributed data
  data_vars <- mvrnorm(n = n, mu = rep(x_mean, p), Sigma = hub_cormat)
  
  # Log normally distributed data
  data_logvars <- exp(data_vars)

  # Structural zeros
  D <- sapply(rep(rev(seq(0.3, 1, length.out=p/4)),4), function(x) rbinom(n, size=1, prob=x)) # per group linear increase in zero-inflation 0-80%
  X <- data_logvars * D
  
  # Extract hubs and indices of true predictors (scenario D)
  hubindex <- cumsum(groupsize) - groupsize + 1 # identify index of hub
  nelem <- c(5,5,5,5)  # number of true predictors in each group
  truepredindex <- c(apply(cbind(hubindex, hubindex + nelem - 1),1, function(x) seq(x[1], x[2], 1)))
  ptrue <- length(truepredindex)
  
  # Generate coeffs and error 
  beta_X <- beta_D <- rep(0, p)
  beta_X[truevarindex] <- c(sapply(nelem, function(x) rev(seq(0, coefmax, length.out = x+1)[-1])))
  beta_D[truevarindex] <-  c(sapply(nelem, function(x) sample(seq(0, coefmax, length.out = x+1)[-1], x)))
  
  # Generate outcome: a controls influence of X and D components
  eps <- rnorm(n, mean = 0, sd = epsstd)
  y <- a * X %*% beta_X + (1-a) * D %*% beta_D  + eps
  
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
               "coef" = data.frame("beta_X" = beta_X, "beta_D" = beta_D),
              "add" = list("adjR2" = adjR2))
  return(out)
}




# ============================ Data analysis =================================
data_analysis <- function(df, data.val){
  
  # Parameter (remove later)
  n = length(df$data_ana$y)
  p = ncol(df$data_ana$x)
  ncv = 10
  nR  = 2
  folds <- sample(rep(1:ncv, ceiling(n/ncv)))[1:n]
  pflist <- list(c(1,2), c(2,1), c(1,3))
  
   # Prepare data
  data.obj <- generate_dataobj(y = df$data_ana$y, x = df$data_ana$x, clinical = NULL)
  
  
  # CV 
  methnames <- c("lasso", "ridge", "lridge", "rlasso", "rgarrote")
  tbl_perf <- data.frame("methods" = rep(methnames, each = 2)[-1],
                         "penalty" = rep(c("combined", "component"), times = 5)[-2],
                         R2 = NA, RMSE = NA, MAE = NA, C = NA, CS = NA)
    
  # --- Lasso
  fit.lasso.cv <- perform_penreg(data.obj, family = "gaussian", alpha1 = 1, cv = ncv, R = nR, penalty = "combined")
  data.val$pred.lasso <- predict_penreg(obj = fit.lasso.cv, newdata = data.val, model = "lasso")
  tbl_perf[1, 4:9] <- eval_performance(pred = data.val$pred.lasso, obs = data.val$y)
  
  # --- Ridge (penalty: combined and component)
  fit.ridge.cv <- perform_penreg(data.obj, family = "gaussian", cv = ncv, R = nR, penalty = "combined")
  data.val$pred.ridge.cb <- predict_penreg(obj = fit.ridge.cv, newdata = data.val, model = "ridge")
  tbl_perf[2, 4:9] <- eval_performance(pred = data.val$pred.ridge.cb, obs = data.val$y)
  
  fit.ridge.cv <- perform_penreg(data.obj, family = "gaussian", cv = ncv, R = nR, penalty = "component", pflist = pflist)
  data.val$pred.ridge.cp <- predict_penreg(obj=fit.ridge.cv, newdata=data.val, model="ridge")
  tbl_perf[3, 4:9] <- eval_performance(pred = data.val$pred.ridge.cp, obs = data.val$y)
  
  
  # --- Lasso-ridge (penalty: combined and component)
  fit.lridge.cv <- perform_lridge(data.obj, family = "gaussian", cv = ncv, R = nR, nlambda = c(10, 10), penalty = "combined")
  data.val$pred.lridge.cb <- predict_lridge(obj = fit.lridge.cv, newdata = data.val)
  tbl_perf[4, 4:9] <- eval_performance(pred = data.val$pred.lridge.cb, obs = data.val$y)
  
  fit.lridge.cv <- perform_lridge(data.obj, family = "gaussian", cv = ncv, R = nR, nlambda = c(10, 10), penalty = "component", pflist = pflist)
  data.val$pred.lridge.cp <- predict_lridge(obj = fit.lridge.cv, newdata = data.val)
  tbl_perf[5, 4:9] <- eval_performance(pred = data.val$pred.lridge.cp, obs = data.val$y)
  
  
  # --- Ridge-lasso (penalty: combined and component)
  fit.rlasso.cv <- perform_rlasso(data.obj, family = "gaussian", cv = ncv, R = nR, nlambda = c(10, 10), penalty = "combined")
  data.val$pred.rlasso.cb <- predict_rlasso(obj = fit.rlasso.cv, newdata = data.val)
  tbl_perf[6, 4:9] <- eval_performance(pred = data.val$pred.rlasso.cb, obs = data.val$y)
  
  fit.rlasso.cv <- perform_rlasso(data.obj, family = "gaussian", cv = ncv, R = nR, nlambda = c(10, 10), penalty = "component", pflist = pflist)
  data.val$pred.rlasso.cp <- predict_rlasso(obj = fit.rlasso.cv, newdata = data.val)
  tbl_perf[7, 4:9] <- eval_performance(pred = data.val$pred.rlasso.cp, obs = data.val$y)
  
  
  # --- Ridge-garrote
  fit.rgarrote.cv <- perform_rgarrote(data.obj, family = "gaussian", cv = ncv, R = nR, nlambda = c(10,10), penalty="combined")
  data.val$pred.rgarrote.cb <- predict_rgarrote(obj = fit.rgarrote.cv, newdata = data.val)
  tbl_perf[8, 4:9] <- eval_performance(pred = data.val$pred.rgarrote.cb, obs = data.val$y)
  
  fit.rgarrote.cv <- perform_rgarrote(data.obj, family = "gaussian", cv = ncv, R = nR, nlambda = c(10,10), penalty="component", pflist = pflist)
  data.val$pred.rgarrote.cp <- predict_rgarrote(obj = fit.rgarrote.cv, newdata = data.val)
  tbl_perf[9, 4:9] <- eval_performance(pred = data.val$pred.rgarrote.cp, obs = data.val$y)

  
  # Merge results
  
  return()
}