#' ==============================================================================
#' Author: MG
#' Date: 09/05/2023
#' Info: Functions only relevant for the Monte Carlo simualtion
#' ==============================================================================


# =========================== Run scenario =====================================

simulate_scenario <- function(scen){
  
  # Run replications
  res_all <-  lapply(1:scen$iter, function(x){
    data_iter <- data_generation(n=scen$n, p=scen$p)
    res_iter <-  data_analysis(df=data_iter)
    data_iter$i <- res_iter$i <- x
    return(list("data_gen"=data_iter, "data_ana"=res_iter)) 
  })
  
  data_gen_all <- do.call(rbind, lapply(res_all, function(x) x[[1]]))
  data_ana_all <- do.call(rbind, lapply(res_all, function(x) x[[2]]))
  
  return(list("data_gen_all"=data_gen_all, "data_ana_all"=data_ana_all))
  
}


# ============================ Data generation =================================
data_generation <- function(n, p, rhomat, power, xmean, epsmean, epsstd, coreps){
  
  # Step 1: Generate Spike at Zero Variables
  # Remove later (parameter input for function)
  n <- 200     # Number of observations
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
data_analysis <- function(df){
  
   # Prepare data
  x <- df$data_ana$x
  colnames(x) <- paste0("x.", 1:ncol(x))
  d <- (x != 0)*1
  colnames(d) <- paste0("d.", 1:ncol(d))
  utmp <- apply(x, 2, function(x) ifelse(x == 0, 0, log2(x)))
  u <- apply(utmp, 2, function(x) ifelse(x == 0, mean(x[x > 0]), x))
  colnames(u) <- paste0("u.", 1:ncol(u))
  global_min <- min(x[x > 0])
  ximp <- apply(x, 2, function(x) log2(ifelse(x == 0, global_min*(1/2), x)))

  data.obj <- list("y" = df$data_ana$y, "x" = ximp, "u" = u, "d" = d, "clinical" = NULL)
  
  # lasso
  
  # ridge
  
  # ridge-garrote
  
  
  # ridge-lasso
  
  
  # lasso-ridge
  
  
  # Merge results
  
  return()
}