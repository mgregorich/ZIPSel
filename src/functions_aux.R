# ============================================================================ #
# Author: Mariella Gregorich
# Date: 02/05/2023
# Info: Functions
# ============================================================================ #


# =============================== GENERAL ======================================

max_proportionZI <- function(vec, prop_zero=0.75){
  # Input: Vector of values
  # Output: Indicate if vector falls below the non-zero threshold (True/False)
  
  out <- (sum(vec == 0)/length(vec)) <= prop_zero
  return(out)
}

find_next_largest <- function(value, vector) {
  larger_values <- vector[vector > value]
  if (length(larger_values) == 0) {
    return(min(vector))
  } else {
    return(min(larger_values))
  }
}

to_factor <- function(x){
  as.factor(as.character(x))
}

to_numeric <- function(x){
  as.numeric(as.character(x))
}

intersect_decimals <- function(x, y, digits = 10) {
  x_rounded <- round(x, digits=digits)
  y_rounded <- round(y, digits = digits)
  return(intersect(x_rounded, y_rounded))
}

partial <- function(f, ...) {
  f_args <- list(...)
  
  function(...) {
    do.call(f, c(f_args, list(...)))
  }
}


# ========================== FIGURES ============================================================

plot_calibration <- function(pred, obs, fig.title = ""){
  df <- data.frame(pred = pred, obs = obs)
  lim_lo <- plyr::round_any(min(df), 10, f = ceiling)
  lim_up <- plyr::round_any(max(df), 10, f = ceiling)
  ggplot(df, aes(x = pred, y = obs)) +
    geom_point() +
    geom_smooth(method = "gam", formula = y ~ x, col = "red")+
    geom_abline(intercept = 0, slope = 1, col = "blue") +
    scale_x_continuous("Predicted", limits = c(lim_lo, lim_up)) +
    scale_y_continuous("Observed", limits = c(lim_lo, lim_up)) +
    ggtitle(fig.title) +
    theme_bw()
}


# ========================== MODELLING ============================================================


generate_simdesign <- function(p, xmean, xstd, ngroups=4){
  

  # Groupwise Hub correlation design filled with Toeplitz
  rhomat <- rbind(c(.8,.2), c(.8,.2), c(.8,.2), c(.2,.1))  
  hub_cormat <- simcor.H(k = ngroups, size = rep(p/ngroups,ngroups),rho = rhomat, power = 1, epsilon = 0.075, eidim = 2)
  if(!matrixcalc::is.positive.definite(hub_cormat)) hub_cormat <- nearPD(hub_cormat, base.matrix = TRUE, keepDiag = TRUE)$mat
  
  # Generate simdata design
  distlist <- rep(list(partial(function(x, meanlog, sdlog) qlnorm(x, meanlog = meanlog, sdlog = sdlog), meanlog = xmean, sdlog = xstd)), nrow(hub_cormat)) 
  dsgn <- simdata::simdesign_norta(cor_target_final = hub_cormat, dist = distlist, transform_initial = data.frame,
                                   names_final = paste0("V",1:nrow(hub_cormat)), seed_initial = 1) 
  return(dsgn)
}

find_apar <- function(scn, dsgn){
    data.tmp <- data_generation(dsgn = dsgn, population = scn$population, n = scn$n, p=scn$p, beta_max = scn$beta_max, a = scn$a, epsstd = scn$epsstd, 
                                propzi = scn$propzi, revzi = scn$revzi, struczero = scn$struczero)
    
    Xb <- data.val$data_gen$X_true %*% data.val$true_coef$beta_X
    Db <- data.val$data_gen$D_struc %*% data.val$true_coef$beta_D

    a_U2D <- (4 - sqrt(8*(var(Xb) / var(Db)))) / (2*(2 - var(Xb) / var(Db)))
    a_UD <- (2 - sqrt(4*(var(Xb) / var(Db)))) / (2*(1 - var(Xb) / var(Db)))
    a <- ifelse(scn$UDdepen %in% "U=2D", a_U2D, a_UD)
    a <- ifelse(scn$UDdepen %in% "U", 1, a)
    return(a)  
}



adjusted_R2 <- function(pred, obs, N, k){
  r2 <- cor(pred,obs, use="pairwise.complete.obs")^2
  r2 <- 1-(((1-r2)*(N-1))/(N-k-1))
  return(r2)
}

c_index <- function(pred, obs){
  c.model <- concreg(data=data.frame(predicted=pred, observed=obs), observed~predicted, npar=TRUE)
  return(1-cindex(c.model))
}


eval_performance <- function(pred, obs){
  r2 <- ifelse(sd(pred)<0.001, NA, cor(pred,obs, use="pairwise.complete.obs")^2)
  CS <- lm(obs ~ pred)$coefficients[2]
 # Cind <- c_index(pred=pred, obs=obs)
  rmspe.val <- sqrt(mean((pred-obs)^2))
  mae.val <- mean(abs(pred-obs))

  res <- data.frame(R2 = r2,
                    RMSPE = rmspe.val,
                    MAE = mae.val,
                  #  C = Cind,
                    CS=CS) 
  return(res)
}

eval_selection <- function(model, penalty=NA, varnames, true_coef, pred_coef, ngroups, p){
  
  # model = attr(list_models[[3]], "class"); penalty = attr(list_models[[3]], "penalty"); varnames = tbl_coef$var;
  # true_coef = df$true_coef; pred_coef = list_models[[3]]$coefficients; ngroups = ngroups; p = p
  
  
  # Initilaization
  groupsize <- rep(p/ngroups, ngroups)
  groupindices <- cbind(cumsum(groupsize)-groupsize+1, cumsum(groupsize))
  betaD <- NA
  if(length(pred_coef[-1])!=nrow(true_coef)){betaD <- pred_coef[(p+2):(2*p+1)]}
  pred_coef <- data.frame("beta_U"=pred_coef[2:(p+1)], "beta_D"=betaD)
  
  # --- Selection per variable
  tbl_varsel <- cbind.data.frame(model, penalty, matrix(NA, nrow=nrow(true_coef), ncol=5))
  colnames(tbl_varsel) <- c("model", "penalty", "varname", "truepos_var", "trueneg_var", "falsepos_var", "falseneg_var")
  for(i in 1:nrow(true_coef)){ 
    truepos_var <- ((true_coef$beta_X[i]!=0 | true_coef$beta_D[i]!=0) & (pred_coef$beta_U[i]!=0))*1
    trueneg_var <- ((true_coef$beta_X[i]==0 & true_coef$beta_D[i]==0) & (pred_coef$beta_U[i]==0))*1
    falsepos_var <- ((true_coef$beta_X[i]==0 & true_coef$beta_D[i]==0) & (pred_coef$beta_U[i]!=0))*1
    falseneg_var <- ((true_coef$beta_X[i]!=0 & true_coef$beta_D[i]!=0) & (pred_coef$beta_U[i]==0))*1
    
    tbl_varsel[i, 3:7] <- c(varnames[i], truepos_var, trueneg_var, falsepos_var, falseneg_var)
  }
  tbl_varsel[,4:7] <- apply(tbl_varsel[,4:7], 2, to_numeric)
  
  # --- Selection per group
  tbl_groupsel <- cbind.data.frame(model, penalty, matrix(NA, nrow=length(groupsize), ncol=6))
  colnames(tbl_groupsel) <- c("model", "penalty", "group", "truepos_any", "truepos_group", "trueneg_group", "falsepos_group", "falseneg_group")
  for(j in 1:ngroups){
    true_group_eff <- any(true_coef$beta_X[groupindices[j,1]:groupindices[j,2]] !=0 | true_coef$beta_D[groupindices[j,1]:groupindices[j,2]] !=0)
    pred_group_eff_U <- any(pred_coef$beta_U[groupindices[j,1]:groupindices[j,2]] !=0)
    pred_group_eff_D <- any(pred_coef$beta_D[groupindices[j,1]:groupindices[j,2]] !=0)
    pred_group_eff <- ifelse(is.na(pred_group_eff_D), pred_group_eff_U, any(pred_group_eff_U, pred_group_eff_D))
    trueposany_group <- (true_group_eff & pred_group_eff)*1
    
    truepos_group <- sum(tbl_varsel[groupindices[j,1]:groupindices[j,2],"truepos_var"])/sum(true_coef[groupindices[j,1]:groupindices[j,2],"beta_X"]!=0)
    trueneg_group <- sum(tbl_varsel[groupindices[j,1]:groupindices[j,2],"trueneg_var"])/sum(true_coef[groupindices[j,1]:groupindices[j,2],"beta_X"]==0)
    falsepos_group <- sum(tbl_varsel[groupindices[j,1]:groupindices[j,2],"falsepos_var"])/sum(true_coef[groupindices[j,1]:groupindices[j,2],"beta_X"]==0)
    falseneg_group <- sum(tbl_varsel[groupindices[j,1]:groupindices[j,2],"falseneg_var"])/sum(true_coef[groupindices[j,1]:groupindices[j,2],"beta_X"]!=0)
    
    tbl_groupsel[j,3:8] <- c(j, trueposany_group, truepos_group, trueneg_group, falsepos_group, falseneg_group)
  }
  
  # --- Selection across all 
  tbl_allsel <- cbind.data.frame(model, penalty, matrix(NA, nrow=1, ncol=4))
  colnames(tbl_allsel) <- c("model", "penalty", "FPDR", "TPDR", "FNDR", "TNDR") 
  tbl_allsel[, "FPDR"] <- ifelse(sum(tbl_varsel$falsepos_var)==0 & (sum(tbl_varsel$truepos_var)+sum(tbl_varsel$falsepos_var)==0), 0, sum(tbl_varsel$falsepos_var)/(sum(tbl_varsel$truepos_var)+sum(tbl_varsel$falsepos_var)))
  tbl_allsel[, "TPDR"] <- ifelse(sum(tbl_varsel$truepos_var)== 0 & (sum(tbl_varsel$truepos_var)+sum(tbl_varsel$falsepos_var))==0, 0, sum(tbl_varsel$truepos_var)/(sum(tbl_varsel$truepos_var)+sum(tbl_varsel$falsepos_var)))
  tbl_allsel[, "FNDR"] <-  ifelse(sum(tbl_varsel$falseneg_var)==0 & (sum(tbl_varsel$trueneg_var)+sum(tbl_varsel$falseneg_var))==0, 0, sum(tbl_varsel$falseneg_var)/(sum(tbl_varsel$trueneg_var)+sum(tbl_varsel$falseneg_var)))
  tbl_allsel[, "TNDR"] <-  ifelse(sum(tbl_varsel$trueneg_var)==0 & (sum(tbl_varsel$trueneg_var)+sum(tbl_varsel$falseneg_var))==0, 0, sum(tbl_varsel$trueneg_var)/(sum(tbl_varsel$trueneg_var)+sum(tbl_varsel$falseneg_var)))
  
  out <- list("var_selection" = tbl_varsel, "group_selection" = tbl_groupsel, "all_selection"=tbl_allsel)
  return(out)
}

generate_dataobj <- function(y, x, clinical = NULL, logtransform = TRUE){
  colnames(x) <- paste0("x.", 1:ncol(x))
  d <- (x != 0)*1
  colnames(d) <- paste0("d.", 1:ncol(x)) 
  if(logtransform){
    u <- apply(x, 2, function(x) ifelse(x == 0, 0, log2(x)))
    meanu <- mean(u[u > 0])
    u[u==0] <- meanu
    colnames(u) <- paste0("u.", 1:ncol(x))
    global_min <- min(log2(x[x > 0]))
    minval <- ifelse(global_min<0, global_min*2, global_min*(1/2))
    ximp <- apply(x, 2, function(x) ifelse(x == 0, minval, log2(x)))   
  }else{
    u <- x
    meanu <- mean(u[u > 0])
    u[u==0] <- meanu
    colnames(u) <- paste0("u.", 1:ncol(x))
    global_min <- min(x[x > 0])
    minval <- ifelse(global_min<0, global_min*2, global_min*(1/2))
    ximp <- apply(x, 2, function(x) ifelse(x == 0, minval, x))    
  }

  data.obj <- list("y" = y, "x" = ximp, "u" = u, "d" = d, "clinical" = clinical)    
  return(data.obj)
}

generate_dataobj_reduced <- function(y, x, clinical = NULL, logtransform = TRUE, n){
  perm_ind <- sort(sample(1:nrow(x), n, replace = FALSE))
  colnames(x) <- paste0("x.", 1:ncol(x))
  d <- (x != 0)*1
  colnames(d) <- paste0("d.", 1:ncol(x)) 
  if(logtransform){
    u <- apply(x, 2, function(x) ifelse(x == 0, 0, log2(x)))
    meanu <- mean(u[u > 0])
    u[u==0] <- meanu
    colnames(u) <- paste0("u.", 1:ncol(x))
    global_min <- min(log2(x[x > 0]))
    minval <- ifelse(global_min<0, global_min*2, global_min*(1/2))
    ximp <- apply(x, 2, function(x) ifelse(x == 0, minval, log2(x)))   
  }else{
    u <- x
    meanu <- mean(u[u > 0])
    u[u==0] <- meanu
    colnames(u) <- paste0("u.", 1:ncol(x))
    global_min <- min(x[x > 0])
    minval <- ifelse(global_min<0, global_min*2, global_min*(1/2))
    ximp <- apply(x, 2, function(x) ifelse(x == 0, minval, x))    
  }
  
  data.obj <- list("y" = y[perm_ind], "x" = ximp[perm_ind,], "u" = u[perm_ind,], "d" = d[perm_ind,], "clinical" = clinical)    
  return(data.obj)
}

simcor.H <- function(k=4, size=c(10,10,10,10), 
                     rho=rbind(c(.8,.2), c(.8,.2), c(.8,.2), c(.2,.2)), power=1,
                     epsilon=.075, eidim=2){
  #' Simulating the Hub Matrix (entries filled in using Toeplitz structure)
  #' Implementation by Hardin et al. (DOI: 10.1214/13-AOAS638)
  #' k is the number of groups
  #' size is a vector of length k specifying the size of each group 
  #' rho is a vector of length k specifying base correlation values
  #' epsilon <- (1-min(rho) - 0.75*min(tau) ) - .01
  #' tau_k = (max(rho_k) - min(rho_k) )/ (size_k -2) 
  #' eidim is the space from which the noise is generated, the smaller the more noise
  #' power = 2 makes the correlations stay high
  #' power = 0.5 makes the correlations descent rapidly
  # k = k; size = groupsize; rho = rbind(c(.8,.2), c(.8,.2), c(.8,.2), c(.2,.2)); power = 1; epsilon = .075; eidim = 2
  
  ndim  <- sum(size)# dim of correlation matrix
  bigcor <- matrix(rep(0, ndim*ndim), ncol=ndim)
  
  ### generating the basic correlation matrix
  
  for (i in 1:(k) ){
    elemsize <- size[i]*(size[i]-1)/2
    corelem <-rho.func(rho[i,1],rho[i,2],power=1.5, p=elemsize) 
    cormat <- matrix(0, ncol=size[i], nrow=size[i])
    cormat[upper.tri(cormat)] <- corelem
    diag(cormat) <- 1
    cormat[lower.tri(cormat)] <- t(cormat)[lower.tri(cormat)]

    if (i==1){bigcor[1:size[1], 1:size[1]] <- cormat}
    if (i!=1){bigcor[(sum(size[1:(i-1)]) + 1):sum(size[1:i]),
                     (sum(size[1:(i-1)]) + 1):sum(size[1:i])] <- cormat}
  }
  if(!isSymmetric.matrix(bigcor)) stop("Not symmetric!")
  diag(bigcor) <- 1 - epsilon
  
  ### adding noise to the correlation matrix
  eivect <- c( )
  for (i in 1:ndim) {
    ei <- runif(eidim, -1, 1)
    eivect <- cbind(eivect, sqrt(epsilon) * ei / sqrt(sum(ei^2) ) )
  }
  
  bigE <- t(eivect) %*% eivect
  cor.nz <- bigcor + bigE
  if(!isSymmetric.matrix(cor.nz)) stop("Not symmetric!")
  
  return(cor.nz)
  }
  

rho.func <- function(r.max, r.min, power,p){
  rhovec <-c()
  rhovec[1] <- 1
  for(i in 2:(p+1)){
    rhovec[i] <- r.max - ((i-2)/(p-2))^power*(r.max-r.min)
  }
  rhovec <- rhovec[-1]
  return(rhovec)
  }


# ============================== METHODS =======================================

perform_penreg <- function(data.obj, penalties = 1, family = "gaussian", penalty = "combined", cv = 10, R = 2, nl1 = 10, alpha1 = 0, 
                           pflist = NULL, split_vars = FALSE, standardize_vars = FALSE){
  # Ridge: data.obj=data.obj; penalties=1; family="gaussian"; penalty="combined"; cv=10; R=2; nl1=10; alpha1=0; split_vars=FALSE
  # Lasso: data.obj=data.obj; penalties=1; family="gaussian"; penalty="combined"; cv=10; R=2; nl1=10; alpha1=1; pflist=NULL; split_vars=FALSE
  
  # Check for misspecifications
  start <- Sys.time()
  if(all(penalty != c("combined", "component"))){ stop("Penalty must be 'combined' or 'component'.") }
  if(is.null(pflist) & penalty == "component"){ stop("pflist must be specified if penalty='component'!") }
  if(penalty == "combined"){ pflist <- list(c(1, 1)) }
  if(all(alpha1 != c(0, 1))){ stop("alpha1 needs to be 0 or 1.") }
 
   # Preliminaries
  glmnet.control(devmax = 1, fdev = 0) # for R2=1 scenario
  x <- data.obj[["x"]]
  u <- data.obj[["u"]]
  d <- data.obj[["d"]]
  y <- data.obj[["y"]]
  if(!is.null(data.obj[["clinical"]])) clinical <- data.frame("y" = y, data.obj[["clinical"]])
  if(alpha1 == 1 & !isTRUE(split_vars)){varmat <- x
  }else if(alpha1 == 0 & isTRUE(split_vars)){
    varmat <- as.matrix(cbind(u, d))
  }else if(alpha1 == 0 & !isTRUE(split_vars)){
    varmat <- x
  }else{
    warning("Splitting X into U and D component not possible if lasso is selected. Continuing with X.")
    split_vars <- FALSE
    varmat <- x
  }
  
  n <- nrow(u)
  k <- ncol(x)
  npf <- length(pflist)
  kclin <- ifelse(!is.null(data.obj[["clinical"]]), ncol(clinical), 0)
  
  # Clinical offset
  if(is.null(data.obj[["clinical"]])){
    clin_offset_coefs <- rep(0, k)
    clin_offset <- rep(0, n)
  }else{
    fit.clin <- glm(y~., data = clinical, family = "gaussian")
    clin_offset_coefs <- fit.clin$coefficients[-1]
    clin_offset <- apply(clinical %>% dplyr::select(-y), 2, to_numeric) %*% clin_offset_coefs        
  }
  
  prederror <- lambdas_pf <- matrix(NA, nrow = npf, ncol = nl1)
  for(p in 1:npf){
    pf <- pflist[[p]]
    pfvector <- rep(pf, each=k)[1:ncol(varmat)]
    
    ## CV model with offset
    lambdas <- cvmerror <- matrix(0, R, nl1) 

    for(outer in 1:R){
      set.seed(outer)
      # CV model with offset
      cv_model <- cv.glmnet(x = varmat, y = y, alpha = alpha1, standardize = standardize_vars, 
                            lambda.min.ratio = 0.0001, offset = clin_offset, penalty.factor = pfvector,
                            nfolds = cv, nlambda = nl1)
      cvmerror[outer, ] <- cv_model$cvm
    }
    prederror[p, ] <- colSums(cvmerror) / R
    lambdas_pf[p, ] <- cv_model$lambda
  }
  index <- which(prederror == min(prederror), arr.ind = TRUE)
  index <- index[1,]
  lambda.min <- lambdas_pf[index]
  pf.min <- pflist[[index[1]]]
  pfvector.min <- rep(pf.min, each=k)[1:ncol(varmat)]
  
  ## Final model
  # CV model with offset
  fit.model <- glmnet(x = varmat, y = y, alpha = alpha1, standardize = standardize_vars, lambda.min.ratio = 0.0001,
                      offset = clin_offset, penalty.factor = pfvector.min, nfolds = cv, nlambda = nl1)
  beta <- coef(fit.model)
  coeffs <- beta[, index[2]]
  fitted.values <- cbind(1, varmat) %*% coeffs + clin_offset
  df.final <- fit.model$df[index[2]]
  dfvars <- ifelse(alpha1 == 0 & isTRUE(split_vars), df.final/2, df.final)
  
  end <- Sys.time()
  res <- list(call = match.call(), family = family, lambda = lambdas_pf, coefficients = coeffs, glmnet.fit.model = fit.model, 
              k = k, kclin = kclin, df = df.final, dfvars = dfvars, cv.pred.err = prederror, se.pred.err = 0, # se error?
              extime = as.numeric (end - start, units = "secs"), 
              fit=list(xmat = varmat, lambda = lambdas_pf, lambda.min = lambda.min, penalty = pfvector.min,
                       clin_offset = clin_offset, clin_offset_coefs = clin_offset_coefs,
                       coefficients = coeffs, beta = beta, index.lambda.min = index, 
                       fitted.values = fitted.values, split_vars = split_vars,
                       standardize_vars = standardize_vars))
  attr(res, "class") <- ifelse(alpha1 == 0, "ridge", "lasso")
  attr(res, "penalty") <- penalty
  return(res)
}

predict_penreg <- function(obj, newdata, type = "link", model = "ridge"){
  x <- newdata[["x"]]
  u <- newdata[["u"]]
  d <- newdata[["d"]]
  
  if(obj$fit$split_vars){
    varmat <- cbind(u, d)
  }else{ varmat <- x }
  
  if(!is.null(newdata[["clinical"]])){
    clinical <- newdata[["clinical"]]
    new_offset <- apply(clinical, 2, to_numeric) %*% obj$fit$clin_offset_coefs
  }else{
    new_offset <- rep(0, nrow(x))
  }
  x <- cbind(1, varmat[ ,names(obj$fit$coefficients[obj$fit$coefficients != 0])[-1]])
  beta <- obj$fit$coefficients[obj$fit$coefficients != 0]
  linpred <- x %*% beta + new_offset
  
  if(type == "response" & obj$family == "binomial"){ linpred <- plogis(linpred) }
  return(linpred)
}


perform_lridge <- function(data.obj, family = "gaussian", nlambda = c(10,10), cv = 10, R = 1, alpha1 = 1, alpha2 = 0, 
                           penalty = "combined", pflist = NULL, split_vars = FALSE, standardize_vars = FALSE){

  # data.obj = data.obj; family = "gaussian"; nlambda = rep(25,2); cv = 10; R = 1; alpha1 = 1; alpha2 = 0;
  # pflist = NULL; penalty = "combined"; split_vars = TRUE; standardize_vars = FALSE

  # Check for misspecifications
  start <- Sys.time()
  if(all(penalty != c("combined", "component"))){ stop("Penalty must be 'combined' or 'component'.") }
  if(is.null(pflist) & penalty=="component"){ stop("pflist must be specified if penalty='component'!") }
  if(!split_vars & penalty=="component"){ stop("Component-specific penalty only valid for split variable.") }
  if(penalty == "combined" & split_vars){ pflist <- list(c(1, 1))
  }else if(penalty == "combined" & !split_vars){pflist <- list(c(1))}
  
  # Preliminaries
  glmnet.control(devmax = 1, fdev = 0) # for R2=1 scenario
  x <- data.obj[["x"]]
  u <- data.obj[["u"]]
  d <- data.obj[["d"]]
  y <- data.obj[["y"]]
  if(!is.null(data.obj[["clinical"]])) clinical <- data.frame("y" = y, data.obj[["clinical"]])
  xmat <- as.matrix(x)  
  if(split_vars){
    varmat <- as.matrix(cbind(u, d))  
  }else{
    varmat <- as.matrix(x)  
  }
  
  n <- nrow(x)
  k <- ncol(x)
  kclin <- ifelse(is.null(data.obj[["clinical"]]), 0, ncol(clinical))
  npf <- length(pflist)
  
  # get number of lambdas
  nl1 <- nlambda[1]
  nl2 <- nlambda[2]

  if(k != ncol(d)) stop("d does not match x in dimensions\n")
  if(n != nrow(d) | n!= length(y)) stop("Not equal sample size in variables\n")
  
  prederror <- prederror2 <- array(NA, c(nl1, nl2, npf))
  for(p in 1:npf){ # penalty loop
    # Penalty
    pf <- pflist[[p]]

    prederr <- prederr2 <- matrix(0, nl1, nl2)
    for(outer in 1:R){ # outer loop for repetitions
      folds <- sample(rep(1:cv, ceiling(n/cv)))[1:n]
      for(inner in 1:cv){ # inner loopf for cv 
        beta <- matrix(0, ncol(varmat) + 1, nl1*nl2)
        rownames(beta) <- c("(Intercept)", colnames(varmat))
        x.train <- xmat[(1:n)[folds != inner], ]
        x.test <- xmat[(1:n)[folds == inner], ]
        var.train <- varmat[(1:n)[folds != inner], ]
        var.test <- varmat[(1:n)[folds == inner], ]
        y.train <- y[(1:n)[folds != inner]]
        y.test <- y[(1:n)[folds == inner]]
        
        # (0) Clinical offset
        if(is.null(data.obj[["clinical"]])){
          clin_offset_coefs <- rep(0, k)
          clin_offset_train <- rep(0, sum(folds!=inner))
        }else{
          c.train <- clinical[(1:n)[folds != inner], ]
          c.test <- clinical[(1:n)[folds == inner], ]
          fit.clin <- glm(y~., data = c.train, family = "gaussian")
          clin_offset_coefs <- fit.clin$coefficients[-1]
          clin_offset_train <- apply(c.train %>% dplyr::select(-y), 2, to_numeric) %*% clin_offset_coefs        
        }

        # (1) Lasso regression
        fit1.lasso <- glmnet(y = y.train, x = x.train, family = family, alpha = alpha1, nlambda = nl1, 
                             offset = clin_offset_train, standardize = standardize_vars, lambda.min.ratio = 0.0001)
        
        # (2) Ridge regression 
        for(i in 2:nl1){   # first lambda picks model with only intercept
          # Include u and d part of selected x variables
          b.lasso <- coef(fit1.lasso)[ ,i]
          nonzero.coefs <- names(b.lasso)[which(b.lasso != 0)][-1]
          nonzero.ud <- c(str_replace_all(nonzero.coefs, "x.", "u."), 
                          str_replace_all(nonzero.coefs, "x.", "d."))
          knz <- length(nonzero.coefs)   # k of non zero peptides
          
          # Ridge model with penalty
          pfvector <- rep(pf, each = knz) # pf should be c(1,1) if penalty="combined
          fit2.ridge <- glmnet(y = y.train, x = var.train[ ,nonzero.ud], alpha = alpha2, nlambda = nl2, offset = clin_offset_train, 
                               penalty.factor = pfvector, standardize = standardize_vars, lambda.min.ratio = 0.0001)
          beta[rownames(coef(fit2.ridge)), (nl2*(i-1)+1):(nl2*(i-1)+nl1)] <- as.matrix(coef(fit2.ridge))          
        } # now we have all nl1*nl2*npf betas               
        
        # validation: compute prediction error
        if(!is.null(data.obj[["clinical"]])){
          clin_offset_test <- apply(c.test %>% dplyr::select(-y), 2, to_numeric) %*% clin_offset_coefs
        }else{clin_offset_test <- rep(0, sum(folds==inner))}
        for(i in 1:nl1){
          for(ii in 1:nl2){
            yhat.test <- cbind(1, var.test) %*% beta[ , nl2*(ii-1)+i] + clin_offset_test
            if(family == "binomial") yhat.test <- plogis(yhat.test)
            prederr[i, ii] <- prederr[i, ii] + mean((y.test - yhat.test)**2)/cv/R
            prederr2[i, ii] <- prederr2[i, ii] + ((mean((y.test - yhat.test)**2))**2)/cv/R
          }
        }
      } # inner loop end
    } # outer loop end
    prederror[ , , p] <- prederr 
    prederror2[ , , p] <- prederr2 
  } # penalty loop end
  
  # Index of minimal cvm
  index <- which(prederror == min(prederror), arr.ind = TRUE)[1,]
  se.prederr <- sqrt(prederror2 - prederror^2)
  
  ## Final betas for all lambda/penalty combinations
  beta <- array(0, c(ncol(varmat) + 1, nl1 * nl2, npf))
  rownames(beta) <- c("(Intercept)", colnames(varmat))
  
  # (0) Clinical offset
  if(is.null(data.obj[["clinical"]])){
    clin_offset_coefs <- rep(0, k)
    clin_offset <- rep(0, n)
  }else{
    fit.clin <- glm(y~., data = clinical, family = "gaussian")
    clin_offset_coefs <- fit.clin$coefficients[-1]
    clin_offset <- apply(clinical %>% dplyr::select(-y), 2, to_numeric) %*% clin_offset_coefs
  }
  
  lambda <- df <- array(NA, c(nl1, nl2, npf))
  for(p in 1:npf){
    # Penalty
    pf <- pflist[[p]]

    # (1) Lasso regression
    fit1.lasso <- glmnet(y = y, x = xmat, family = family, alpha = alpha1, lambda.min.ratio = 0.0001, 
                         nlambda = nl1, offset = clin_offset, standardize = standardize_vars)
    rownames(lambda) <- fit1.lasso$lambda
    
    # (2) Ridge regression with selected variables (continuous and binary part)
    for(i in 2:nl2){
      # Include u and d part of selected vars
      b.lasso <- coef(fit1.lasso)[ ,i]
      nonzero.coefs <- names(b.lasso)[which(b.lasso != 0)][-1]
      nonzero.ud <- c(str_replace_all(nonzero.coefs, "x.", "u."), 
                      str_replace_all(nonzero.coefs, "x.", "d."))
      knz <- length(nonzero.coefs)
      
      # Ridge model with penalty
      pfvector <- rep(pf, each = knz)
      fit2.ridge <- glmnet(y = y, x = varmat[ , nonzero.ud], alpha = alpha2, nlambda = nl2, lambda.min.ratio = 0.0001, 
                                offset = clin_offset, penalty.factor = pfvector, standardize = standardize_vars)
      lambda[i, , p] <- fit2.ridge$lambda
      df[i, ,p] <- fit2.ridge$df
      beta[rownames(coef(fit2.ridge)), (nl2 * (i - 1) + 1):(nl2 * (i - 1) + nl1), p] <- as.matrix(coef(fit2.ridge))        
    } # now we have all nl1*nl2*npf beta vectors  
  }
  
  ## Return
  # lasso and ridge lambda
  lambda.min <- c(as.numeric(rownames(lambda)[index[1]]), as.numeric(lambda[index[1], index[2], index[3]]))
  names(lambda.min) <- c("lasso", "ridge")
  df.final <- df[index[1],index[2],index[3]]
  dfvars <- ifelse(split_vars, df.final/2, df.final)
  
  # Coefs for best lambda
  coeffs <- beta[ , nl2 * (index[1] - 1) + index[2], index[3]]

  # component-specific penalty factor
  pf.min <- pflist[[index[3]]]
  names(pf.min) <- "pf X"
  if(split_vars){names(pf.min) <- c("pf U", "pf D")}
  
  # fitted values
  fitted.values <- cbind(1, varmat) %*% coeffs + clin_offset
  
  end <- Sys.time()
  res <- list(call = match.call(), family = family, lambda = lambda, coefficients = coeffs, glmnet.fit1.lasso = fit1.lasso, glmnet.fit2.ridge = fit2.ridge, 
              k = k, kclin = kclin, df = df.final,  dfvars = dfvars, cv.pred.err = prederr,se.pred.err = se.prederr, 
              extime = as.numeric (end - start, units = "secs"),
              fit = list(xmat = xmat, varmat = varmat, lambda = lambda, lambda.min = lambda.min, pf.min = pf.min,
                         clin_offset = clin_offset, clin_offset_coefs = clin_offset_coefs,
                         coefficients = coeffs, beta = beta, 
                         index.lambda.min = c(index[1], index[2]), 
                         fitted.values = fitted.values, split_vars = split_vars, standardize_vars = standardize_vars)) 
  attr(res,"class") <- "lasso-ridge"
  attr(res,"penalty")<- penalty
  return(res)
}

predict_lridge <- function(obj, newdata, type = "link"){
  x <- newdata[["x"]]
  u <- newdata[["u"]]
  d <- newdata[["d"]]
  if(obj$fit$split_vars){
    varmat <- cbind(u, d)
  }else{ varmat <- x }
  
  if(!is.null(newdata[["clinical"]])){
    clinical <- data.frame(newdata[["clinical"]])
    new_offset <- apply(clinical,2,to_numeric)%*% obj$fit$clin_offset_coefs
  }else{
    new_offset <- rep(0, nrow(x))
  }
  x <-cbind(1, varmat[,names(obj$fit$coefficients[obj$fit$coefficients != 0])[-1]])
  beta <- obj$fit$coefficients[obj$fit$coefficients != 0]
  linpred <- x %*% beta + new_offset
  
  if(type == "response" & obj$family == "binomial") linpred <- plogis(linpred)
  return(linpred)
}



perform_rlasso <- function(data.obj, family = "gaussian", nlambda = c(10, 10), cv = 10, R = 2, 
                           alpha1 = 0, alpha2 = 1, penalty = "combined", pflist = NULL, split_vars = FALSE,
                           standardize_vars = FALSE){
  
  # data.obj = data.obj; nlambda=rep(10,2); cv=10; R=1; alpha1=0; alpha2=1; family="gaussian";
  # penalty="combined"; pflist=list(c(1,2), c(2,1)); split_vars = TRUE; standardize_vars=FALSE

  # Check for misspecifications
  start <- Sys.time()
  if(all(penalty != c("combined", "component"))){ stop("Penalty must be 'combined' or 'component'.") }
  if(is.null(pflist) & penalty == "component"){ stop("pflist must be specified if penalty='component'!") }
  if(!split_vars & penalty=="component"){ stop("Component-specific penalty only valid for split variable.") }
  if(penalty == "combined" & split_vars){ pflist <- list(c(1, 1))
  }else if(penalty == "combined" & !split_vars){pflist <- list(c(1))}

  # Preliminaries
  glmnet.control(devmax = 1, fdev = 0) # for R2=1 scenario
  x <- data.obj[["x"]]
  u <- data.obj[["u"]]
  d <- data.obj[["d"]]
  y <- data.obj[["y"]]  
  if(split_vars){
    varmat <- as.matrix(cbind(u, d))  
  }else{
    varmat <- as.matrix(x)  
  }
  if(!is.null(data.obj[["clinical"]])) clinical <- data.frame("y" = y, data.obj[["clinical"]])
  
  n <- nrow(x)
  k <- ncol(x)
  kclin <- ifelse(is.null(data.obj[["clinical"]]), 0, ncol(clinical))
  npf <- length(pflist)
  
  # get number of lambdas
  nl1 <- nlambda[1]
  nl2 <- nlambda[2]
  
  if(k != ncol(d)) stop("d does not match x in dimensions\n")
  if(n != nrow(d) | n != length(y)) stop("Not equal sample size in variables\n")
  
  prederror <- prederror2 <- rsquare_pf <- array(NA, c(nl1, nl2, npf))
  for(p in 1:npf){
    # Penalty
    pf <- pflist[[p]]
    pfvector <- rep(pf, each = k)
    
    prederr <- rsquare <- prederr2 <- matrix(0, nl1, nl2)
    for(outer in 1:R){
      folds <- sample(rep(1:cv, ceiling(n/cv)))[1:n]
      for(inner in 1:cv){
        beta <- matrix(0, ncol(varmat) + 1, nl1 * nl2)
        rownames(beta) <- c("(Intercept)", colnames(varmat))
        var.train <- varmat[(1:n)[folds != inner], ]
        var.test <- varmat[(1:n)[folds == inner], ]
        y.train <- y[(1:n)[folds != inner]]
        y.test <- y[(1:n)[folds == inner]]
        
        # (0) Clinical offset
        if(is.null(data.obj[["clinical"]])){
          clin_offset_coefs <- rep(0, k)
          clin_offset_train <- rep(0, sum(folds!=inner))
        }else{
          c.train <- clinical[(1:n)[folds != inner], ]
          c.test <- clinical[(1:n)[folds == inner], ]
          fit.clin <- glm(y~., data = c.train, family = "gaussian")
          clin_offset_coefs <- fit.clin$coefficients[-1]
          clin_offset_train <- apply(c.train %>% dplyr::select(-y), 2, to_numeric) %*% clin_offset_coefs        
        }
        
        # (1) Ridge regression
        fit1.ridge <- glmnet(y = y.train, x = var.train, family = family, alpha = alpha1, nlambda = nl1, offset = clin_offset_train, 
                             penalty.factor = pfvector, standardize = standardize_vars, lambda.min.ratio = 0.0001)
        
        # (2) Lasso regression 
        for(i in 1:nl1){
          # Penalty
          b.ridge <- abs(c(coef(fit1.ridge)[-1,i]))
          if(split_vars){
            penalty_x <- rep(b.ridge[1:k] + b.ridge[(k+1):(2*k)], times = 2)
          }else{penalty_x <- b.ridge}
          penalty_lasso <- 1 / c(penalty_x)^1
          
          # Lasso: lower.limits ensures positiveness of coeffs
          fit2.lasso <- glmnet(y = y.train, x = var.train, alpha = alpha2, nlambda = nl2, offset = clin_offset_train, 
                               penalty.factor = penalty_lasso, standardize = standardize_vars, lambda.min.ratio = 0.0001)
          beta[rownames(coef(fit2.lasso)), (nl2 * (i - 1) + 1):(nl2 * (i - 1) + nl1)] <- as.matrix(coef(fit2.lasso))
        } # now we have all nl1*nl2 betas               
        
        # validation: compute prediction error
        if(!is.null(data.obj[["clinical"]])){
          clin_offset_test <- apply(c.test %>% dplyr::select(-y), 2, to_numeric) %*% clin_offset_coefs
        }else{clin_offset_test <- rep(0, sum(folds==inner))}
        for(i in 1:nl1){
          for(ii in 1:nl2){
            yhat.test <- cbind(1, var.test) %*% beta[ , nl2 * (i - 1) + ii] + clin_offset_test
            if(family == "binomial") yhat.test <- plogis(yhat.test)
            prederr[i, ii] <- prederr[i, ii] + mean((y.test - yhat.test)**2)/cv/R
            prederr2[i, ii] <- prederr2[i, ii] + ((mean((y.test - yhat.test)**2))**2)/cv/R
            if(sd(yhat.test) != 0) rsquare[i, ii] <- cor(y.test, yhat.test)^2
            else rsquare[i, ii] <- 0
          }
        }
      } # inner loop end
    } # outer loop end
    prederror[,, p] <- prederr 
    prederror2[,, p] <- prederr2 
    rsquare_pf[,, p] <- rsquare
  }
  
  # Index of minimal cvm
  index <- which(prederror == min(prederror), arr.ind = TRUE)[1,] 
  se.prederror <- sqrt(prederror2 - prederror^2)
  
  ## Final model
  # (0) Clinical offset
  if(is.null(data.obj[["clinical"]])){
    clin_offset_coefs <- rep(0, k)
    clin_offset <- rep(0, n)
  }else{
    fit.clin <- glm(y~., data = clinical, family = "gaussian")
    clin_offset_coefs <- fit.clin$coefficients[-1]
    clin_offset <- apply(clinical %>% dplyr::select(-y), 2, to_numeric) %*% clin_offset_coefs
  }
  
  beta <- array(0, c(ncol(varmat) + 1, nl1 * nl2, npf))
  rownames(beta) <- c("(Intercept)", colnames(varmat))
  lambda <- df <- array(NA, c(nl1, nl2, npf))
  for(p in 1:npf){
    # Penalty
    pf <- pflist[[p]]
    pfvector <- rep(pf, each = k)
    
    # (1) Ridge regression
    fit1.ridge <- glmnet(y = y, x = varmat, family = family, alpha = alpha1, nlambda = nl1, offset = clin_offset,
                         penalty.factor = pfvector, standardize = standardize_vars, lambda.min.ratio = 0.0001)
    rownames(lambda) <- fit1.ridge$lambda
    
    # (2) Positive lasso 
    for(i in 1:nl1){
      # Penalty
      b.ridge <- abs(c(coef(fit1.ridge)[-1, i]))
      if(split_vars){
        penalty_x <- rep(b.ridge[1:k] + b.ridge[(k+1):(2*k)],2)
      }else{penalty_x <- b.ridge}
      penalty_lasso <- 1 / c(penalty_x)^1
      
      fit2.lasso <- glmnet(y = y, x = varmat, family = family, alpha = alpha2, nlambda = nl1, offset = clin_offset, 
                           penalty.factor = penalty_lasso, standardize = standardize_vars, lambda.min.ratio = 0.0001)
      lambda[i, , p] <- fit2.lasso$lambda
      df[i, 1:length(fit2.lasso$df), p] <- apply(as.matrix(coef(fit2.lasso)),2, function(x) length(unique(str_remove_all( names(which(x!=0)), "x.|d."))[-1]))
      beta[rownames(coef(fit2.lasso)),(nl2 * (i - 1) + 1):(nl2 * (i - 1) + nl1), p] <- as.matrix(coef(fit2.lasso))
    } # now we have all nl1*nl2*npf beta vectors               
  }
  ## Return
  # lasso and ridge lambda
  lambda.min <- c(as.numeric(rownames(lambda)[index[1]]), as.numeric(lambda[index[1], index[2], index[3]]))
  names(lambda.min)<-c("ridge","lasso")
  df.final <-  df[index[1],index[2],index[3]]

  # Coeffs for best lambda
  coeffs <- beta[, nl2 * (index[1] - 1) + index[2], index[3]]

  # component-specific penalty factor
  pf.min <- pflist[[index[3]]]
  names(pf.min) <- "pf X"
  if(split_vars){names(pf.min) <- c("pf U", "pf D")}
  
  # fitted values
  fitted.values <- cbind(1, varmat) %*% coeffs + clin_offset
  end <- Sys.time()
  
  res <- list(call = match.call(), family = family, lambda = lambda, coefficients = coeffs, glmnet.fit1.ridge = fit1.ridge, glmnet.fit.alasso = fit2.lasso, 
            k = k, kclin = kclin, df = df.final,  dfvars = df.final, cv.pred.err = prederror,
            cv.rsquare = rsquare_pf, se.pred.err = se.prederror, extime = as.numeric (end - start, units = "secs"),
            fit = list(xmat = varmat, lambda = lambda, lambda.min = lambda.min, pf.min = pf.min,
                     clin_offset = clin_offset, clin_offset_coefs = clin_offset_coefs,
                     coefficients = coeffs, beta = beta, # coefs ... coefs oft best lambda combination, beta.. all coefs for all lambdas
                     index.lambda.min = c(index[1], index[2]), 
                     fitted.values = fitted.values, split_vars = split_vars, standardize_vars = standardize_vars))
  attr(res,"class") <- "ridge-lasso"
  attr(res,"penalty") <- penalty

  return(res)
}

predict_rlasso <- function(obj, newdata, type = "link", split_vars = TRUE){
  x <- newdata[["x"]]
  u <- newdata[["u"]]
  d <- newdata[["d"]]
  
  if(split_vars){
    varmat <- as.matrix(cbind(u, d))  
  }else{
    varmat <- as.matrix(x)  
  }
  
  if(!is.null(newdata[["clinical"]])){
    clinical <- data.frame(newdata[["clinical"]])
    new_offset <- apply(clinical,2,to_numeric)%*% obj$fit$clin_offset_coefs
  }else{
    new_offset <- rep(0, nrow(x))
  }
  X <- cbind(1, varmat)
  beta <- obj$fit$coefficients
  linpred <- X %*% beta + new_offset
  
  if(type == "response" & obj$family == "binomial") linpred <- plogis(linpred)
  return(linpred)
}


perform_rgarrote <- function(data.obj, family = "gaussian", nlambda = c(10,10), cv = 10, R = 2, alpha1 = 0, penalty="combined", 
                             pflist=NULL, split_vars = FALSE, standardize_vars = FALSE){
  # data.obj=data.obj; family="gaussian"; nlambda=c(10,10); cv=10; R=2; alpha1=0; alpha2=1; penalty="component"; pflist=list(c(1,2), c(2,1))
  # data.obj=data.obj; family="gaussian"; nlambda=c(100,100); cv=10; R=1; alpha1=0; alpha2=1; penalty="combined";
  # pflist=NULL; split_vars = TRUE
  
  # Check for misspecifications
  start <- Sys.time()
  if(all(penalty != c("combined", "component"))){ stop("Penalty must be 'combined' or 'component'.") }
  if(is.null(pflist) & penalty == "component"){ stop("pflist must be specified if penalty='component'!") }
  if(!split_vars & penalty=="component"){ stop("Component-specific penalty only valid for split variable.") }
  if(penalty == "combined" & split_vars){ pflist <- list(c(1, 1))
  }else if(penalty == "combined" & !split_vars){pflist <- list(c(1))}
  
  # Preliminaries
  glmnet.control(devmax = 1, fdev = 0) # for R2=1 scenario
  alpha2 = 1
  x <- data.obj[["x"]]
  u <- data.obj[["u"]]
  d <- data.obj[["d"]]
  y <- data.obj[["y"]]
  if(split_vars | (alpha1 == 0 & split_vars))  varmat <- as.matrix(cbind(u, d)) else varmat <- as.matrix(x)  
  if(!is.null(data.obj[["clinical"]])) clinical <- data.frame("y" = y, data.obj[["clinical"]])
  
  n <- nrow(u)
  k <- ncol(u)
  kclin <- ifelse(is.null(data.obj[["clinical"]]), 0, ncol(clinical))
  npf <- length(pflist)
  
  # get number of lambdas
  nl1 <- nlambda[1]
  nl2 <- nlambda[2]
  
  if(k != ncol(d)) stop("d does not match x in dimensions\n")
  if(n != nrow(d) | n != length(y)) stop("Not equal sample size in variables\n")
  
  prederror <- prederror2 <- rsquare_pf <- array(NA, c(nl1, nl2, npf))
  for(p in 1:npf){
    # Penalty
    pf <- pflist[[p]]
    pfvector <- rep(pf, each = k)
    
    prederr <- rsquare <- prederr2 <- matrix(0, nl1, nl2)
    for(outer in 1:R){
      folds <- sample(rep(1:cv, ceiling(n/cv)))[1:n]
      for(inner in 1:cv){
        beta_inner <- matrix(0, ncol(varmat) + 1, nl1 * nl2)
        var.train <- varmat[(1:n)[folds != inner], ]
        var.test <- varmat[(1:n)[folds == inner], ]
        y.train <- y[(1:n)[folds != inner]]
        y.test <- y[(1:n)[folds == inner]]
        
        # (0) Clinical offset
        if(is.null(data.obj[["clinical"]])){
          clin_offset_coefs <- rep(0, k)
          clin_offset_train <- rep(0, sum(folds!=inner))
        }else{
          c.train <- clinical[(1:n)[folds != inner], ]
          c.test <- clinical[(1:n)[folds == inner], ]
          fit.clin <- glm(y~., data = c.train, family = "gaussian")
          clin_offset_coefs <- fit.clin$coefficients[-1]
          clin_offset_train <- apply(c.train %>% dplyr::select(-y), 2, to_numeric) %*% clin_offset_coefs        
        }
        
        # (1) Initial regression
        fit1.init <- glmnet(y = y.train, x = var.train, family = family, alpha = alpha1, standardize = standardize_vars,
                            nlambda = nl1, penalty.factor = pfvector, offset = clin_offset_train, lambda.min.ratio = 0.0001)
        
        # (2) Lasso regression with restriction of positive coefficients for non-negative shrinkage factors
        for(i in 2:nl1){
          # X garrote = X*beta_ridge
          beta.init <- coef(fit1.init)[, i]
          beta.init2 <- rep(beta.init[-1], each = dim(var.train)[1]) # for faster columnwise matrix-vector multiplication
          XB.ridge <-  var.train*beta.init2
          if(split_vars){varmat.gar <- XB.ridge[, seq(1, k, 1)] +  XB.ridge[, seq(k + 1, 2*k,1)]
          }else{varmat.gar <- XB.ridge[, seq(1, k, 1)]}
          
          # Positive Lasso: lower.limits ensures positiveness of coeffs
          fit2.garrote <- glmnet(y = y.train, x = varmat.gar, family = family, alpha = alpha2, standardize = standardize_vars, 
                                 lower.limits = 0, nlambda = nl2, offset = clin_offset_train, lambda.min.ratio = 0.0001)
          beta.rgarrote <- coef(fit2.garrote)
          # (3) Garrote coefficients
          for(ii in 1:nl2){
            if(split_vars){
              beta_inner[,nl2*(i-1)+ii]<-c(beta.rgarrote[1, ii], # intercept
                                           rep(beta.rgarrote[2:(k + 1), ii],2) * beta.init[2:(2 * k + 1)]) # d             
            }else{
              beta_inner[,nl2*(i-1)+ii]<-c(beta.rgarrote[1, ii], # intercept
                                           beta.rgarrote[2:(k + 1), ii] * beta.init[2:(k + 1)]) # x             
            }
          }
        } # now we have all nl1*nl2 betas               
        
        # validation: compute prediction error
        if(!is.null(data.obj[["clinical"]])){
          clin_offset_test <- apply(c.test %>% dplyr::select(-y), 2, to_numeric) %*% clin_offset_coefs
        }else{clin_offset_test <- rep(0, sum(folds==inner))}
        for(i in 1:nl1){
          for(ii in 1:nl2){
            yhat.test <- cbind(1, var.test) %*% beta_inner[, nl2*(i - 1) + ii] + clin_offset_test
            if(family == "binomial") yhat.test <- plogis(yhat.test)
            prederr[i, ii] <- prederr[i, ii] + mean((y.test - yhat.test)**2) / cv / R
            prederr2[i, ii] <- prederr2[i, ii] + ((mean((y.test - yhat.test)**2))**2) / cv / R
            if(sd(yhat.test) != 0) rsquare[i, ii] <- cor(y.test, yhat.test)^2
            else rsquare[i, ii] <- 0
          }
        }
      } # inner loop end
    } # outer loop end
    prederror[,, p] <- prederr 
    prederror2[,, p] <- prederr2 
    rsquare_pf[,, p] <- rsquare
  } 
  
  # Index of minimal cvm
  index <- which(prederror == min(prederror), arr.ind = TRUE)[1,]
  se.prederror <- sqrt(prederror2 - prederror^2)
  
  ## Final model
  # (0) Clinical offset
  if(is.null(data.obj[["clinical"]])){
    clin_offset_coefs <- rep(0, k)
    clin_offset <- rep(0, n)
  }else{
    fit.clin <- glm(y~., data = clinical, family = "gaussian")
    clin_offset_coefs <- fit.clin$coefficients[-1]
    clin_offset <- apply(clinical %>% dplyr::select(-y), 2, to_numeric) %*% clin_offset_coefs
  }
  
  beta <- array(0, c(ncol(varmat) + 1, nl1 * nl2, npf)) # all coefs for lambdas
  rownames(beta) <- c("(Intercept)", colnames(varmat))
  lambda <- df <- array(NA, c(nl1, nl2, npf))
  for(p in 1:npf){
    # Penalty
    pf <- pflist[[p]]
    pfvector <- rep(pf, each = k)
    
    # (1) Ridge regression
    fit1.init <- glmnet(y = y, x = varmat, family = family, alpha = alpha1, nlambda = nl1, standardize = standardize_vars, 
                        offset = clin_offset, penalty.factor  = pfvector, lambda.min.ratio = 0.0001)
    rownames(lambda) <- fit1.init$lambda
    
    # (2) Positive lasso 
    fit2.garrote <- vector(mode = "list", length = nl1)
    for(i in 2:nl1){
      # X garrote = X*beta_ridge
      beta.init <- coef(fit1.init)[, i]
      beta.init2 <- rep(beta.init[-1], each=dim(varmat)[1]) # for faster columnwise matrix-vector multiplication
      XB.ridge <-  varmat*beta.init2
      if(split_vars){varmat.gar <- XB.ridge[, seq(1, k, 1)] +  XB.ridge[, seq(k + 1, 2*k,1)]
      }else{varmat.gar <- XB.ridge[, seq(1, k, 1)]}
      
      fit2.garrote <- glmnet(y = y, x = varmat.gar, family = family, alpha = alpha2, lower.limits = 0, 
                             nlambda = nl2, offset = clin_offset, standardize = standardize_vars, lambda.min.ratio = 0.0001)
      beta.rgarrote <- coef(fit2.garrote)
      lambda[i, ,p] <- fit2.garrote$lambda
      df[i, ,p] <- fit2.garrote$df
      for(ii in 1:nl2){
        if(split_vars){
          beta[,nl2*(i-1)+ii, p] <- c(beta.rgarrote[1, ii], # intercept
                                      rep(beta.rgarrote[2:(k + 1), ii],2) * beta.init[2:(2 * k + 1)]) # d  
        }else{
          beta[,nl2*(i-1)+ii, p] <- c(beta.rgarrote[1, ii], # intercept
                                      beta.rgarrote[2:(k + 1), ii] * beta.init[2:(k + 1)]) # x             
        }
      }
    } # now we have all nl1*nl2 beta vectors               
  }
  
  ## Return
  # ridge and garrote lambda
  lambda.min <- c(as.numeric(rownames(lambda)[index[1]]), as.numeric(lambda[index[1], index[2], index[3]]))
  model.init <- ifelse(alpha1==0, "ridge", "lasso")
  names(lambda.min)<-c(model.init,"garrote")
  df.final <- df[index[1],index[2],index[3]]
  
  # Coefficients of best model
  coeffs <- beta[, nl2 * (index[1] - 1) + index[2], index[3]]
  
  # component-specific penalty factor
  pf.min <- pflist[[index[3]]]
  names(pf.min) <- "pf X"
  if(split_vars){names(pf.min) <- c("pf U", "pf D") }
  
  # fitted values
  fitted.values <- cbind(1, varmat) %*% beta[, nl2 * (index[1] - 1) + index[2], index[3]] + clin_offset
  
  end <- Sys.time()
  res <- list(call = match.call(), family = family, lambda = lambda, coefficients = coeffs, glmnet.fit.init = fit1.init, 
              glmnet.fit.garrote = fit2.garrote, k = k, kclin = kclin, df = df.final, dfvars = df.final, cv.pred.err = prederr,
              cv.rsquare = rsquare, se.pred.err = se.prederror, extime = as.numeric (end - start, units = "secs"), 
              fit = list(varmat = varmat, varmat.gar = varmat.gar, lambda = lambda, lambda.min = lambda.min, pf.min = pf.min,
                         clin_offset = clin_offset, clin_offset_coefs = clin_offset_coefs,
                         coefficients = coeffs, beta = beta, 
                         index.lambda.min = c(index[1], index[2], index[3]), 
                         fitted.values = fitted.values, split_vars = split_vars, standardize_vars = standardize_vars))
  attr(res, "class") <- ifelse(alpha1==0, "ridge-garrote", "lasso-garrote")
  attr(res, "penalty") <- penalty
  
  return(res)
}


predict_rgarrote <- function(obj, newdata, lambda = "lambda.min", type = "link"){
  x <- newdata[["x"]]
  u <- newdata[["u"]]
  d <- newdata[["d"]]
  if(obj$fit$split_vars){
    varmat <- cbind(u, d)
  }else{ varmat <- x }
  
  
  if(!is.null(newdata[["clinical"]])){
    clinical <- data.frame(newdata[["clinical"]])
    new_offset <- apply(clinical,2,to_numeric)%*% obj$fit$clin_offset_coefs
  }else{
    new_offset <- rep(0, nrow(varmat))
  }
  x <- cbind(1, varmat)
  beta <- obj$fit$coefficients
  linpred <- x %*% beta + new_offset

  if(type == "response" & obj$family == "binomial") linpred <- plogis(linpred)
  return(linpred)
}
