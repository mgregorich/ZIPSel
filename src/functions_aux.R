# ============================================================================ #
# Author: Mariella Gregorich
# Date: 02/05/2023
# Info: Functions
# ============================================================================ #


# =============================== GENERAL ======================================

find_next_largest <- function(value, vector) {
  larger_values <- vector[vector > value]
  if (length(larger_values) == 0) {
    return(min(vector))
  } else {
    return(min(larger_values))
  }
}

intersect_decimals <- function(x, y, digits = 10) {
  x_rounded <- round(x, digits=digits)
  y_rounded <- round(y, digits = digits)
  return(intersect(x_rounded, y_rounded))
}

max_proportionZI <- function(vec, prop_zero=0.75){
  # Input: Vector of values
  # Output: Indicate if vector falls below the non-zero threshold (True/False)
  
  out <- (sum(vec == 0)/length(vec)) <= prop_zero
  return(out)
}

partial <- function(f, ...) {
  f_args <- list(...)
  
  function(...) {
    do.call(f, c(f_args, list(...)))
  }
}

to_factor <- function(x){
  as.factor(as.character(x))
}

to_numeric <- function(x){
  as.numeric(as.character(x))
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
  rmspe.val <- sqrt(mean((pred-obs)^2))
  mae.val <- mean(abs(pred-obs))
  
  res <- data.frame(R2 = r2,
                    RMSPE = rmspe.val,
                    MAE = mae.val,
                    CS=CS) 
  return(res)
}

eval_selection <- function(model, varnames, true_coef, pred_coef, ngroups, p){
  
  # Initialization
  groupsize <- rep(p/ngroups, ngroups)
  groupindices <- cbind(cumsum(groupsize)-groupsize+1, cumsum(groupsize))
  betaD <- NA
  if(length(pred_coef[-1])!=nrow(true_coef)){betaD <- pred_coef[(p+2):(2*p+1)]}
  pred_coef <- data.frame("beta_U"=pred_coef[2:(p+1)], "beta_D"=betaD)
  
  # --- Selection per variable
  tbl_varsel <- cbind.data.frame(model,  matrix(NA, nrow=nrow(true_coef), ncol=5))
  colnames(tbl_varsel) <- c("model", "varname", "truepos_var", "trueneg_var", "falsepos_var", "falseneg_var")
  for(i in 1:nrow(true_coef)){ 
    truepos_var <- ((true_coef$beta_X[i]!=0 | true_coef$beta_D[i]!=0) & (pred_coef$beta_U[i]!=0))*1
    trueneg_var <- ((true_coef$beta_X[i]==0 & true_coef$beta_D[i]==0) & (pred_coef$beta_U[i]==0))*1
    falsepos_var <- ((true_coef$beta_X[i]==0 & true_coef$beta_D[i]==0) & (pred_coef$beta_U[i]!=0))*1
    falseneg_var <- ((true_coef$beta_X[i]!=0 & true_coef$beta_D[i]!=0) & (pred_coef$beta_U[i]==0))*1
    
    tbl_varsel[i, 2:6] <- c(varnames[i], truepos_var, trueneg_var, falsepos_var, falseneg_var)
  }
  tbl_varsel[,3:6] <- apply(tbl_varsel[,3:6], 2, to_numeric)
  
  # --- Selection per group
  tbl_groupsel <- cbind.data.frame(model, matrix(NA, nrow=length(groupsize), ncol=6))
  colnames(tbl_groupsel) <- c("model", "group", "truepos_any", "truepos_group", "trueneg_group", "falsepos_group", "falseneg_group")
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
    
    tbl_groupsel[j,2:7] <- c(j, trueposany_group, truepos_group, trueneg_group, falsepos_group, falseneg_group)
  }
  
  # --- Selection across all 
  tbl_allsel <- cbind.data.frame(model, matrix(NA, nrow=1, ncol=4))
  colnames(tbl_allsel) <- c("model", "FPDR", "TPDR", "FNDR", "TNDR") 
  tbl_allsel[, "FPDR"] <- ifelse(sum(tbl_varsel$falsepos_var)==0 & (sum(tbl_varsel$truepos_var)+sum(tbl_varsel$falsepos_var)==0), 0, sum(tbl_varsel$falsepos_var)/(sum(tbl_varsel$truepos_var)+sum(tbl_varsel$falsepos_var)))
  tbl_allsel[, "TPDR"] <- ifelse(sum(tbl_varsel$truepos_var)== 0 & (sum(tbl_varsel$truepos_var)+sum(tbl_varsel$falsepos_var))==0, 0, sum(tbl_varsel$truepos_var)/(sum(tbl_varsel$truepos_var)+sum(tbl_varsel$falsepos_var)))
  tbl_allsel[, "FNDR"] <-  ifelse(sum(tbl_varsel$falseneg_var)==0 & (sum(tbl_varsel$trueneg_var)+sum(tbl_varsel$falseneg_var))==0, 0, sum(tbl_varsel$falseneg_var)/(sum(tbl_varsel$trueneg_var)+sum(tbl_varsel$falseneg_var)))
  tbl_allsel[, "TNDR"] <-  ifelse(sum(tbl_varsel$trueneg_var)==0 & (sum(tbl_varsel$trueneg_var)+sum(tbl_varsel$falseneg_var))==0, 0, sum(tbl_varsel$trueneg_var)/(sum(tbl_varsel$trueneg_var)+sum(tbl_varsel$falseneg_var)))
  
  out <- list("var_selection" = tbl_varsel, "group_selection" = tbl_groupsel, "all_selection"=tbl_allsel)
  return(out)
}


find_apar <- function(scn, dsgn){
  data.tmp <- data_generation(dsgn = dsgn, OGM = scn$OGM, n = scn$n, p=scn$p, beta_max = scn$beta_max, a = scn$a, epsstd = scn$epsstd, 
                              propzi = scn$propzi, revzi = scn$revzi, struczero = scn$struczero)
  
  Xb <- data.val$data_gen$X_true %*% data.val$true_coef$beta_X
  Db <- data.val$data_gen$D_struc %*% data.val$true_coef$beta_D
  
  a_U2D <- (4 - sqrt(8*(var(Xb) / var(Db)))) / (2*(2 - var(Xb) / var(Db)))
  a_UD <- (2 - sqrt(4*(var(Xb) / var(Db)))) / (2*(1 - var(Xb) / var(Db)))
  a <- ifelse(scn$UDdepen %in% "U=2D", a_U2D, a_UD)
  a <- ifelse(scn$UDdepen %in% "U", 1, a)
  return(a)  
}

generate_dataobj <- function(y, z, clinical = NULL, logtransform = TRUE){
  d <- (z != 0)*1
  n <- nrow(z)
  colnames(d) <- paste0("d.", 1:ncol(z)) 
  if(logtransform){ # Applied example
    sdz <- apply(z, 2, function(zcol) sd(log2(zcol[zcol>0])))
    
    u <- apply(z, 2, function(zcol) ifelse(zcol == 0, 0, log2(zcol)))
    meanu <- apply(u, 2, function(ucol) mean(ucol[ucol>0]))
    u[u == 0] <- meanu[col(u)][u == 0]
    
    us <- sapply(1:ncol(u), function(j) u[,j] / sdz[j])
    colnames(us) <- paste0("u.", 1:ncol(us))
    
    x <- apply(z, 2, function(zcol) ifelse(zcol == 0, 0, log2(zcol))) 
    xs <- sapply(1:ncol(x), function(j) x[,j] / sdz[j])
    xs[xs==0] <- (1/2)* min(xs[xs>0])
    colnames(xs) <- paste0("x.", 1:ncol(xs))
    
    data.obj <- list("y" = y, "x" = xs, "u" = us, "d" = d, "clinical" = clinical)    
  }else{ # Simulation study
    u <- z
    meanu <- mean(u[u > 0])
    u[u==0] <- meanu
    colnames(u) <- paste0("u.", 1:ncol(z))
    global_min <- min(z[z > 0])
    minval <- ifelse(global_min<0, global_min*2, global_min*(1/2))
    
    x <- z
    x[x==0] <- minval
    colnames(x) <- paste0("x.", 1:ncol(x))
    
    data.obj <- list("y" = y, "x" = x, "u" = u, "d" = d, "clinical" = clinical)    
  }
  
  return(data.obj)
}


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

rho.func <- function(r.max, r.min, power,p){
  #' Implementation by Hardin et al. (DOI: 10.1214/13-AOAS638)
  rhovec <-c()
  rhovec[1] <- 1
  for(i in 2:(p+1)){
    rhovec[i] <- r.max - ((i-2)/(p-2))^power*(r.max-r.min)
  }
  rhovec <- rhovec[-1]
  return(rhovec)
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
  

# ============================== METHODS =======================================

perform_penreg <- function(data.obj, family = "gaussian", cv = 10, R = 10, nl1 = 100, alpha1 = 0, 
                           split_vars = FALSE, standardize_vars = FALSE, lmin.ratio = 0.0001){

  # Check for misspecifications
  start <- Sys.time()
  if(all(alpha1 != c(0, 1))){ stop("alpha1 needs to be 0 or 1.") }
 
   # Preliminaries
  glmnet.control(devmax = 0.999, fdev = 1e-05) # for R2=1 scenario
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
  
  ## CV model with offset
  cvm_error <- matrix(0, ncol=nl1, nrow=R) 
  cvm_stderr <- matrix(0, ncol=nl1, nrow=R) 
  lambdas <- NULL
  for(outer in 1:R){
    set.seed(outer)
    # CV model with offset
    cv_model <- cv.glmnet(x = varmat, y = y, alpha = alpha1, standardize = standardize_vars, 
                          offset = clin_offset, lambda.min.ratio = lmin.ratio,
                          nfolds = cv, nlambda = nl1)
    cvm_error[outer,] <-  c(cv_model$cvm, rep(cv_model$cvm[length(cv_model$cvm)], nl1-length(cv_model$cvm))) 
    cvm_stderr[outer,] <-  c(cv_model$cvsd, rep(cv_model$cvsd[length(cv_model$cvsd)], nl1-length(cv_model$cvsd))) 
    
    lambdas <- union(lambdas, cv_model$lambda)
  }
  lambda_seq <- rep(NA, nl1)
  lambda_seq[1:length(lambdas)] <- lambdas
  colnames(cvm_error) <- colnames(cvm_stderr) <- lambda_seq
  cvmerror <- colSums(cvm_error)
  index <- which(cvmerror == min(cvmerror), arr.ind = TRUE)
  index <- index[1]
  lambda.min <- lambdas[index]

  ## Final model
  # CV model with offset
  fit.model <- glmnet(x = varmat, y = y, alpha = alpha1, standardize = standardize_vars, 
                      offset = clin_offset, nfolds = cv, nlambda = nl1, lambda.min.ratio=lmin.ratio)
  beta <- coef(fit.model)
  coeffs <- beta[, index]
  fitted.values <- cbind(1, varmat) %*% coeffs + clin_offset
  df.final <- fit.model$df[index[2]]
  dfvars <- length(unique(str_remove_all(names(which(coeffs!=0))[-1], "x.|u.|d.")))
  
  end <- Sys.time()
  res <- list(call = match.call(), family = family, nlambda = nl1, coefficients = coeffs, glmnet.fit.model = fit.model, 
              k = k, kclin = kclin, df = df.final, dfvars = dfvars, cv.pred.err = cvm_error, cv.pred.sd = cvm_stderr,
              split_vars = split_vars, standardize_vars = standardize_vars,
              extime = as.numeric (end - start, units = "secs"), 
              fit=list(xmat = varmat, lambda = lambdas, lambda.min = lambda.min, 
                       clin_offset = clin_offset, clin_offset_coefs = clin_offset_coefs,
                       coefficients = coeffs, beta = beta, index.lambda.min = index, 
                       fitted.values = fitted.values))
  attr(res, "class") <- ifelse(alpha1 == 0, "ridge", "lasso")
  return(res)
}

predict_penreg <- function(obj, newdata, type = "link", model = "ridge"){
  x <- newdata[["x"]]
  u <- newdata[["u"]]
  d <- newdata[["d"]]
  
  if(obj$split_vars){
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


perform_lridge <- function(data.obj, family = "gaussian", nlambda = c(100,100), cv = 10, R = 10,  
                           split_vars = FALSE, standardize_vars = FALSE, lmin.ratio = 0.0001){
  
  # Preliminaries
  start <- Sys.time()
  glmnet.control(devmax = 1, fdev = 0) # for R2=1 scenario
  x <- data.obj[["x"]]
  u <- data.obj[["u"]]
  d <- data.obj[["d"]]
  y <- data.obj[["y"]]
  if(!is.null(data.obj[["clinical"]])) clinical <- data.frame("y" = y, data.obj[["clinical"]])
  xmat <- as.matrix(x)  
  if(split_vars) varmat <- as.matrix(cbind(u, d)) else varmat <- as.matrix(x)  
  
  n <- nrow(x)
  k <- ncol(x)
  kclin <- ifelse(is.null(data.obj[["clinical"]]), 0, ncol(clinical))

  # get number of lambdas
  nl1 <- nlambda[1]
  nl2 <- nlambda[2]
  
  if(k != ncol(d)) stop("d does not match x in dimensions\n")
  if(n != nrow(d) | n!= length(y)) stop("Not equal sample size in variables\n")
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
      fit1.lasso <- glmnet(y = y.train, x = x.train, family = family, alpha = 1, standardize = standardize_vars,
                           nlambda = nl1, offset = clin_offset_train, lambda.min.ratio=lmin.ratio)
      
      # (2) Ridge regression 
      for(i in 2:nl1){   # first lambda picks model with only intercept
        # Include u and d part of selected x variables
        b.lasso <- coef(fit1.lasso)[ ,i]
        nonzero.coefs <- names(b.lasso)[which(b.lasso != 0)][-1]
        knz <- length(nonzero.coefs)   # k of non zero peptides
        if(split_vars){nonzero.coefs <- c(str_replace_all(nonzero.coefs, "x.", "u."), str_replace_all(nonzero.coefs, "x.", "d."))          }
        
        # Ridge model with penalty
        if(length(nonzero.coefs)>1){
          fit2.ridge <- glmnet(y = y.train, x = var.train[ ,nonzero.coefs], alpha = 0, nlambda = nl2, offset = clin_offset_train,
                               standardize = standardize_vars, lambda.min.ratio=lmin.ratio)
          beta[rownames(coef(fit2.ridge)), (nl2*(i-1)+1):(nl2*(i-1)+nl1)] <- as.matrix(coef(fit2.ridge))          
        }} # now we have all nl1*nl2*npf betas               
      
      # validation: compute prediction error
      if(!is.null(data.obj[["clinical"]])){
        clin_offset_test <- apply(c.test %>% dplyr::select(-y), 2, to_numeric) %*% clin_offset_coefs
      }else{clin_offset_test <- rep(0, sum(folds==inner))}
      yhat.test <- apply(beta, 2, function(bcol) cbind(1, var.test) %*% bcol + clin_offset_test)
      if(family == "binomial") yhat.test <- plogis(yhat.test)
      squared_diff <- (y.test - yhat.test)^2
      prederr <- prederr + matrix(colMeans(squared_diff), ncol = nl2, byrow = TRUE) / cv / R
      prederr2 <- prederr2 + matrix((colMeans(squared_diff))^2, ncol = nl2, byrow = TRUE) / cv / R
    } # inner loop end
  } # outer loop end

  
  # Index of minimal cvm
  index <- which(prederr == min(prederr), arr.ind = TRUE)[1,]
  se.prederr <- sqrt(prederr2 - prederr^2)
  
  ## Final betas for all lambda/penalty combinations
  beta <- matrix(0, nrow=ncol(varmat) + 1, ncol=nl1 * nl2)
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
  
  lambda <- df <- matrix(NA, nrow=nl1, ncol=nl2)
  # (1) Lasso regression
  fit1.lasso <- glmnet(y = y, x = xmat, family = family, alpha = 1, 
                       nlambda = nl1, offset = clin_offset, standardize = standardize_vars, lambda.min.ratio=lmin.ratio)
  rownames(lambda) <- fit1.lasso$lambda
  
  # (2) Ridge regression with selected variables (continuous and binary part)
  for(i in 2:nl1){
    # Include u and d part of selected vars
    b.lasso <- coef(fit1.lasso)[ ,i]
    nonzero.coefs <- names(b.lasso)[which(b.lasso != 0)][-1]
    knz <- length(nonzero.coefs)   # k of non zero peptides
    if(split_vars){nonzero.coefs <- c(str_replace_all(nonzero.coefs, "x.", "u."), str_replace_all(nonzero.coefs, "x.", "d."))          }
    
    # Ridge model with penalty
    if(length(nonzero.coefs)>1){
    fit2.ridge <- glmnet(y = y, x = varmat[ , nonzero.coefs], alpha = 0, nlambda = nl2, offset = clin_offset, 
                         standardize = standardize_vars, lambda.min.ratio=lmin.ratio)
    lambda[i, ] <- fit2.ridge$lambda
    df[i, ] <- fit2.ridge$df
    beta[rownames(coef(fit2.ridge)), (nl2 * (i - 1) + 1):(nl2 * (i - 1) + nl1) ] <- as.matrix(coef(fit2.ridge))          
  }}  
  rownames(beta) <- ifelse(rep(isTRUE(split_vars),dim(varmat)[2]+1), c("(Intercept)", colnames(u), colnames(d)), c("(Intercept)", colnames(x)))

  ## Return
  # lasso and ridge lambda
  lambda.min <- c(as.numeric(rownames(lambda)[index[1]]), as.numeric(lambda[index[1], index[2]]))
  names(lambda.min) <- c("lasso", "ridge")
  df.final <- df[index[1],index[2]]

  # Coefs for best lambda
  coeffs <- beta[ , nl2 * (index[1] - 1) + index[2]]
  names(coeffs) <- ifelse(rep(isTRUE(split_vars),dim(varmat)[2]+1), c("(Intercept)", colnames(u), colnames(d)), c("(Intercept)", colnames(x)))
  dfvars <- length(unique(str_remove_all(names(which(coeffs!=0))[-1], "x.|u.|d.")))
  
  # fitted values
  fitted.values <- cbind(1, varmat) %*% coeffs + clin_offset
  
  end <- Sys.time()
  res <- list(call = match.call(), family = family, lambda = lambda, coefficients = coeffs, glmnet.fit1.lasso = fit1.lasso, glmnet.fit2.ridge = fit2.ridge, 
              k = k, kclin = kclin, df = df.final,  dfvars = dfvars, cv.pred.err = prederr, se.pred.err = se.prederr, 
              extime = as.numeric (end - start, units = "secs"),
              fit = list(xmat = xmat, varmat = varmat, lambda = lambda, lambda.min = lambda.min,
                         clin_offset = clin_offset, clin_offset_coefs = clin_offset_coefs,
                         coefficients = coeffs, beta = beta, 
                         index.lambda.min = index, 
                         fitted.values = fitted.values, split_vars = split_vars)) 
  attr(res,"class") <- "lasso-ridge"
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

perform_rlasso <- function(data.obj, family = "gaussian", nlambda = c(100, 100), cv = 10, R = 10, 
                           split_vars = FALSE, standardize_vars = FALSE, lmin.ratio = 0.000001){
  
  # Preliminaries
  start <- Sys.time()
  glmnet.control(devmax = 1, fdev = 0) # for R2=1 scenario
  x <- data.obj[["x"]]
  u <- data.obj[["u"]]
  d <- data.obj[["d"]]
  y <- data.obj[["y"]]  
  if(split_vars) varmat <- as.matrix(cbind(u, d)) else varmat <- as.matrix(x)  
  if(!is.null(data.obj[["clinical"]])) clinical <- data.frame("y" = y, data.obj[["clinical"]])
  
  n <- nrow(x)
  k <- ncol(x)
  kclin <- ifelse(is.null(data.obj[["clinical"]]), 0, ncol(clinical))

  # get number of lambdas
  nl1 <- nlambda[1]
  nl2 <- nlambda[2]
  
  if(k != ncol(d)) stop("d does not match x in dimensions\n")
  if(n != nrow(d) | n != length(y)) stop("Not equal sample size in variables\n")
  
  prederr <- prederr2 <- matrix(0, nl1, nl2)
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
      fit1.ridge <- glmnet(y = y.train, x = var.train, family = family, alpha = 0, nlambda = nl1, offset = clin_offset_train, 
                           standardize = standardize_vars, lambda.min.ratio=lmin.ratio)
      
      # (2) Lasso regression 
      coef.ridge <- abs(coef(fit1.ridge))[-1,]
      if(split_vars) coef.ridge <- rbind(coef.ridge[1:k,] + coef.ridge[(k+1):(2*k),], coef.ridge[1:k,] + coef.ridge[(k+1):(2*k),])
      penalty_mat <- 1/coef.ridge
      for(i in 1:nl1){
        penalty_lasso <- penalty_mat[,i]
        fit2.lasso <- glmnet(y = y.train, x = var.train, alpha = 1, nlambda = nl2, offset = clin_offset_train, 
                             standardize = standardize_vars, penalty.factor = penalty_lasso, lambda.min.ratio=lmin.ratio)
        beta[rownames(coef(fit2.lasso)), (nl2 * (i - 1) + 1):(nl2 * (i - 1) + nl1)] <- as.matrix(coef(fit2.lasso))
      } # now we have all nl1*nl2 betas               
      
      # validation: compute prediction error
      if(!is.null(data.obj[["clinical"]])){
        clin_offset_test <- apply(c.test %>% dplyr::select(-y), 2, to_numeric) %*% clin_offset_coefs
      }else{clin_offset_test <- rep(0, sum(folds==inner))}
      yhat.test <- apply(beta, 2, function(bcol) cbind(1, var.test) %*% bcol + clin_offset_test)
      if(family == "binomial") yhat.test <- plogis(yhat.test)
      squared_diff <- (y.test - yhat.test)^2
      prederr <- prederr + matrix(colMeans(squared_diff), ncol = nl2, byrow = TRUE) / cv / R
      prederr2 <- prederr2 + matrix((colMeans(squared_diff))^2, ncol = nl2, byrow = TRUE) / cv / R
    } # inner loop end
  } # outer loop end
  
  # Index of minimal cvm
  index <- which(prederr == min(prederr), arr.ind = TRUE)[1,] 
  se.prederror <- sqrt(prederr2 - prederr^2)
  
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
  
  beta <- matrix(0, nrow=ncol(varmat) + 1, ncol=nl1 * nl2)
  rownames(beta) <- c("(Intercept)", colnames(varmat))
  lambda <- df <- matrix(NA, ncol=nl1, nrow=nl2)

  # (1) Ridge regression
  fit1.ridge <- glmnet(y = y, x = varmat, family = family, alpha = 0, nlambda = nl1, offset = clin_offset,
                       standardize = standardize_vars, lambda.min.ratio=lmin.ratio)
  rownames(lambda) <- fit1.ridge$lambda
  
  # (2) Lasso 
  coef.ridge <- abs(as.matrix(coef(fit1.ridge)))[-1,]
  if(split_vars) coef.ridge <- rbind(coef.ridge[1:k,] + coef.ridge[(k+1):(2*k),], coef.ridge[1:k,] + coef.ridge[(k+1):(2*k),])
  penalty_mat <- 1/coef.ridge
  for(i in 1:nl1){
    penalty_lasso <- penalty_mat[,i]
    fit2.lasso <- glmnet(y = y, x = varmat, family = family, alpha = 1, nlambda = nl1, offset = clin_offset, 
                         penalty.factor = penalty_lasso, standardize = standardize_vars, lambda.min.ratio=lmin.ratio)
    lambda[i, ] <- fit2.lasso$lambda
    df[i, ] <- apply(as.matrix(coef(fit2.lasso)),2, function(x) length(unique(str_remove_all( names(which(x!=0)), "u.|x.|d."))[-1]))
    beta[rownames(coef(fit2.lasso)),(nl2 * (i - 1) + 1):(nl2 * (i - 1) + nl1)] <- as.matrix(coef(fit2.lasso))
  } # now we have all nl1*nl2 beta vectors               
  
  ## Return
  # lasso and ridge lambda
  lambda.min <- c(as.numeric(rownames(lambda)[index[1]]), as.numeric(lambda[index[1], index[2]]))
  names(lambda.min) <- c("ridge","lasso")
  df.final <-  df[index[1],index[2]]

  # Coeffs for best lambda
  coeffs <- beta[, nl2 * (index[1] - 1) + index[2]]
  dfvars <- length(unique(str_remove_all(names(which(coeffs!=0))[-1], "x.|u.|d.")))
  
  # fitted values
  fitted.values <- cbind(1, varmat) %*% coeffs + clin_offset
  end <- Sys.time()

  res <- list(call = match.call(), family = family, lambda = lambda, coefficients = coeffs, glmnet.fit1.ridge = fit1.ridge, glmnet.fit.alasso = fit2.lasso, 
            k = k, kclin = kclin, df = df.final,  dfvars = dfvars, cv.pred.err = prederr,
            se.pred.err = se.prederror, extime = as.numeric (end - start, units = "secs"),
            fit = list(xmat = varmat, lambda = lambda, lambda.min = lambda.min, 
                     clin_offset = clin_offset, clin_offset_coefs = clin_offset_coefs,
                     coefficients = coeffs, beta = beta, 
                     index.lambda.min = c(index[1], index[2]), 
                     fitted.values = fitted.values, split_vars = split_vars, standardize_vars = standardize_vars))
  attr(res,"class") <- "ridge-lasso"
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

perform_rgarrote <- function(data.obj, family = "gaussian", nlambda = c(100, 100), cv = 10, R = 10, split_vars = FALSE, 
                             standardize_vars = FALSE, lmin.ratio = 0.000001){

  # Preliminaries
  start <- Sys.time()
  glmnet.control(devmax = 1, fdev = 0) # for R2=1 scenario
  x <- data.obj[["x"]]
  u <- data.obj[["u"]]
  d <- data.obj[["d"]]
  y <- data.obj[["y"]]
  if(split_vars)  varmat <- as.matrix(cbind(u, d)) else varmat <- as.matrix(x)  
  if(!is.null(data.obj[["clinical"]])) clinical <- data.frame("y" = y, data.obj[["clinical"]])
  
  n <- nrow(u)
  k <- ncol(u)
  kclin <- ifelse(is.null(data.obj[["clinical"]]), 0, ncol(clinical))

  # get number of lambdas
  nl1 <- nlambda[1]
  nl2 <- nlambda[2]
  
  if(k != ncol(d)) stop("d does not match x in dimensions\n")
  if(n != nrow(d) | n != length(y)) stop("Not equal sample size in variables\n")
  
  prederr <- prederr2 <- matrix(0, nl1, nl2)
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
      fit1.init <- glmnet(y = y.train, x = var.train, family = family, alpha = 0, standardize = standardize_vars,
                          nlambda = nl1, offset = clin_offset_train, lambda.min.ratio=lmin.ratio)
      
      # (2) Lasso regression with restriction of positive coefficients for non-negative shrinkage factors
      for(i in 2:nl1){
        # X garrote = X*beta_ridge
        beta.init <- coef(fit1.init)[, i]
        XB.ridge <-  var.train*rep(beta.init[-1], each = dim(var.train)[1]) 
        if(split_vars){varmat.gar <- XB.ridge[, seq(1, k, 1)] +  XB.ridge[, seq(k + 1, 2*k,1)]
        }else{varmat.gar <- XB.ridge[, seq(1, k, 1)]}
        
        # Positive Lasso: lower.limits ensures positiveness of coeffs
        fit2.garrote <- glmnet(y = y.train, x = varmat.gar, family = family, alpha = 1, standardize = standardize_vars, 
                               lower.limits = 0, nlambda = nl2, offset = clin_offset_train, lambda.min.ratio=lmin.ratio)
        beta.rgarrote <- coef(fit2.garrote)
        # (3) Garrote coefficients
        if(split_vars){
          beta_inner[ ,(nl2*(i-1)+1):(nl2*(i-1)+nl1)] <- as.matrix(rbind(beta.rgarrote[1, ], # intercept
                                                               beta.init[-1] * rbind(beta.rgarrote[-1,], beta.rgarrote[-1,])))
        }else{
          beta_inner[ ,(nl2*(i-1)+1):(nl2*(i-1)+nl1)] <- as.matrix(rbind(beta.rgarrote[1, ], # intercept
                                                               beta.init[-1] * rbind(beta.rgarrote[-1,])))
        }
      } # now we have all nl1*nl2 betas               
      
      # validation: compute prediction error
      if(!is.null(data.obj[["clinical"]])){
        clin_offset_test <- apply(c.test %>% dplyr::select(-y), 2, to_numeric) %*% clin_offset_coefs
      }else{clin_offset_test <- rep(0, sum(folds==inner))}
      yhat.test <- apply(beta_inner, 2, function(bcol) cbind(1, var.test) %*% bcol + clin_offset_test)
      if(family == "binomial") yhat.test <- plogis(yhat.test)
      squared_diff <- (y.test - yhat.test)^2
      prederr <- prederr + matrix(colMeans(squared_diff), ncol = nl2, byrow = TRUE) / cv / R
      prederr2 <- prederr2 + matrix((colMeans(squared_diff))^2, ncol = nl2, byrow = TRUE) / cv / R
    } # inner loop end
  } # outer loop end

  # Index of minimal cvm
  index <- which(prederr == min(prederr), arr.ind = TRUE)[1,]
  se.prederror <- sqrt(prederr2 - prederr^2)
  
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
  
  beta <- matrix(0, nrow=ncol(varmat) + 1, ncol=nl1 * nl2) # all coefs for lambdas
  rownames(beta) <- c("(Intercept)", colnames(varmat))
  lambda <- df <- matrix(NA, nrow=nl1, ncol=nl2)

  # (1) Ridge regression
  fit1.init <- glmnet(y = y, x = varmat, family = family, alpha = 0, nlambda = nl1, standardize = standardize_vars, 
                      offset = clin_offset, lambda.min.ratio=lmin.ratio)
  rownames(lambda) <- fit1.init$lambda
  
  # (2) Positive lasso 
  fit2.garrote <- vector(mode = "list", length = nl1)
  for(i in 2:nl1){
    # X garrote = X*beta_ridge
    beta.init <- coef(fit1.init)[, i]
    XB.ridge <-  varmat*rep(beta.init[-1], each=dim(varmat)[1]) 
    if(split_vars){varmat.gar <- XB.ridge[, seq(1, k, 1)] +  XB.ridge[, seq(k + 1, 2*k,1)]
    }else{varmat.gar <- XB.ridge[, seq(1, k, 1)]}
    
    fit2.garrote <- glmnet(y = y, x = varmat.gar, family = family, alpha = 1, lower.limits = 0, 
                           nlambda = nl2, offset = clin_offset, standardize = standardize_vars, lambda.min.ratio=lmin.ratio)
    beta.rgarrote <- coef(fit2.garrote)
    lambda[i, ] <- fit2.garrote$lambda
    df[i, ] <- fit2.garrote$df
    if(split_vars){
      beta[ ,(nl2*(i-1)+1):(nl2*(i-1)+nl1)] <- as.matrix(rbind(beta.rgarrote[1, ], # intercept
                                                           beta.init[-1] * rbind(beta.rgarrote[-1,], beta.rgarrote[-1,])))
    }else{
      beta[ ,(nl2*(i-1)+1):(nl2*(i-1)+nl1)] <- as.matrix(rbind(beta.rgarrote[1, ], # intercept
                                                           beta.init[-1] * rbind(beta.rgarrote[-1,])))
    }
  }               
  
  ## Return
  # Ridge and garrote lambda
  lambda.min <- c(as.numeric(rownames(lambda)[index[1]]), as.numeric(lambda[index[1], index[2]]))
  names(lambda.min)<-c("ridge","garrote")
  df.final <- df[index[1],index[2]]
  
  # Coefficients of best model
  coeffs <- beta[, nl2 * (index[1] - 1) + index[2]]
  dfvars <- length(unique(str_remove_all(names(which(coeffs!=0))[-1], "x.|u.|d.")))

  # Fitted values
  fitted.values <- cbind(1, varmat) %*% coeffs + clin_offset
  
  end <- Sys.time()
  res <- list(call = match.call(), family = family, lambda = lambda, coefficients = coeffs, glmnet.fit.init = fit1.init, 
              glmnet.fit.garrote = fit2.garrote, k = k, kclin = kclin, df = df.final, dfvars = dfvars, cv.pred.err = prederr,
              se.pred.err = se.prederror, extime = as.numeric (end - start, units = "secs"), 
              fit = list(varmat = varmat, varmat.gar = varmat.gar, lambda = lambda, lambda.min = lambda.min, 
                         clin_offset = clin_offset, clin_offset_coefs = clin_offset_coefs,
                         coefficients = coeffs, beta = beta, index.lambda.min = index, 
                         fitted.values = fitted.values, split_vars = split_vars, standardize_vars = standardize_vars))
  attr(res, "class") <-  "ridge-garrote"
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
