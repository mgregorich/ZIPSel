# ============================================================================ #
# Author: Mariella Gregorich
# Date: 02/05/2023
# Info: Functions
# ============================================================================ #


# =============================== GENERAL ======================================

assess_nonzeros <- function(vec, prop_nonzero=0.75){
  # Input: Vector of values
  # Output: Indicate if vector falls below the non-zero threshold (True/False)
  
  out <- (sum(vec > 0)/length(vec)) >= prop_nonzero
  return(out)
}

to_factor <- function(x){
  as.factor(as.character(x))
}

to_numeric <- function(x){
  as.numeric(as.character(x))
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
  r2 <- cor(pred,obs, use="pairwise.complete.obs")^2
  CS <- lm(obs ~ pred)$coefficients[2]
  Cind <- c_index(pred=pred, obs=obs)
  rmse.val <- sqrt(mean((pred-obs)^2))
  mae.val <- mean(abs(pred-obs))

  res <- data.frame(R2 = r2,
                    RMSE = rmse.val,
                    MAE = mae.val,
                    C = Cind,
                    CalbSlope=CS) 
  return(res)
}



simcor.H <- function(k=6, size=c(10,5,8,7,15,50), 
                     rho=rbind(c(.9,.7), c(.7,.7), c(.7,.2), c(.5,.3), c(.9,.85), c(.3,.2)), power=1,
                     epsilon=.08, eidim=2){
  #' Simulating the Hub Matrix (entries filled in using Toeplitz structure)
  #' Implementation by Hardin et al. (DOI: 10.1214/13-AOAS638)
  #' k is the number of groups
  #' size is a vector of length k specifying the size of each group 
  #' rho is a vector of length k specifying base correlation values
  #' epsilon <- (1-min(rho) - 0.75*min(tau) ) - .01
  #' tau_k = (max(rho_k) - min(rho_k) )/ (size_k -2) 
  #' eidim is the space from which the noise is generated, the smaller the more noise
  #' power = 2 makes the correlations stay high
  #'power = 0.5 makes the correlations descent rapidly
  
  ndim <- sum(size)# dim of correlation matrix
  bigcor<- matrix(rep(0, ndim*ndim), ncol=ndim)
  
  ### generating the basic correlation matrix
  for (i in 1:(k) ){
    
    cor <- toeplitz(rho.func(rho[i,1],rho[i,2],power,size[i]) )
    
    if (i==1){bigcor[1:size[1], 1:size[1]] <- cor}
    if (i!=1){bigcor[(sum(size[1:(i-1)]) + 1):sum(size[1:i]),
                     (sum(size[1:(i-1)]) + 1):sum(size[1:i])] <- cor}
  }
  diag(bigcor) <- 1 - epsilon
  
  return(bigcor)
  
  ### adding noise to the correlation matrix
  # eivect <- c( )
  # for (i in 1:ndim) {
  #   ei <- runif(eidim, -1, 1)
  #   eivect <- cbind(eivect, sqrt(epsilon) * ei / sqrt(sum(ei^2) ) )
  # }
  # 
  # 
  # bigE <- t(eivect) %*% eivect
  # cor.nz <- bigcor + bigE
  # cor.nz
}

# ========================== NONNEGATIVE GARROTE ===============================
# Code by Georg Heinze, modified by MG

protogarrote<-function(data.obj, center.interaction.x=0, scale.interaction.x=1, penalties=1, family="gaussian", nlambda=c(11,11), target.L1norm=20, cv=10, outer=2, alpha1=0, alpha2=0.99){
  require(glmnet)
  x<-data.obj[["x"]]
  d<-data.obj[["d"]]
  y<-data.obj[["y"]]
  clinical<-data.obj[["clinical"]]
  interaction.x<-data.obj[["interaction.x"]]
  
  n<-nrow(x)
  k<-ncol(x)
  
  if(k != ncol(d)) stop("d does not match x in dimensions\n")
  if(n != nrow(d) | n!= length(y)) stop("Not equal sample size in variables\n")
  
  if(!is.null(interaction.x)) {
    int.x<-(interaction.x - center.interaction.x)/scale.interaction.x
    clinical <- cbind(clinical, interaction.x)  
    fit.int<-TRUE
  } else {
    fit.int<-FALSE
  }
  
  kclin <- ncol(clinical)
  
  xmat <- as.matrix(cbind(x, d, clinical))  
  if(fit.int) xmat <- as.matrix(cbind(x,  d, x*int.x, d*int.x, clinical))
  if(length(penalties)==1) penalties<-rep(penalties, ncol(xmat))
  # get lambdas
  nl1 <- nlambda[1]
  nl2 <- nlambda[2]
  
  prederr <- rsquare <- prederr2 <- matrix(0, nl1, nl2)
  for(outer.iter in 1:outer){
    folds <- rep(1:cv, each=ceiling(n/10))
    folds <- sample(folds)[1:n]
    for(inner in 1:cv){
      beta<-matrix(0,ncol(xmat)+1,nl1*nl2)
      devel <- (1:n)[folds!=inner]
      test <- (1:n)[folds==inner]
      x.devel<-xmat[devel,]
      x.test <- xmat[test,]
      y.devel <- y[devel]
      y.test <- y[test]
      
      fit1 <- glmnet(y=y.devel, x=x.devel, family=family, alpha=alpha1, nlambda=nl1, penalty.factor =penalties)
      if(length(fit1$lambda)<nl1) stop("In CV loop ", inner, ", number of evaluated lambdas (", length(fit1$lambda), ") smaller than ", nl1, ". Try nlambda[1] = ",length(fit1$lambda),".\n")
      # second step: positive lasso: easier and safer in a loop  
      for(i in 1:nl1){
        #          lambda1<-lambda$lambda1[i]
        beta1 <- coef(fit1)[,i]
        if(!fit.int) {
          partial.eta <- sapply(1:k, function(j) x.devel[,j]*beta1[1+j] + x.devel[,k+j]*beta1[1+k+j])
          partial.clinical <- matrix(sapply(1:ncol(clinical), function(j) x.devel[,2*k+j] * beta1[2*k+j+1]), nrow(x.devel), kclin,byrow=FALSE)
        }
        if(fit.int){
          partial.eta <- sapply(1:k, function(j) x.devel[,j]*beta1[1+j] 
                                + x.devel[,k+j]*beta1[k+1+j] 
                                + x.devel[,2*k+j]*beta1[2*k+1+j] 
                                + x.devel[,3*k+j]*beta1[3*k+1+j])
          partial.clinical <- matrix(sapply(1:ncol(clinical), function(j) x.devel[,4*k+j] * beta1[4*k+j+1]), nrow(x.devel), kclin,byrow=FALSE)
        }
        xmat2 <- as.matrix(cbind(partial.eta, partial.clinical))
        fit2 <- glmnet(y=y.devel, x=xmat2, family=family, alpha=alpha2, lower.limits=0, standardize=FALSE, nlambda=nl2)
        for(ii in 1:nl2){
          if(!fit.int){
            beta[,nl2*(i-1)+ii]<-c(coef(fit2)[1,ii], # intercept
                                   coef(fit2)[2:(k+1),ii]*coef(fit1)[2:(k+1),i], # x
                                   coef(fit2)[2:(k+1),ii]*coef(fit1)[(k+2):(2*k+1),i], # d
                                   coef(fit2)[(k+2):(k+1+kclin),ii]*coef(fit1)[(2*k+2):(2*k+1+kclin),i]) # clinical 
          }
          if(fit.int){
            beta[,nl2*(i-1)+ii]<-c(coef(fit2)[1,ii], # intercept
                                   coef(fit2)[2:(k+1),ii]*coef(fit1)[2:(k+1),i], # x
                                   coef(fit2)[2:(k+1),ii]*coef(fit1)[(k+2):(2*k+1),i], # d
                                   coef(fit2)[2:(k+1),ii]*coef(fit1)[(2*k+2):(3*k+1),i], # x*int.x
                                   coef(fit2)[2:(k+1),ii]*coef(fit1)[(3*k+2):(4*k+1),i], # d*int.x
                                   coef(fit2)[(k+2):(k+1+kclin),ii]*coef(fit1)[(4*k+2):(4*k+1+kclin),i]) # clinical 
          }
          
        }
      } # now we have all nl1*nl2 betas               
      
      # validation: compute prediction error
      
      for(i in 1:nl1){
        for(ii in 1:nl2){
          yhat.test <- cbind(1,x.test) %*% beta[,nl2*(i-1)+ii]
          if(family=="binomial") yhat.test <- plogis(yhat.test)
          prederr[i,ii] <- prederr[i,ii]+mean((y.test-yhat.test)**2)/outer/cv
          prederr2[i,ii] <- prederr2[i,ii]+((mean((y.test-yhat.test)**2))**2)/outer/cv
          if(sd(yhat.test)!=0) rsquare[i,ii] <- cor(y.test, yhat.test)^2
          else rsquare[i,ii] <- 0
        }
      }
    }
  }
  index<-which(prederr==min(prederr), arr.ind=TRUE)
  if(length(index)>2) index<-tail(index,1)
  se.prederr<-sqrt(prederr2 - prederr**2)
  
  # now fit the garrote with all lambda1, lambda2, compute final betas and return object (including which lambda is best)
  
  beta<-matrix(0,ncol(xmat)+1,nl1*nl2)
  lambda <- matrix(0, nl1, nl2)
  # first step: ridge
  fit1 <- glmnet(y=y, x=xmat, family=family, alpha=alpha1, nlambda=nl1, penalty.factor =penalties)
  rownames(lambda) <- fit1$lambda
  
  # second step: positive lasso: easier and safer in a loop  
  fit2 <- vector(mode = "list", length = nl1)
  for(i in 1:nl1){
    #          lambda1<-lambda$lambda1[i]
    beta1 <- coef(fit1)[,i]
    if(!fit.int) {
      partial.eta <- sapply(1:k, function(j) xmat[,j]*beta1[1+j] + xmat[,k+j]*beta1[1+k+j])
      partial.clinical <- matrix(sapply(1:ncol(clinical), function(j) xmat[,2*k+j] * beta1[2*k+j+1]), nrow(xmat), kclin,byrow=FALSE)
    }
    if(fit.int){
      partial.eta <- sapply(1:k, function(j) xmat[,j]*beta1[1+j] 
                            + xmat[,k+j]*beta1[k+1+j] 
                            + xmat[,2*k+j]*beta1[2*k+1+j] 
                            + xmat[,3*k+j]*beta1[3*k+1+j])
      partial.clinical <- matrix(sapply(1:ncol(clinical), function(j) xmat[,4*k+j] * beta1[4*k+j+1]), nrow(xmat), kclin,byrow=FALSE)
    }
    xmat2 <- as.matrix(cbind(partial.eta, partial.clinical))
    fit2[[i]] <- glmnet(y=y, x=xmat2, family=family, alpha=alpha2, lower.limits=0, standardize=FALSE, nlambda=nl2)
    lambda[i,]<-fit2[[i]]$lambda
    if(all(c(i, ii)==index)){
      xmat2.save <- xmat2
    }
    for(ii in 1:nl2){
      if(!fit.int){
        beta[,nl2*(i-1)+ii]<-c(coef(fit2[[i]])[1,ii], # intercept
                               coef(fit2[[i]])[2:(k+1),ii]*coef(fit1)[2:(k+1),i], # x
                               coef(fit2[[i]])[2:(k+1),ii]*coef(fit1)[(k+2):(2*k+1),i], # d
                               coef(fit2[[i]])[(k+2):(k+1+kclin),ii]*coef(fit1)[(2*k+2):(2*k+1+kclin),i]) # clinical 
      }
      if(fit.int){
        beta[,nl2*(i-1)+ii]<-c(coef(fit2[[i]])[1,ii], # intercept
                               coef(fit2[[i]])[2:(k+1),ii]*coef(fit1)[2:(k+1),i], # x
                               coef(fit2[[i]])[2:(k+1),ii]*coef(fit1)[(k+2):(2*k+1),i], # d
                               coef(fit2[[i]])[2:(k+1),ii]*coef(fit1)[(2*k+2):(3*k+1),i], # x*int.x
                               coef(fit2[[i]])[2:(k+1),ii]*coef(fit1)[(3*k+2):(4*k+1),i], # d*int.x
                               coef(fit2[[i]])[(k+2):(k+1+kclin),ii]*coef(fit1)[(4*k+2):(4*k+1+kclin),i]) # clinical 
      }
    }
  } # now we have all nl1*nl2 beta vectors               
  coeffs<-beta[,nl2*(index[1]-1)+index[2]]
  if(!fit.int) {
    names(coeffs)<-rownames(beta)<-c("(Intercept)", colnames(x), colnames(d), colnames(clinical))
    
  }
  if(fit.int){
    names(coeffs)<-rownames(beta)<-c("(Intercept)", colnames(x), colnames(d), paste("i*", colnames(x), sep=""), paste("i*", colnames(d), sep=""), colnames(clinical))
  }
  df <- matrix(unlist(lapply(fit2, function(X) X$df)),nrow=nrow(lambda), ncol=ncol(lambda),byrow=TRUE)
  L1norm <- matrix(unlist(lapply(fit2, function(X) apply(coef(X)[-1,],2,sum))),nrow=nrow(lambda), ncol=ncol(lambda),byrow=TRUE)
  index.lambda.L1norm=c(nl1-1,which(L1norm[nl1-1,]==head(L1norm[nl1-1,][L1norm[nl1-1,]>=target.L1norm],1)))  # second smallest lambda1 hard-coded
  lambda.L1norm=c(as.numeric(rownames(lambda)[index.lambda.L1norm[1]]),lambda[index.lambda.L1norm[1],index.lambda.L1norm[2]])  
  lambda.min=c(as.numeric(rownames(lambda))[index[1]], lambda[index[1],index[2]])
  names(lambda.L1norm)<-names(lambda.min)<-c("ridge","garrote")
  
  res<-list(call=match.call(), fit.int=fit.int, family=family, lambda=lambda, coefficients=beta, glmnet.fit1=fit1, glmnet.fit2=fit2, k=k, kclin=kclin, df=df, L1norm=L1norm, cv.pred.err=prederr,
            cv.rsquare=rsquare, se.pred.err=se.prederr, center.interaction.x=center.interaction.x, scale.interaction.x=scale.interaction.x,
            fit=list(xmat=xmat, xmat2=xmat2, lambda=lambda, lambda.min=lambda.min, 
                     coefficients=coeffs, beta=beta, 
                     coeffs.L1norm=beta[,nl2*(index.lambda.L1norm[1]-1)+index.lambda.L1norm[2]],
                     index.lambda.min=c(index[1], index[2]), 
                     index.lambda.L1norm=index.lambda.L1norm, 
                     lambda.L1norm=lambda.L1norm,
                     fitted.values=cbind(1,xmat) %*% beta[,nl2*(index[1]-1)+index[2]],
                     fitted.values.L1norm=cbind(1,xmat) %*% beta[,nl2*(index.lambda.L1norm[1]-1)+index.lambda.L1norm[2]]))
  attr(res,"class")<-"protogarrotte"
  return(res)
}



predict.protogarrote<-function(obj, newdata, lambda="lambda.min", type="link"){
  if(!obj$fit.int){
    x<-newdata[["x"]]
    d<-newdata[["d"]]
    clinical <- newdata[["clinical"]]
    xmat<-cbind(x,d,clinical)
    if(lambda=="lambda.min"){
      x<-cbind(1, xmat[,names(obj$fit$coefficients[obj$fit$coefficients !=0])[-1]])
      beta <- obj$fit$coefficients[obj$fit$coefficients !=0]
      eta <- x %*% beta
    } else if (lambda=="L1norm"){
      x<-cbind(1, xmat[,names(obj$fit$coeffs.L1norm[obj$fit$coeffs.L1norm !=0])[-1]])
      beta <- obj$fit$coeffs.L1norm[obj$fit$coeffs.L1norm !=0]
      eta <- x %*% beta
    }
    else if (lambda=="all"){
      x <- cbind(1, xmat[,names(obj$fit$coefficients)[-1]])
      eta <- x %*% obj$coefficients
    }
  }
  if(obj$fit.int){
    # search for 'x' and 'd' in newdata:
    x<-newdata[["x"]]
    d<-newdata[["d"]]
    interaction.x <- (newdata[["interaction.x"]] - obj$center.interaction.x)/obj$scale.interaction.x
    
    ix<-x * interaction.x
    id<-d * interaction.x
    
    xmat<-cbind(x, d, ix, id, newdata[["clinical"]], newdata[["interaction.x"]])
    
    if(lambda=="lambda.min"){
      x<-cbind(1, xmat)
      beta <- coefficients.protogarrote(obj)
      eta <- x %*% beta
    } else if (lambda=="L1norm"){
      x<-cbind(1, xmat)
      beta <- obj$fit$coeffs.L1norm
      eta <- x %*% beta
    } else if (lambda=="all"){
      x <- cbind(1, xmat)
      eta <- x %*% coefficients.protogarrote(obj)
    }
  }
  if(type=="response" & obj$family=="binomial") eta <- plogis(eta)
  return(eta)
}


coefficients.protogarrote<-function(obj, lambda="lambda.min"){
  if(lambda=="lambda.min") beta <- obj$fit$coefficients
  if(lambda=="L1norm") beta <- obj$fit$coeffs.L1norm
  if(lambda=="all") beta <- obj$fit$beta
  return(beta)
}  

plot.protogarrote<-function(obj, highlight, jitter=20){
  if(missing(highlight)) highlight<-which(obj$cv.pred.err==min(obj$cv.pred.err), arr.ind=TRUE)
  yrange<-range(obj$cv.pred.err+obj$se.pred.err, obj$cv.pred.err-obj$se.pred.err)
  log.lambda1 <- log(as.numeric(rownames(obj$lambda)))
  dll1<-((1:ncol(obj$lambda))-(ncol(obj$lambda)+1)/2) *diff(log.lambda1)[1]/jitter
  xrange<-range(log(as.numeric(rownames(obj$lambda))))+range(dll1)
  plot(log.lambda1+dll1[1],obj$cv.pred.err[,1], type="l", col="grey", ylab="Prediction error", xlab="log lambda1", xlim=xrange, ylim=yrange)
  for(i in 2:ncol(obj$lambda)) lines(log.lambda1+dll1[i],obj$cv.pred.err[,i], type="l", col="grey")
  for(i in 1:ncol(obj$lambda)) {
    for(j in 1:nrow(obj$lambda)) {
      lines(rep(log.lambda1[j],2)+dll1[i],c(obj$cv.pred.err[j,i]-obj$se.pred.err[j,i], obj$cv.pred.err[j,i]+obj$se.pred.err[j,i]), type="l", col="grey")
    }
  }
  lines(rep(log.lambda1[highlight[1]],2)+dll1[highlight[2]],c(obj$cv.pred.err[highlight]-obj$se.pred.err[highlight], obj$cv.pred.err[highlight]+obj$se.pred.err[highlight]), type="l", col="black", lwd=1.5)
  lines(log.lambda1,obj$cv.pred.err[,highlight[2]], type="l", col="black", lwd=1.5)
}

    
plot.coefficients.protogarrote<-function(obj, order="none", scale=c(1,1,1,1), plot=TRUE, lambda="lambda.min"){
  # obj: a protogarrote object
  # this function shows the proteomics- associated coefficients
  # obj=fit.garrote; order="none"; scale=c(1,1,1,1); plot=TRUE; lambda="lambda.min"
  beta<-coefficients.protogarrote(obj, lambda=lambda)
  beta.hot<-beta[beta!=0]
  beta.u<-beta.hot[substr(names(beta.hot),1,2)=="u."]
  beta.d<-beta.hot[substr(names(beta.hot),1,2)=="d."]
  set.u<-substr(names(beta.u),3,10)
  set.d<-substr(names(beta.d),3,10)
  set.u<-union(set.u, set.d)

  nam.u<-paste("u.", set.u, sep="")
  nam.d<-paste("d.", set.u, sep="")
  beta.u <-beta[nam.u]
  beta.d <- beta[nam.d]
  beta.u[is.na(beta.u)]<-0
  beta.d[is.na(beta.d)]<-0
  if(order == "u") ord<-order(beta.u)
  if(order == "d") ord<-order(beta.d)
  if(order == "none") ord<-1:length(beta.u)
  beta.u <- beta.u[ord]
  beta.d <- beta.d[ord]
  xrange <- range(beta.u*scale[1], beta.d*scale[2])
  if(plot){
    plot(beta.u*scale[1], 1:length(beta.u), type="o", xlim=xrange, xlab="beta (scaled)", ylab="Coefficient")
    points(beta.d*scale[2], 1:length(beta.d), type="o", lty=2, col="red")
    for(i in 1:length(beta.u)) lines(xrange, c(i,i), lty=3, col="grey")
    legend("bottomright", pch=c("o","o"), lty=c(1,2), col=c("black","red"), legend=c("U","D"))
    abline(v=0, col="grey")
  }else{
  return(list(beta=cbind(beta.u, beta.d), selected=set.u))}
  
}  

  
# ============================== METHODS =======================================

perform_ridge <- function(data.obj, penalties=1, family="gaussian", nlambda=c(10,10), target.L1norm=20, cv=10, R=2, alpha1=0, alpha2=1){
  u<-data.obj[["u"]]
  d<-data.obj[["d"]]
  y<-data.obj[["y"]]
  clinical<-data.frame("y"=y, data.obj[["clinical"]])
  udmat <- as.matrix(cbind(u, d))  
  
  n<-nrow(u)
  k<-ncol(u)
  kclin <- ncol(clinical)
  
  pfvector <- rep(1, 2*k)
  if(penalty=="component"){pfvector <- rep(pf, each=k)}
  
  prederr <- rsquare <- prederr2 <- matrix(0, nl1, nl2)
  for(outer in 1:R){
    folds <- sample(rep(1:cv, ceiling(n/cv)))[1:n]
    
    # (0) Clinical offset
    fit.clin <- glm(y~., data=clinical, family="gaussian")
    clin_offset_coefs <- fit.clin$coefficients[-1]
    clin_offset <- apply(c.train %>% select(-y),2,to_numeric)%*% clin_offset_coefs
    
    # (1) CV Ridge model with offset
    cv_ridge <- cv.glmnet(x=udmat %>% select(-c(log2eGFR, age, sex, preds)) %>% as.matrix(), y=udmat$log2eGFR, 
                          alpha = 0, standardize=FALSE, offset = clin_offset, penalty.factor=pfvector)
    fit.ridge <- cv_ridge$glmnet.fit
  }
}


perform_rgarrote<-function(data.obj, penalties=1, family="gaussian", nlambda=c(10,10), cv=10, R=2, alpha1=0, alpha2=1){
  # data.obj=data.obj; penalties=1; family="gaussian"; nlambda=c(10,10);
  # target.L1norm=20; cv=10; R=2; alpha1=0; alpha2=1

  u<-data.obj[["u"]]
  d<-data.obj[["d"]]
  y<-data.obj[["y"]]
  clinical<-data.frame("y"=y, data.obj[["clinical"]])
  udmat <- as.matrix(cbind(u, d))  
  
  n<-nrow(u)
  k<-ncol(u)
  kclin <- ncol(clinical)
  
  # get number of lambdas
  nl1 <- nlambda[1]
  nl2 <- nlambda[2]
  
  if(k != ncol(d)) stop("d does not match x in dimensions\n")
  if(n != nrow(d) | n!= length(y)) stop("Not equal sample size in variables\n")
  if(length(penalties)==1) penalties<-rep(penalties, 2*k)
  
  prederr <- rsquare <- prederr2 <- matrix(0, nl1, nl2)
  for(outer in 1:R){
    folds <- sample(rep(1:cv, ceiling(n/cv)))[1:n]
    for(inner in 1:cv){
      beta<-matrix(0,ncol(udmat)+1,nl1*nl2)
      ud.train <- udmat[(1:n)[folds!=inner],]
      ud.test <- udmat[(1:n)[folds==inner],]
      c.train <- clinical[(1:n)[folds!=inner],]
      c.test <- clinical[(1:n)[folds==inner],]
      y.train <- y[(1:n)[folds!=inner]]
      y.test <- y[(1:n)[folds==inner]]
      
      # (0) Clinical offset
      fit.clin <- glm(y~., data=c.train, family="gaussian")
      clin_offset_coefs <- fit.clin$coefficients[-1]
      clin_offset_train <- apply(c.train %>% select(-y),2,to_numeric)%*% clin_offset_coefs
      
      # (1) Ridge regression
      fit1.ridge <- glmnet(y=y.train, x=ud.train, family=family, standardize = FALSE,
                           alpha=alpha1, nlambda=nl1, penalty.factor =penalties, offset=clin_offset_train)
      
      # (2) Lasso regression with restriction of positive coefficients for non-negative shrinkage factors
      for(i in 1:nl1){
        # X garrote = X*beta_ridge
        B.ridge <- diag(coef(fit1.ridge)[-1,i])
        XB.ridge <- ud.train %*% B.ridge
        udmat.gar <- XB.ridge[,seq(1,k,1)]+ XB.ridge[,seq(k+1, 2*k,1)]
        
        # Positive Lasso: lower.limits ensures positiveness of coeffs
        fit2.garrote <- glmnet(y=y.train, x=udmat.gar, family=family, 
                               alpha=alpha2, lower.limits=0, standardize=FALSE, nlambda=nl2, offset = clin_offset_train)
        
        # (3) Garrote coefficients
        for(ii in 1:nl2){
          beta[,nl2*(i-1)+ii]<-c(coef(fit2.garrote)[1,ii], # intercept
                                 coef(fit2.garrote)[2:(k+1),ii]*coef(fit1.ridge)[2:(k+1),i], # x
                                 coef(fit2.garrote)[2:(k+1),ii]*coef(fit1.ridge)[(k+2):(2*k+1),i]) # d
        }
      } # now we have all nl1*nl2 betas               
      
      # validation: compute prediction error
      clin_offset_test <- apply(c.test %>% select(-y),2,to_numeric)%*% clin_offset_coefs
      for(i in 1:nl1){
        for(ii in 1:nl2){
          yhat.test <- cbind(1,ud.test) %*% beta[,nl2*(i-1)+ii] +clin_offset_test
          if(family=="binomial") yhat.test <- plogis(yhat.test)
          prederr[i,ii] <- prederr[i,ii]+mean((y.test-yhat.test)**2)/cv/R
          prederr2[i,ii] <- prederr2[i,ii]+((mean((y.test-yhat.test)**2))**2)/cv/R
          if(sd(yhat.test)!=0) rsquare[i,ii] <- cor(y.test, yhat.test)^2
          else rsquare[i,ii] <- 0
        }
      }
    } # inner loop end
  } # outer loop end
  
  
  index<-which(prederr==min(prederr), arr.ind=TRUE) 
  if(length(index)>2) index<-tail(index,1)
  se.prederr<-sqrt(prederr2 - prederr^2)

  ## Final model
  # (0) Clinical offset
  fit.clin <- glm(y~., data=clinical, family="gaussian")
  clin_offset_coefs <- fit.clin$coefficients[-1]
  clin_offset <- apply(clinical %>% select(-y),2,to_numeric)%*% clin_offset_coefs
  
  # (1) Ridge regression
  lambda <- matrix(0, nl1, nl2)
  fit1.ridge <- glmnet(y=y, x=udmat, family=family, alpha=alpha1, nlambda=nl1, penalty.factor =penalties, standardize = FALSE, offset=clin_offset)
  rownames(lambda) <- fit1.ridge$lambda
  
  # (2) Positive lasso 
  beta<-matrix(0,ncol(udmat)+1,nl1*nl2)
  fit.lasso <- vector(mode = "list", length = nl1)
  for(i in 1:nl1){
    # X garrote = X*beta_ridge
    B.ridge <- diag(coef(fit1.ridge)[-1,i])
    XB.ridge <- udmat %*% B.ridge
    udmat.gar <- XB.ridge[,seq(1,k,1)]+ XB.ridge[,seq(k+1, 2*k,1)]
    
    fit.lasso[[i]] <- glmnet(y=y, x=udmat.gar, family=family, alpha=alpha2, 
                             lower.limits=0, standardize=FALSE, nlambda=nl2, offset = clin_offset)
    lambda[i,] <- fit.lasso[[i]]$lambda
    if(all(c(i, ii)==index)){
      udmat.gar.save <- udmat.gar
    }
    for(ii in 1:nl2){
      beta[,nl2*(i-1)+ii]<-c(coef(fit.lasso[[i]])[1,ii], # intercept
                             coef(fit.lasso[[i]])[2:(k+1),ii]*coef(fit1.ridge)[2:(k+1),i], # x
                             coef(fit.lasso[[i]])[2:(k+1),ii]*coef(fit1.ridge)[(k+2):(2*k+1),i]) # d 
      
      
    }
  } # now we have all nl1*nl2 beta vectors               
  coeffs<-beta[,nl2*(index[1]-1)+index[2]]
  names(coeffs)<-rownames(beta)<-c("(Intercept)", colnames(u), colnames(d))
  
  # Return
  df <- matrix(unlist(lapply(fit.lasso, function(X) X$df)),nrow=nrow(lambda), ncol=ncol(lambda),byrow=TRUE)
  lambda.min=c(as.numeric(rownames(lambda))[index[1]], lambda[index[1],index[2]])
  fitted.values <- cbind(1,udmat) %*% beta[,nl2*(index[1]-1)+index[2]] + clin_offset

  res<-list(call=match.call(), family=family, lambda=lambda, coefficients=beta, glmnet.fit.ridge=fit1.ridge, 
            glmnet.fit.lasso=fit.lasso, k=k, kclin=kclin, df=df, cv.pred.err=prederr,
            cv.rsquare=rsquare, se.pred.err=se.prederr, 
            fit=list(udmat=udmat, udmat.gar=udmat.gar, lambda=lambda, lambda.min=lambda.min, 
                     clin_offset=clin_offset, clin_offset_coefs=clin_offset_coefs,
                     coefficients=coeffs, beta=beta, 
                     index.lambda.min=c(index[1], index[2]), 
                     fitted.values=fitted.values))
  attr(res,"class")<-"ridge-garrote"
  return(res)
}


predict_rgarrote<-function(obj, newdata, lambda="lambda.min", type="link"){
  u<-newdata[["u"]]
  d<-newdata[["d"]]
  clinical <- newdata[["clinical"]]
  udmat<-cbind(u,d)
  
  new_offset <- apply(clinical,2,to_numeric)%*% obj$fit$clin_offset_coefs
  if(lambda=="lambda.min"){
    x<-cbind(1, udmat[,names(obj$fit$coefficients[obj$fit$coefficients !=0])[-1]])
    beta <- obj$fit$coefficients[obj$fit$coefficients !=0]
    linpred <- x %*% beta + new_offset
  } else if (lambda=="L1norm"){
    x<-cbind(1, udmat[,names(obj$fit$coeffs.L1norm[obj$fit$coeffs.L1norm !=0])[-1]])
    beta <- obj$fit$coeffs.L1norm[obj$fit$coeffs.L1norm !=0]
    linpred <- x %*% beta + new_offset
  }
  else if (lambda=="all"){
    x <- cbind(1, udmat[,names(obj$fit$coefficients)[-1]])
    linpred <- x %*% obj$coefficients + new_offset
  }
  if(type=="response" & obj$family=="binomial") linpred <- plogis(linpred)
  return(linpred)
}

perform_rlasso <- function(data.obj, family="gaussian", nlambda=c(10,10), cv=10, R=2, alpha1=0, alpha2=1){
  
 # data.obj <- list("y" = log2(data_mos$eGFR), "clinical"=c, "x"=ximp, "d"=d); nlambda=c(11,11); cv=10; R=2; alpha1=0; alpha2=1; family="gaussian"
  x<-data.obj[["x"]]
  u<-data.obj[["u"]]
  d<-data.obj[["d"]]
  y<-data.obj[["y"]]
  clinical<-data.frame("y"=y, data.obj[["clinical"]])
  xmat <- as.matrix(x)  
  udmat <- as.matrix(cbind(u, d))  
  
  n<-nrow(x)
  k<-ncol(x)
  kclin <- ncol(clinical)
  
  # get number of lambdas
  nl1 <- nlambda[1]
  nl2 <- nlambda[2]
  
  if(k != ncol(d)) stop("d does not match x in dimensions\n")
  if(n != nrow(d) | n!= length(y)) stop("Not equal sample size in variables\n")

  prederr <- rsquare <- prederr2 <- matrix(0, nl1, nl2)
  for(outer in 1:R){
    folds <- sample(rep(1:cv, ceiling(n/cv)))[1:n]
    for(inner in 1:cv){
      beta<-matrix(0,ncol(xmat)+1,nl1*nl2)
      x.train <- xmat[(1:n)[folds!=inner],]
      x.test <- xmat[(1:n)[folds==inner],]
      ud.train <- udmat[(1:n)[folds!=inner],]
      ud.test <- udmat[(1:n)[folds==inner],]
      y.train <- y[(1:n)[folds!=inner]]
      y.test <- y[(1:n)[folds==inner]]
      c.train <- clinical[(1:n)[folds!=inner],]
      c.test <- clinical[(1:n)[folds==inner],]
      
      # (0) Clinical offset
      fit.clin <- glm(y~., data=c.train, family="gaussian")
      clin_offset_coefs <- fit.clin$coefficients[-1]
      clin_offset_train <- apply(c.train %>% select(-y),2,to_numeric)%*% clin_offset_coefs
      
      # (1) Ridge regression
      fit1.ridge <- glmnet(y=y.train, x=ud.train, family=family, alpha=alpha1, nlambda=nl1, offset = clin_offset_train)
      
      # (2) Lasso regression 
      for(i in 1:nl1){
        # Penalty
        b.ridge <- c(coef(fit1.ridge)[-1,i])
        penalty_peptides <- abs(b.ridge[1:k])+ abs(b.ridge[(k+1):(2*k)])
        penalty <- 1 / c(penalty_peptides)^1
        
        # Lasso: lower.limits ensures positiveness of coeffs
        fit2.lasso <- glmnet(y=y.train, x=x.train, alpha = alpha2, nlambda = nl2, penalty.factor = penalty, offset = clin_offset_train)
        b.alasso <- coef(fit2.lasso)
      
        beta[1:(ncol(xmat)+1),(nl2*(i-1)+1):(nl2*(i-1)+nl1)]<-as.matrix(b.alasso)
        
      } # now we have all nl1*nl2 betas               
      
      # validation: compute prediction error
      clin_offset_test <- apply(c.test %>% select(-y),2,to_numeric)%*% clin_offset_coefs
      for(i in 1:nl1){
        for(ii in 1:nl2){
          yhat.test <- cbind(1,x.test) %*% beta[,nl2*(i-1)+ii] + clin_offset_test
          if(family=="binomial") yhat.test <- plogis(yhat.test)
          prederr[i,ii] <- prederr[i,ii]+mean((y.test-yhat.test)**2)/cv/R
          prederr2[i,ii] <- prederr2[i,ii]+((mean((y.test-yhat.test)**2))**2)/cv/R
          if(sd(yhat.test)!=0) rsquare[i,ii] <- cor(y.test, yhat.test)^2
          else rsquare[i,ii] <- 0
        }
      }
    } # inner loop end
  } # outer loop end
  
  
  index<-which(prederr==min(prederr), arr.ind=TRUE) 
  if(length(index)>2) index<-tail(index,1)
  se.prederr<-sqrt(prederr2 - prederr^2)
  
  ## Final model
  # (0) Clinical offset
  fit.clin <- glm(y~., data=clinical, family="gaussian")
  clin_offset_coefs <- fit.clin$coefficients[-1]
  clin_offset <- apply(clinical %>% select(-y),2,to_numeric)%*% clin_offset_coefs
  
  # (1) Ridge regression
  lambda <- matrix(0, nl1, nl2)
  fit1.ridge <- glmnet(y=y, x=udmat, family=family, alpha=alpha1, nlambda=nl1)
  rownames(lambda) <- fit1.ridge$lambda
  
  # (2) Positive lasso 
  beta<-matrix(0,ncol(xmat)+1,nl1*nl2)
  fit2.alasso <- vector(mode = "list", length = nl1)
  for(i in 1:nl1){
    # Penalty
    b.ridge <- c(coef(fit1.ridge)[-1,i])
    penalty_peptides <- abs(b.ridge[1:k])+ abs(b.ridge[(k+1):(2*k)])
    penalty <- 1 / c(penalty_peptides)^1
    
    fit2.alasso[[i]] <- glmnet(y=y, x=xmat, family=family, alpha=alpha2, nlambda=nl2, offset = clin_offset)
    lambda[i,]<-fit2.alasso[[i]]$lambda
    b.alasso <- coef(fit2.alasso[[i]])
    
    beta[1:(ncol(xmat)+1),(nl2*(i-1)+1):(nl2*(i-1)+nl1)] <- as.matrix(b.alasso)
  } # now we have all nl1*nl2 beta vectors               
  coeffs<-beta[,nl2*(index[1]-1)+index[2]]
  names(coeffs)<-rownames(beta)<-c("(Intercept)", colnames(x))
  
  # Return
  df <- matrix(unlist(lapply(fit2.alasso, function(X) X$df)),nrow=nrow(lambda), ncol=ncol(lambda),byrow=TRUE)
  lambda.min=c(as.numeric(rownames(lambda))[index[1]], lambda[index[1],index[2]])
  names(lambda.min)<-c("ridge","alasso")
  fitted.values <- cbind(1,xmat) %*% beta[,nl2*(index[1]-1)+index[2]] + clin_offset
  
  res<-list(call=match.call(), family=family, lambda=lambda, coefficients=beta, glmnet.fit1.ridge=fit1.ridge, glmnet.fit.alasso=fit2.alasso, 
            k=k, kclin=kclin, df=df,  cv.pred.err=prederr,
            cv.rsquare=rsquare, se.pred.err=se.prederr, 
            fit=list(xmat=xmat, lambda=lambda, lambda.min=lambda.min, 
                     clin_offset=clin_offset, clin_offset_coefs=clin_offset_coefs,
                     coefficients=coeffs, beta=beta, 
                     index.lambda.min=c(index[1], index[2]), 
                     fitted.values=fitted.values))
  attr(res,"class")<-"ridge-lasso"
  return(res)
}

predict_rlasso<-function(obj, newdata, type="link"){
  x<-newdata[["x"]]
  clinical <- newdata[["clinical"]]
  xmat <- as.matrix(x)  

  new_offset <- apply(clinical,2,to_numeric)%*% obj$fit$clin_offset_coefs
  x<-cbind(1, xmat[,names(obj$fit$coefficients[obj$fit$coefficients !=0])[-1]])
  beta <- obj$fit$coefficients[obj$fit$coefficients !=0]
  linpred <- x %*% beta + new_offset

  if(type=="response" & obj$family=="binomial") linpred <- plogis(linpred)
  return(linpred)
}

perform_lridge <- function(data.obj, family="gaussian", nlambda=c(10,10), cv=10, R=2, alpha1=1, alpha2=0, 
                           penalty_type="combined", pflist=NULL){
  
  # data.obj=data.obj; family="gaussian"; nlambda=c(11,11); cv=10; R=2; alpha1=1; alpha2=0; pflist=NULL; penalty_type="combined"
  x<-data.obj[["x"]]
  u<-data.obj[["u"]]
  d<-data.obj[["d"]]
  y<-data.obj[["y"]]
  clinical<-data.frame("y"=y, data.obj[["clinical"]])
  xmat <- as.matrix(x)  
  udmat <- as.matrix(cbind(u, d))  
  
  n<-nrow(x)
  k<-ncol(x)
  kclin <- ncol(clinical)
  
  # get number of lambdas
  nl1 <- nlambda[1]
  nl2 <- nlambda[2]
  if(is.null(pflist) & penalty_type=="independent"){stop("pflist must be specified if penalty_type='independent'!")}
  ptype <- ifelse(penalty_type=="combined", list(c(1,1)), pflist)
  
  if(k != ncol(d)) stop("d does not match x in dimensions\n")
  if(n != nrow(d) | n!= length(y)) stop("Not equal sample size in variables\n")
  
  prederr <- rsquare <- prederr2 <- matrix(0, nl1, nl2)
  for(outer in 1:R){
    folds <- sample(rep(1:cv, ceiling(n/cv)))[1:n]
    for(inner in 1:cv){
      beta<-matrix(0,ncol(udmat)+1,nl1*nl2)
      rownames(beta) <- c("(Intercept)", colnames(udmat))
      
      x.train <- xmat[(1:n)[folds!=inner],]
      x.test <- xmat[(1:n)[folds==inner],]
      ud.train <- udmat[(1:n)[folds!=inner],]
      ud.test <- udmat[(1:n)[folds==inner],]
      c.train <- clinical[(1:n)[folds!=inner],]
      c.test <- clinical[(1:n)[folds==inner],]
      y.train <- y[(1:n)[folds!=inner]]
      y.test <- y[(1:n)[folds==inner]]
      
      # (0) Clinical offset
      fit.clin <- glm(y~., data=c.train, family="gaussian")
      clin_offset_coefs <- fit.clin$coefficients[-1]
      clin_offset_train <- apply(c.train %>% select(-y),2,to_numeric)%*% clin_offset_coefs
      
      # (1) Lasso regression
      fit1.lasso <- glmnet(y=y.train, x=ud.train, family=family, alpha=alpha1, nlambda=nl1, offset = clin_offset_train)
      
      # (2) Ridge regression 
      for(i in 1:nl1){
        # Include x and d part of selected vars
        b.lasso <- coef(fit1.lasso)[,i]
        nonzero.coefs <- names(b.lasso)[which(b.lasso!=0)]
        nonzero.peps <- str_remove_all(nonzero.coefs, "x.")[-1]
        if(length(nonzero.peps)>0){
          ud.train2 <- ud.train[,str_detect(colnames(ud.train), pattern=paste0(nonzero.peps, collapse="|"))]
          fit2.ridge <- glmnet(y=y.train, x=ud.train2, alpha = alpha2, nlambda = nl2, offset = clin_offset_train)
          b.ridge <- coef(fit2.ridge)
          beta[rownames(b.ridge),(nl2*(i-1)+1):(nl2*(i-1)+nl1)]<-as.matrix(b.ridge)          
        }
      } # now we have all nl1*nl2 betas               
      
      # validation: compute prediction error
      clin_offset_test <- apply(c.test %>% select(-y),2,to_numeric)%*% clin_offset_coefs
      for(i in 1:nl1){
        for(ii in 1:nl2){
          yhat.test <- cbind(1,ud.test) %*% beta[,nl2*(i-1)+ii] + clin_offset_test
          if(family=="binomial") yhat.test <- plogis(yhat.test)
          prederr[i,ii] <- prederr[i,ii]+mean((y.test-yhat.test)**2)/cv/R
          prederr2[i,ii] <- prederr2[i,ii]+((mean((y.test-yhat.test)**2))**2)/cv/R
          if(sd(yhat.test)!=0) rsquare[i,ii] <- cor(y.test, yhat.test)^2
          else rsquare[i,ii] <- 0
        }
      }
    } # inner loop end
  } # outer loop end
  
  
  index<-which(prederr==min(prederr), arr.ind=TRUE) 
  if(length(index)>2) index<-tail(index,1)
  se.prederr<-sqrt(prederr2 - prederr^2)
  
  # Final betas for all lambda combinations and return object (including which lambda is best)
  beta<-matrix(0,ncol(udmat)+1,nl1*nl2)
  rownames(beta) <- c("(Intercept)", colnames(udmat))
  
  # (0) Clinical offset
  fit.clin <- glm(y~., data=clinical, family="gaussian")
  clin_offset_coefs <- fit.clin$coefficients[-1]
  clin_offset <- apply(clinical %>% select(-y),2,to_numeric)%*% clin_offset_coefs
  
  lambda <- matrix(0, nl1, nl2)
  # (1) Lasso regression
  fit1.lasso <- glmnet(y=y, x=xmat, family=family, alpha=alpha1, nlambda=nl1, offset = clin_offset)
  rownames(lambda) <- fit1.lasso$lambda
  
  # (2) Ridge regression with selected variables (continuous and binary part)
  fit2.ridge <- vector(mode = "list", length = nl1)
  for(i in 1:nl1){
    # Include x and d part of selected vars
    b.lasso <- coef(fit1.lasso)[,i]
    nonzero.coefs <- names(b.lasso)[which(b.lasso!=0)]
    nonzero.peps <- str_remove_all(nonzero.coefs, "x.")[-1]
    
    if(length(nonzero.peps)>0){
      udmat2 <- udmat[,str_detect(colnames(udmat), pattern=paste0(nonzero.peps, collapse="|"))]
      fit2.ridge[[i]] <- glmnet(y=y, x=udmat2, alpha = alpha2, nlambda = nl2, offset = clin_offset)
      b.ridge <- coef(fit2.ridge[[i]])
      lambda[i,]<-fit2.ridge[[i]]$lambda
      beta[rownames(b.ridge),(nl2*(i-1)+1):(nl2*(i-1)+nl1)]<-as.matrix(b.ridge)          
    }else{
      fit2.ridge[[i]]$df <- rep(NA, nl2)
    }
  } # now we have all nl1*nl2 beta vectors               
  coeffs<-beta[,nl2*(index[1]-1)+index[2]]
  names(coeffs)<-rownames(beta)<-c("(Intercept)", colnames(u), colnames(d))
  
  # Return
  df <- matrix(unlist(lapply(fit2.ridge, function(X) X$df)),nrow=nrow(lambda), ncol=ncol(lambda),byrow=TRUE)
  lambda.min=c(as.numeric(rownames(lambda))[index[1]], lambda[index[1],index[2]])
  names(lambda.min)<-c("lasso","ridge")
  fitted.values <- cbind(1,udmat) %*% beta[,nl2*(index[1]-1)+index[2]] + clin_offset
                 
  res<-list(call=match.call(), family=family, lambda=lambda, coefficients=beta, glmnet.fit1.lasso=fit1.lasso, glmnet.fit2.ridge=fit2.ridge, 
            k=k, kclin=kclin, df=df,  cv.pred.err=prederr,
            cv.rsquare=rsquare, se.pred.err=se.prederr, 
            fit=list(xmat=xmat, udmat=udmat, lambda=lambda, lambda.min=lambda.min, 
                     clin_offset=clin_offset, clin_offset_coefs=clin_offset_coefs,
                     coefficients=coeffs, beta=beta, 
                     index.lambda.min=c(index[1], index[2]), 
                     fitted.values=fitted.values)) 
  attr(res,"class")<-"lasso-ridge"
  return(res)
}

predict_lridge<-function(obj, newdata, type="link"){
  x<-newdata[["x"]]
  u<-newdata[["u"]]
  d<-newdata[["d"]]
  clinical <- newdata[["clinical"]]
  xmat<-cbind(u,d)
  
  new_offset <- apply(clinical,2,to_numeric)%*% obj$fit$clin_offset_coefs
  x <-cbind(1, xmat[,names(obj$fit$coefficients[obj$fit$coefficients !=0])[-1]])
  beta <- obj$fit$coefficients[obj$fit$coefficients !=0]
  linpred <- x %*% beta + new_offset

  if(type=="response" & obj$family=="binomial") linpred <- plogis(linpred)
  return(linpred)
}


  
  
  
