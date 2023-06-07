# ============================================================================ #
# Author: Mariella Gregorich
# Date: 02/05/2023
# Info: Functions
# ============================================================================ #


# =============================== GENERAL ======================================

assess_nonzeros <- function(vec, thresh=0.75){
  # Input: Vector of values
  # Output: Indicate if vector falls below the non-zero threshold (True/False)
  
  out <- (sum(vec > 0)/length(vec)) >= thresh
  return(out)
}

to_factor <- function(x){
  as.character(as.factor(x))
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

# ========================== NETWORK-BASED GROUP PENALIZATION =========================================

# Obtention of the adjacency matrix using WGCNA
# @param data matrix of features  
# Descr: power adjacency matrix; power is chosen st network approximate scale-free topology
GetAdjacencyWGCNA=function(data){
  # data = xt
  powers = c(c(1:10), seq(from = 12, to=20, by=2))
  R_sq<-pickSoftThreshold(data, powerVector = powers, corFnc = cor, verbose = 1, RsquaredCut=0.8, moreNetworkConcepts = T)
  power <- min(R_sq$fitIndices[which.max(R_sq$fitIndices[,2]),1],R_sq$powerEstimate,min(which(R_sq$fitIndices[,2]>0.5)), na.rm=T) 
  adjacency = adjacency(data, power = power,type='unsigned')
  
  results<-list(adjacency,power,R_sq)
  return(results)
}

# Clustering analysis based on the dynamic tree cut algorithm
# @param adjacency adjacency matrix of the network studied  
# @param nodes vector of names of the features in the network
GetClusters=function(adjacency,nodes, data1, method) {
  
  
  # Turn adjacency into topological overlap
  TOM = WGCNA::TOMsimilarity(adjacency, verbose=0);
  dissTOM = 1-TOM
  # Call the hierarchical clustering function
  geneTree = hclust(as.dist(dissTOM), method = "average");
  
  # We like large modules, so we set the minimum module size relatively high:
  minModuleSize = 10;
  # :
  dynamicMods = cutreeDynamic(dendro = geneTree, distM = dissTOM,
                              deepSplit = 2, pamRespectsDendro = FALSE,
                              minClusterSize = minModuleSize, verbose=0);
  clusters = as.numeric(dynamicMods)
  nodes=data.frame(cbind(nodes,clusters)) # grey - no cluster
  return(nodes)
}


NetPredCV <- function(fold, data1, pheno1, data.val, pheno.val, nfoldsIn, seed){
  # fold=2; method1=method[1]; method2=method[2]; data1=dataCKD[-folds[[fold]],]; pheno1=phenotype[-folds[[fold]]]; data.val=dataCKD[folds[[fold]],]; pheno.val=phenotype[folds[[fold]]]; nfoldsIn=3; seed=100
  
  nodes=c(1:ncol(data1[,network_vars]))
  coef.cv=predict.fit=modules=ncluster=NULL
  
  # 1. WGCNA
  adjacency=GetAdjacencyWGCNA(data1[,network_vars])

  # 2. CLUSTERING
  clusters=GetClusters(adjacency[[1]],nodes, data1, method[1])
  modules<-cbind(modules,clusters[,2])
  ncluster<-c(ncluster,length(unique(clusters[,2])))
  clusters0=as.numeric(clusters[,2])
  clusters0[which(clusters0==0)]=max(clusters0)+1
  result_list <- list()
  clusters0 <- c(rep(max(clusters0)+1, length(clinical_vars)), clusters0)
  if(!is.null(add_vars)){clusters0 <- c(rep(max(clusters0)+1, length(add_vars)), clusters0)}
  res_clust <- data.frame(var=colnames(dataCKD), cluster=clusters0)
  
  # 3. MODEL
  fit.gglasso=gglasso(as.matrix(data1),pheno1,nlam=150,lambda.factor=0.01,loss="logit",group=clusters0)
  cv.fit.gglasso=cv.gglasso(as.matrix(data1),pheno1,lambda.factor=0.01,nlam=150,loss="logit",group=clusters0, nfolds=nfoldsIn)
  predict_val <- unlist(predict(fit.gglasso,as.matrix(data.val),s=cv.fit.gglasso$lambda.min,type="link"))
  predict.fit<-(exp(predict_val)/(1+exp(predict_val)))
  coef.cv<- exp(coef(fit.gglasso,as.matrix(data.val),s=cv.fit.gglasso$lambda.min)[-1])
  cv.val <- rms::val.prob(predict.fit, pheno.val, pl=F)
  
  result_list <- list("fold"=fold, "predict.fit"=predict.fit, "coef.cv"=coef.cv, "model.val"=cv.val,
                      "adjmatrix"=adjacency[[1]], "adjpower"=adjacency[[2]], "R_sq"=adjacency[[3]],
                      "clusters"=res_clust, "fit.model"=fit.gglasso, "cv.fit.penal"=cv.fit.gglasso$lambda.min )
    
  
  
  return(result_list)
}

# ========================== NONNEGATIVE GARROTE ===============================
# Code by Georg Heinze, modified by MG

perform_garrote<-function(data.obj, penalties=1, family="gaussian", nlambda=c(11,11), 
                       target.L1norm=20, cv=10, R=2, alpha1=0, alpha2=1){
  # data.obj=data.obj; center.interaction.x=0; scale.interaction.x=1; penalties=1; family="gaussian"; nlambda=c(11,11);
  # target.L1norm=20; cv=c(10,10); outer=2; alpha1=0; alpha2=0.99

  x<-data.obj[["x"]]
  d<-data.obj[["d"]]
  y<-data.obj[["y"]]
  clinical<-data.obj[["clinical"]]
  xmat <- as.matrix(cbind(x, d, clinical))  
  
  n<-nrow(x)
  k<-ncol(x)
  kclin <- ncol(clinical)
  
  # get number of lambdas
  nl1 <- nlambda[1]
  nl2 <- nlambda[2]
  
  if(k != ncol(d)) stop("d does not match x in dimensions\n")
  if(n != nrow(d) | n!= length(y)) stop("Not equal sample size in variables\n")
  if(length(penalties)==1) penalties<-rep(penalties, 2*k+kclin)
  
  prederr <- rsquare <- prederr2 <- matrix(0, nl1, nl2)
  for(outer in 1:R){
    folds <- sample(rep(1:cv, ceiling(n/cv)))[1:n]
    for(inner in 1:cv){
      beta<-matrix(0,ncol(xmat)+1,nl1*nl2)
      x.train <- xmat[(1:n)[folds!=inner],]
      x.test <- xmat[(1:n)[folds==inner],]
      y.train <- y[(1:n)[folds!=inner]]
      y.test <- y[(1:n)[folds==inner]]
      
      # (1) Ridge regression
      fit1.ridge <- glmnet(y=y.train, x=x.train, family=family, alpha=alpha1, nlambda=nl1, penalty.factor =penalties)
      
      # (2) Lasso regression with restriction of positive coefficients for non-negative shrinkage factors
      for(i in 1:nl1){
        # X garrote = X*beta_ridge
        B.ridge <- diag(coef(fit1.ridge)[-1,i])
        XB.ridge <- x.train %*% B.ridge
        xmat.gar <- cbind(XB.ridge[,seq(1,k,1)]+ XB.ridge[,seq(k+1, 2*k,1)], XB.ridge[,seq(2*k+1, ncol(XB.ridge),1)])
        
        # Positive Lasso: lower.limits ensures positiveness of coeffs
        fit2.garrote <- glmnet(y=y.train, x=xmat.gar, family=family, alpha=alpha2, lower.limits=0, standardize=FALSE, nlambda=nl2)

        # (3) Garrote coefficients
        for(ii in 1:nl2){
          beta[,nl2*(i-1)+ii]<-c(coef(fit2.garrote)[1,ii], # intercept
                                 coef(fit2.garrote)[2:(k+1),ii]*coef(fit1.ridge)[2:(k+1),i], # x
                                 coef(fit2.garrote)[2:(k+1),ii]*coef(fit1.ridge)[(k+2):(2*k+1),i], # d
                                 coef(fit2.garrote)[(k+2):(k+1+kclin),ii]*coef(fit1.ridge)[(2*k+2):(2*k+1+kclin),i]) # clinical 
        }
      } # now we have all nl1*nl2 betas               
      
      # validation: compute prediction error
      for(i in 1:nl1){
        for(ii in 1:nl2){
          yhat.test <- cbind(1,x.test) %*% beta[,nl2*(i-1)+ii]
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
  
  # Final betas for all lmabda combinations and return object (including which lambda is best)
  beta<-matrix(0,ncol(xmat)+1,nl1*nl2)
  lambda <- matrix(0, nl1, nl2)
  # (1) Ridge regression
  fit1.ridge <- glmnet(y=y, x=xmat, family=family, alpha=alpha1, nlambda=nl1, penalty.factor =penalties)
  rownames(lambda) <- fit1.ridge$lambda
  
  # (2) Positive lasso 
  fit.lasso <- vector(mode = "list", length = nl1)
  for(i in 1:nl1){
    # X garrote = X*beta_ridge
    B.ridge <- diag(coef(fit1.ridge)[-1,i])
    XB.ridge <- xmat %*% B.ridge
    xmat.gar <- cbind(XB.ridge[,seq(1,k,1)]+ XB.ridge[,seq(k+1, 2*k,1)], XB.ridge[,seq(2*k+1, ncol(XB.ridge),1)])
  
    fit.lasso[[i]] <- glmnet(y=y, x=xmat.gar, family=family, alpha=alpha2, lower.limits=0, standardize=FALSE, nlambda=nl2)
    lambda[i,]<-fit.lasso[[i]]$lambda
    if(all(c(i, ii)==index)){
      xmat.gar.save <- xmat.gar
    }
    for(ii in 1:nl2){
      beta[,nl2*(i-1)+ii]<-c(coef(fit.lasso[[i]])[1,ii], # intercept
                             coef(fit.lasso[[i]])[2:(k+1),ii]*coef(fit1.ridge)[2:(k+1),i], # x
                             coef(fit.lasso[[i]])[2:(k+1),ii]*coef(fit1.ridge)[(k+2):(2*k+1),i], # d
                             coef(fit.lasso[[i]])[(k+2):(k+1+kclin),ii]*coef(fit1.ridge)[(2*k+2):(2*k+1+kclin),i]) # clinical 
      

    }
  } # now we have all nl1*nl2 beta vectors               
  coeffs<-beta[,nl2*(index[1]-1)+index[2]]
  names(coeffs)<-rownames(beta)<-c("(Intercept)", colnames(x), colnames(d), colnames(clinical))

  # Return
  df <- matrix(unlist(lapply(fit.lasso, function(X) X$df)),nrow=nrow(lambda), ncol=ncol(lambda),byrow=TRUE)
  L1norm <- matrix(unlist(lapply(fit.lasso, function(X) apply(coef(X)[-1,],2,sum))),nrow=nrow(lambda), ncol=ncol(lambda),byrow=TRUE)
  index.lambda.L1norm=c(nl1-1,which(L1norm[nl1-1,]==head(L1norm[nl1-1,][L1norm[nl1-1,]>=target.L1norm],1)))  # second smallest lambda1 hard-coded
  lambda.L1norm=c(as.numeric(rownames(lambda)[index.lambda.L1norm[1]]),lambda[index.lambda.L1norm[1],index.lambda.L1norm[2]])  
  lambda.min=c(as.numeric(rownames(lambda))[index[1]], lambda[index[1],index[2]])
  names(lambda.L1norm)<-names(lambda.min)<-c("ridge","garrote")
  
  res<-list(call=match.call(), family=family, lambda=lambda, coefficients=beta, glmnet.fit.ridge=fit1.ridge, glmnet.fit.lasso=fit.lasso, k=k, kclin=kclin, df=df, L1norm=L1norm, cv.pred.err=prederr,
            cv.rsquare=rsquare, se.pred.err=se.prederr, 
            fit=list(xmat=xmat, xmat.gar=xmat.gar, lambda=lambda, lambda.min=lambda.min, 
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


predict_garrote<-function(obj, newdata, lambda="lambda.min", type="link"){
    x<-newdata[["x"]]
    d<-newdata[["d"]]
    clinical <- newdata[["clinical"]]
    xmat<-cbind(x,d,clinical)
    if(lambda=="lambda.min"){
      x<-cbind(1, xmat[,names(obj$fit$coefficients[obj$fit$coefficients !=0])[-1]])
      beta <- obj$fit$coefficients[obj$fit$coefficients !=0]
      linpred <- x %*% beta
    } else if (lambda=="L1norm"){
      x<-cbind(1, xmat[,names(obj$fit$coeffs.L1norm[obj$fit$coeffs.L1norm !=0])[-1]])
      beta <- obj$fit$coeffs.L1norm[obj$fit$coeffs.L1norm !=0]
      linpred <- x %*% beta
    }
      else if (lambda=="all"){
        x <- cbind(1, xmat[,names(obj$fit$coefficients)[-1]])
        linpred <- x %*% obj$coefficients
    }
  if(type=="response" & obj$family=="binomial") linpred <- plogis(linpred)
  return(linpred)
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
  beta.x<-beta.hot[substr(names(beta.hot),1,2)=="x."]
  beta.d<-beta.hot[substr(names(beta.hot),1,2)=="d."]
  set.x<-substr(names(beta.x),3,10)
  set.d<-substr(names(beta.d),3,10)
  set.u<-union(set.x, set.d)

  nam.x<-paste("x.", set.u, sep="")
  nam.d<-paste("d.", set.u, sep="")
  beta.x <-beta[nam.x]
  beta.d <- beta[nam.d]
  beta.x[is.na(beta.x)]<-0
  beta.d[is.na(beta.d)]<-0
  if(order == "x") ord<-order(beta.x)
  if(order == "d") ord<-order(beta.d)
  if(order == "none") ord<-1:length(beta.x)
  beta.x <- beta.x[ord]
  beta.d <- beta.d[ord]
  xrange <- range(beta.x*scale[1], beta.d*scale[2])
  if(plot){
    plot(beta.x*scale[1], 1:length(beta.x), type="o", xlim=xrange, xlab="beta (scaled)", ylab="Coefficient")
    points(beta.d*scale[2], 1:length(beta.d), type="o", lty=2, col="red")
    for(i in 1:length(beta.x)) lines(xrange, c(i,i), lty=3, col="grey")
    legend("bottomright", pch=c("o","o"), lty=c(1,2), col=c("black","red"), legend=c("X","D"))
    abline(v=0, col="grey")
  }
  return(list(beta=cbind(beta.x, beta.d), selected=set.u))
  
}  

  
# ============================== METHODS =======================================

perform_garrote<-function(data.obj, penalties=1, family="gaussian", nlambda=c(11,11), 
                          target.L1norm=20, cv=10, R=2, alpha1=0, alpha2=1){
  # data.obj=data.obj; center.interaction.x=0; scale.interaction.x=1; penalties=1; family="gaussian"; nlambda=c(11,11);
  # target.L1norm=20; cv=c(10,10); outer=2; alpha1=0; alpha2=0.99
  
  x<-data.obj[["x"]]
  d<-data.obj[["d"]]
  y<-data.obj[["y"]]
  clinical<-data.obj[["clinical"]]
  xmat <- as.matrix(cbind(x, d, clinical))  
  
  n<-nrow(x)
  k<-ncol(x)
  kclin <- ncol(clinical)
  
  # get number of lambdas
  nl1 <- nlambda[1]
  nl2 <- nlambda[2]
  
  if(k != ncol(d)) stop("d does not match x in dimensions\n")
  if(n != nrow(d) | n!= length(y)) stop("Not equal sample size in variables\n")
  if(length(penalties)==1) penalties<-rep(penalties, 2*k+kclin)
  
  prederr <- rsquare <- prederr2 <- matrix(0, nl1, nl2)
  for(outer in 1:R){
    folds <- sample(rep(1:cv, ceiling(n/cv)))[1:n]
    for(inner in 1:cv){
      beta<-matrix(0,ncol(xmat)+1,nl1*nl2)
      x.train <- xmat[(1:n)[folds!=inner],]
      x.test <- xmat[(1:n)[folds==inner],]
      y.train <- y[(1:n)[folds!=inner]]
      y.test <- y[(1:n)[folds==inner]]
      
      # (1) Ridge regression
      fit1.ridge <- glmnet(y=y.train, x=x.train, family=family, alpha=alpha1, nlambda=nl1, penalty.factor =penalties)
      
      # (2) Lasso regression with restriction of positive coefficients for non-negative shrinkage factors
      for(i in 1:nl1){
        # X garrote = X*beta_ridge
        B.ridge <- diag(coef(fit1.ridge)[-1,i])
        XB.ridge <- x.train %*% B.ridge
        xmat.gar <- cbind(XB.ridge[,seq(1,k,1)]+ XB.ridge[,seq(k+1, 2*k,1)], XB.ridge[,seq(2*k+1, ncol(XB.ridge),1)])
        
        # Positive Lasso: lower.limits ensures positiveness of coeffs
        fit2.garrote <- glmnet(y=y.train, x=xmat.gar, family=family, alpha=alpha2, lower.limits=0, standardize=FALSE, nlambda=nl2)
        
        # (3) Garrote coefficients
        for(ii in 1:nl2){
          beta[,nl2*(i-1)+ii]<-c(coef(fit2.garrote)[1,ii], # intercept
                                 coef(fit2.garrote)[2:(k+1),ii]*coef(fit1.ridge)[2:(k+1),i], # x
                                 coef(fit2.garrote)[2:(k+1),ii]*coef(fit1.ridge)[(k+2):(2*k+1),i], # d
                                 coef(fit2.garrote)[(k+2):(k+1+kclin),ii]*coef(fit1.ridge)[(2*k+2):(2*k+1+kclin),i]) # clinical 
        }
      } # now we have all nl1*nl2 betas               
      
      # validation: compute prediction error
      for(i in 1:nl1){
        for(ii in 1:nl2){
          yhat.test <- cbind(1,x.test) %*% beta[,nl2*(i-1)+ii]
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
  
  # Final betas for all lmabda combinations and return object (including which lambda is best)
  beta<-matrix(0,ncol(xmat)+1,nl1*nl2)
  lambda <- matrix(0, nl1, nl2)
  # (1) Ridge regression
  fit1.ridge <- glmnet(y=y, x=xmat, family=family, alpha=alpha1, nlambda=nl1, penalty.factor =penalties)
  rownames(lambda) <- fit1.ridge$lambda
  
  # (2) Positive lasso 
  fit.lasso <- vector(mode = "list", length = nl1)
  for(i in 1:nl1){
    # X garrote = X*beta_ridge
    B.ridge <- diag(coef(fit1.ridge)[-1,i])
    XB.ridge <- xmat %*% B.ridge
    xmat.gar <- cbind(XB.ridge[,seq(1,k,1)]+ XB.ridge[,seq(k+1, 2*k,1)], XB.ridge[,seq(2*k+1, ncol(XB.ridge),1)])
    
    fit.lasso[[i]] <- glmnet(y=y, x=xmat.gar, family=family, alpha=alpha2, lower.limits=0, standardize=FALSE, nlambda=nl2)
    lambda[i,]<-fit.lasso[[i]]$lambda
    if(all(c(i, ii)==index)){
      xmat.gar.save <- xmat.gar
    }
    for(ii in 1:nl2){
      beta[,nl2*(i-1)+ii]<-c(coef(fit.lasso[[i]])[1,ii], # intercept
                             coef(fit.lasso[[i]])[2:(k+1),ii]*coef(fit1.ridge)[2:(k+1),i], # x
                             coef(fit.lasso[[i]])[2:(k+1),ii]*coef(fit1.ridge)[(k+2):(2*k+1),i], # d
                             coef(fit.lasso[[i]])[(k+2):(k+1+kclin),ii]*coef(fit1.ridge)[(2*k+2):(2*k+1+kclin),i]) # clinical 
      
      
    }
  } # now we have all nl1*nl2 beta vectors               
  coeffs<-beta[,nl2*(index[1]-1)+index[2]]
  names(coeffs)<-rownames(beta)<-c("(Intercept)", colnames(x), colnames(d), colnames(clinical))
  
  # Return
  df <- matrix(unlist(lapply(fit.lasso, function(X) X$df)),nrow=nrow(lambda), ncol=ncol(lambda),byrow=TRUE)
  L1norm <- matrix(unlist(lapply(fit.lasso, function(X) apply(coef(X)[-1,],2,sum))),nrow=nrow(lambda), ncol=ncol(lambda),byrow=TRUE)
  index.lambda.L1norm=c(nl1-1,which(L1norm[nl1-1,]==head(L1norm[nl1-1,][L1norm[nl1-1,]>=target.L1norm],1)))  # second smallest lambda1 hard-coded
  lambda.L1norm=c(as.numeric(rownames(lambda)[index.lambda.L1norm[1]]),lambda[index.lambda.L1norm[1],index.lambda.L1norm[2]])  
  lambda.min=c(as.numeric(rownames(lambda))[index[1]], lambda[index[1],index[2]])
  names(lambda.L1norm)<-names(lambda.min)<-c("ridge","garrote")
  
  res<-list(call=match.call(), family=family, lambda=lambda, coefficients=beta, glmnet.fit.ridge=fit1.ridge, glmnet.fit.lasso=fit.lasso, k=k, kclin=kclin, df=df, L1norm=L1norm, cv.pred.err=prederr,
            cv.rsquare=rsquare, se.pred.err=se.prederr, 
            fit=list(xmat=xmat, xmat.gar=xmat.gar, lambda=lambda, lambda.min=lambda.min, 
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


predict_garrote<-function(obj, newdata, lambda="lambda.min", type="link"){
  x<-newdata[["x"]]
  d<-newdata[["d"]]
  clinical <- newdata[["clinical"]]
  xmat<-cbind(x,d,clinical)
  if(lambda=="lambda.min"){
    x<-cbind(1, xmat[,names(obj$fit$coefficients[obj$fit$coefficients !=0])[-1]])
    beta <- obj$fit$coefficients[obj$fit$coefficients !=0]
    linpred <- x %*% beta
  } else if (lambda=="L1norm"){
    x<-cbind(1, xmat[,names(obj$fit$coeffs.L1norm[obj$fit$coeffs.L1norm !=0])[-1]])
    beta <- obj$fit$coeffs.L1norm[obj$fit$coeffs.L1norm !=0]
    linpred <- x %*% beta
  }
  else if (lambda=="all"){
    x <- cbind(1, xmat[,names(obj$fit$coefficients)[-1]])
    linpred <- x %*% obj$coefficients
  }
  if(type=="response" & obj$family=="binomial") linpred <- plogis(linpred)
  return(linpred)
}

perform_alasso <- function(data.obj, family="gaussian", nlambda=c(11,11), cv=10, R=2, alpha1=0, alpha2=1){
  
 # data.obj <- list("y" = log2(data_mos$eGFR), "clinical"=cbind("age"=data_mos$Age, "sex"=data_mos$Sex), "x"=xt, "d"=d)
  x<-data.obj[["x"]]
  d<-data.obj[["d"]]
  y<-data.obj[["y"]]
  clinical<-data.obj[["clinical"]]
  xmat <- as.matrix(cbind(x, d, clinical))  
  
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
      y.train <- y[(1:n)[folds!=inner]]
      y.test <- y[(1:n)[folds==inner]]
      
      # (1) Ridge regression
      fit1.ridge <- glmnet(y=y.train, x=x.train, family=family, alpha=alpha1, nlambda=nl1)
      
      # (2) Lasso regression 
      for(i in 1:nl1){
        # Penalty
        b.ridge <- c(coef(fit1.ridge)[-1,i])
        penalty_peptides <- abs(b.ridge[1:k])+ abs(b.ridge[(k+1):(2*k)])
        penalty_clinical <- abs(b.ridge[((2*k+1):length(b.ridge))])
        penalty <- 1 / c(penalty_clinical, rep(penalty_peptides, 2))^1
        
        # Lasso: lower.limits ensures positiveness of coeffs
        fit2.lasso <- glmnet(y=y.train, x=x.train, alpha = alpha2, nlambda = nl2, penalty.factor = penalty)
        b.alasso <- coef(fit2.lasso)
      
        beta[1:(ncol(xmat)+1),(nl2*(i-1)+1):(nl2*(i-1)+nl1)]<-as.matrix(b.alasso)
        
      } # now we have all nl1*nl2 betas               
      
      # validation: compute prediction error
      for(i in 1:nl1){
        for(ii in 1:nl2){
          yhat.test <- cbind(1,x.test) %*% beta[,nl2*(i-1)+ii]
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
  
  # Final betas for all lamda combinations and return object (including which lambda is best)
  beta<-matrix(0,ncol(xmat)+1,nl1*nl2)
  lambda <- matrix(0, nl1, nl2)
  # (1) Ridge regression
  fit1.ridge <- glmnet(y=y, x=xmat, family=family, alpha=alpha1, nlambda=nl1)
  rownames(lambda) <- fit1.ridge$lambda
  
  # (2) Positive lasso 
  fit2.alasso <- vector(mode = "list", length = nl1)
  for(i in 1:nl1){
    # Penalty
    b.ridge <- c(coef(fit1.ridge)[-1,i])
    penalty_peptides <- abs(b.ridge[1:k])+ abs(b.ridge[(k+1):(2*k)])
    penalty_clinical <- abs(b.ridge[((2*k+1):length(b.ridge))])
    penalty <- 1 / c(penalty_clinical, rep(penalty_peptides, 2))^1
    
    fit2.alasso[[i]] <- glmnet(y=y, x=xmat, family=family, alpha=alpha2, nlambda=nl2)
    lambda[i,]<-fit2.alasso[[i]]$lambda
    b.alasso <- coef(fit2.alasso[[i]])
    
    beta[1:(ncol(xmat)+1),(nl2*(i-1)+1):(nl2*(i-1)+nl1)] <- as.matrix(b.alasso)
  } # now we have all nl1*nl2 beta vectors               
  coeffs<-beta[,nl2*(index[1]-1)+index[2]]
  names(coeffs)<-rownames(beta)<-c("(Intercept)", colnames(x), colnames(d), colnames(clinical))
  
  # Return
  df <- matrix(unlist(lapply(fit2.alasso, function(X) X$df)),nrow=nrow(lambda), ncol=ncol(lambda),byrow=TRUE)
  lambda.min=c(as.numeric(rownames(lambda))[index[1]], lambda[index[1],index[2]])
  names(lambda.min)<-c("ridge","alasso")
  
  res<-list(call=match.call(), family=family, lambda=lambda, coefficients=beta, glmnet.fit1.ridge=fit1.ridge, glmnet.fit.alasso=fit2.alasso, 
            k=k, kclin=kclin, df=df,  cv.pred.err=prederr,
            cv.rsquare=rsquare, se.pred.err=se.prederr, 
            fit=list(xmat=xmat, lambda=lambda, lambda.min=lambda.min, 
                     coefficients=coeffs, beta=beta, 
                     index.lambda.min=c(index[1], index[2]), 
                     fitted.values=cbind(1,xmat) %*% beta[,nl2*(index[1]-1)+index[2]]))
  attr(res,"class")<-"adaptive lasso"
  return(res)
}

predict_alasso<-function(obj, newdata, lambda="lambda.min", type="link"){
  x<-newdata[["x"]]
  d<-newdata[["d"]]
  clinical <- newdata[["clinical"]]
  xmat<-cbind(x,d,clinical)
  if(lambda=="lambda.min"){
    x<-cbind(1, xmat[,names(obj$fit$coefficients[obj$fit$coefficients !=0])[-1]])
    beta <- obj$fit$coefficients[obj$fit$coefficients !=0]
    linpred <- x %*% beta
  } else if (lambda=="L1norm"){
    x<-cbind(1, xmat[,names(obj$fit$coeffs.L1norm[obj$fit$coeffs.L1norm !=0])[-1]])
    beta <- obj$fit$coeffs.L1norm[obj$fit$coeffs.L1norm !=0]
    linpred <- x %*% beta
  }
  else if (lambda=="all"){
    x <- cbind(1, xmat[,names(obj$fit$coefficients)[-1]])
    linpred <- x %*% obj$coefficients
  }
  if(type=="response" & obj$family=="binomial") linpred <- plogis(linpred)
  return(linpred)
}

perform_lridge <- function(data.obj, family="gaussian", nlambda=c(11,11), cv=10, R=2, alpha1=1, alpha2=0){
  
  # data.obj=data.obj; family="gaussian"; nlambda=c(11,11); cv=10; R=2; alpha1=1; alpha2=0
  x<-data.obj[["x"]]
  d<-data.obj[["d"]]
  y<-data.obj[["y"]]
  clinical<-data.obj[["clinical"]]
  xmat <- as.matrix(cbind(x, clinical))  
  xmat_full <- as.matrix(cbind(x, d, clinical))  
  
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
      beta<-matrix(0,ncol(xmat_full)+1,nl1*nl2)
      rownames(beta) <- c("(Intercept)", colnames(xmat_full))
      x.train <- xmat[(1:n)[folds!=inner],]
      x.test <- xmat[(1:n)[folds==inner],]
      xf.train <- xmat_full[(1:n)[folds!=inner],]
      xf.test <- xmat_full[(1:n)[folds==inner],]
      y.train <- y[(1:n)[folds!=inner]]
      y.test <- y[(1:n)[folds==inner]]
      
      # (1) Lasso regression
      fit1.lasso <- glmnet(y=y.train, x=x.train, family=family, alpha=alpha1, nlambda=nl1)
      
      # (2) Ridge regression 
      for(i in 1:nl1){
        # Include x and d part of selected vars
        b.lasso <- coef(fit1.lasso)[,i]
        nonzero.coefs <- names(b.lasso)[which(b.lasso!=0)]
        nonzero.peps <- str_remove_all(nonzero.coefs, "x.")[-1]
        if(length(nonzero.peps)>0){
          xs.train <- xf.train[,str_detect(colnames(xf.train), pattern=paste0(nonzero.peps, collapse="|"))]
          fit2.ridge <- glmnet(y=y.train, x=xs.train, alpha = alpha2, nlambda = nl2)
          b.ridge <- coef(fit2.ridge)
          beta[rownames(b.ridge),(nl2*(i-1)+1):(nl2*(i-1)+nl1)]<-as.matrix(b.ridge)          
        }
      } # now we have all nl1*nl2 betas               
      
      # validation: compute prediction error
      for(i in 1:nl1){
        for(ii in 1:nl2){
          yhat.test <- cbind(1,xf.test) %*% beta[,nl2*(i-1)+ii]
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
  beta<-matrix(0,ncol(xmat_full)+1,nl1*nl2)
  rownames(beta) <- c("(Intercept)", colnames(xmat_full))
  lambda <- matrix(0, nl1, nl2)
  # (1) Lasso regression
  fit1.lasso <- glmnet(y=y, x=xmat, family=family, alpha=alpha1, nlambda=nl1)
  rownames(lambda) <- fit1.lasso$lambda
  
  # (2) Ridge regression with selected variables (continuous and binary part)
  fit2.ridge <- vector(mode = "list", length = nl1)
  for(i in 1:nl1){
    # Include x and d part of selected vars
    b.lasso <- coef(fit1.lasso)[,i]
    nonzero.coefs <- names(b.lasso)[which(b.lasso!=0)]
    nonzero.peps <- str_remove_all(nonzero.coefs, "x.")[-1]
    
    if(length(nonzero.peps)>0){
      xmat_sub <- xmat_full[,str_detect(colnames(xmat_full), pattern=paste0(nonzero.peps, collapse="|"))]
      fit2.ridge[[i]] <- glmnet(y=y, x=xmat_sub, alpha = alpha2, nlambda = nl2)
      b.ridge <- coef(fit2.ridge[[i]])
      lambda[i,]<-fit2.ridge[[i]]$lambda
      beta[rownames(b.ridge),(nl2*(i-1)+1):(nl2*(i-1)+nl1)]<-as.matrix(b.ridge)          
    }
    else{
      fit2.ridge[[i]]$df <- rep(NA, nl2)
    }
  } # now we have all nl1*nl2 beta vectors               
  coeffs<-beta[,nl2*(index[1]-1)+index[2]]
  names(coeffs)<-rownames(beta)<-c("(Intercept)", colnames(x), colnames(d), colnames(clinical))
  
  # Return
  df <- matrix(unlist(lapply(fit2.ridge, function(X) X$df)),nrow=nrow(lambda), ncol=ncol(lambda),byrow=TRUE)
  lambda.min=c(as.numeric(rownames(lambda))[index[1]], lambda[index[1],index[2]])
  names(lambda.min)<-c("lasso","ridge")
  
  res<-list(call=match.call(), family=family, lambda=lambda, coefficients=beta, glmnet.fit1.lasso=fit1.lasso, glmnet.fit2.ridge=fit2.ridge, 
            k=k, kclin=kclin, df=df,  cv.pred.err=prederr,
            cv.rsquare=rsquare, se.pred.err=se.prederr, 
            fit=list(xmat=xmat_full, lambda=lambda, lambda.min=lambda.min, 
                     coefficients=coeffs, beta=beta, 
                     index.lambda.min=c(index[1], index[2]), 
                     fitted.values=cbind(1,xmat_full) %*% beta[,nl2*(index[1]-1)+index[2]]))
  attr(res,"class")<-"lasso-ridge"
  return(res)
}

  
  
  
