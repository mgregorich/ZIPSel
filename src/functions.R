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
  mae.val <- mean(pred-obs)
  cil.val <- mean(obs)- mean(pred)
  
  res <- data.frame(R2 = r2,
                    RMSE = rmse.val,
                    MAE = mae.val,
                    C = Cind,
                    CalbSlope=CS,
                    CalbinLarge = cil.val) 
  return(res)
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
# Code by Georg Heinze

protogarrote<-function(data.obj, center.interaction.x=0, scale.interaction.x=1, penalties=1, family="gaussian", nlambda=c(11,11), 
                       target.L1norm=20, cv=10, outer=2, alpha1=0, alpha2=0.99){
  # data.obj=data.obj; center.interaction.x=0; scale.interaction.x=1; penalties=1; family="gaussian"; nlambda=c(11,11);
  # target.L1norm=20; cv=10; outer=2; alpha1=0; alpha2=0.99

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
        # lambda1<-lambda$lambda1[i]
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
  beta.x<-beta.hot[substr(names(beta.hot),1,2)=="x."]
  beta.d<-beta.hot[substr(names(beta.hot),1,2)=="d."]
  set.x<-substr(names(beta.x),3,10)
  set.d<-substr(names(beta.d),3,10)
  set.u<-union(set.x, set.d)
  if(!obj$fit.int){
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
  if(obj$fit.int){
    nam.ix<-paste("i*x.", set.u, sep="")
    nam.id<-paste("i*d.", set.u, sep="")
    beta.ix <- beta[nam.ix]
    beta.id <- beta[nam.id]
#    set.ix<-substr(names(beta.ix),5,12)
#    set.id<-substr(names(beta.id),5,12)
#    set.u<-union(union(union(set.x, set.d), set.ix), set.id)
    nam.x<-paste("x.", set.u, sep="")
    nam.d<-paste("d.", set.u, sep="")
    nam.ix<-paste("i*x.", set.u, sep="")
    nam.id<-paste("i*d.", set.u, sep="")
    beta.x <-beta[nam.x]
    beta.d <- beta[nam.d]
    beta.ix <-beta[nam.ix]
    beta.id <- beta[nam.id]
    beta.x[is.na(beta.x)]<-0
    beta.d[is.na(beta.d)]<-0
    beta.ix[is.na(beta.ix)]<-0
    beta.id[is.na(beta.id)]<-0
    if(order == "x") ord<-order(beta.x)
    if(order == "d") ord<-order(beta.d)
    if(order =="i*x") ord<-order(beta.ix)
    if(order =="i*d") ord<-order(beta.id)
    if(order == "none") ord<-1:length(beta.x)
    beta.x <- beta.x[ord]
    beta.d <- beta.d[ord]
    beta.ix <- beta.ix[ord]
    beta.id <- beta.id[ord]
    xrange <- range(beta.x*scale[1], beta.d*scale[2], beta.ix*scale[3], beta.id*scale[4])
    if(plot){
      plot(beta.x*scale[1], 1:length(beta.x), type="o", xlim=xrange, xlab="beta (scaled)", ylab="Coefficient")
      points(beta.d*scale[2], 1:length(beta.d), type="o", lty=2, col="red")
      points(beta.ix*scale[3], 1:length(beta.ix), type="o", lty=2, col="blue")
      points(beta.id*scale[4], 1:length(beta.id), type="o", lty=2, col="green")
      for(i in 1:length(beta.x)) lines(xrange, c(i,i), lty=3, col="grey")
      legend("bottomright", pch=c("o","o"), lty=c(1,2), col=c("black","red", "blue", "green"), legend=c("X","D", "I*X", "I*D"))
      abline(v=0, col="grey")
    }
    return(list(beta=cbind(beta.x, beta.d, beta.ix, beta.id), selected=set.u))
  }
}  

  


  
  
