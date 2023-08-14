rm(list=ls())

# Packages
pacman::p_load(simdata)

# ===============================================================================
# Functions for correlation matrix (not relevant for the question)
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
# =============================================================================


# =============================================================
# Example for simdata question 
# =============================================================

# === Define simdata design
# Parameters 
xmean <- 0
xstd <- 0.5
p <- 100
rhomat <- rbind(c(.8,.2), c(.8,.2), c(.8,.2), c(.4,.2))

# Correlation matrix
hub_cormat <- simcor.H(k = 4, size = rep(p/4,4), rho = rhomat, power = 1, epsilon = 0.075, eidim = 2)
if(!matrixcalc::is.positive.definite(hub_cormat)) hub_cormat <- nearPD(hub_cormat, base.matrix = TRUE, keepDiag = TRUE)$mat


# Simdesign
partial <- function(f, ...) {
  f_args <- list(...)
  
  function(...) {
    do.call(f, c(f_args, list(...)))
  }
}
distlist <- rep(list(partial(function(x, meanlog, sdlog) qlnorm(x, meanlog = meanlog, sdlog = sdlog),
                             meanlog = xmean, sdlog = xstd)), nrow(hub_cormat))
dsgn = simdata::simdesign_norta(cor_target_final = hub_cormat, dist = distlist, 
                                transform_initial = data.frame,
                                names_final = paste0("V",1:nrow(hub_cormat)), seed_initial = 1) 


# ==== Simdesign works when xmean and xstd specified in global environment
simulate_logdata <- function(n, xmean, xstd){
  data_logvars = simdata::simulate_data(dsgn, n, seed = 2)
  return(data_logvars)
}

xmean <- 0
xstd <- 0
n <- 100
simulate_logdata(n, xmean, xstd)

# ==== Simdata does not work when specified in function?
rm(list=setdiff(ls(), c("dsgn", "simulate_logdata")))

scn <- data.frame(n =100, xmean = 0, xstd = 0.5)

simulate_logdata(n=scn$n, xmean=scn$xmean, xstd=scn$xstd)

#' Frage: wie kann ich in meiner Funktion simulate_logdata simulate_data abändern,
#' sodass die Funktion simdata::simulate_data() die Parameter hat?