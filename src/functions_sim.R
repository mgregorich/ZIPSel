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
data_generation <- function(n, p){
  
  # Remove later
  # n<-300 
  # p<-100 
  
  # Clinical variables
  age <- floor(runif(n, 20, 81))
  sex <- rbinom(n, size=1, prob=0.5)
  
  # Zero-inflated predictors
  prop.nonzero<-runif(p,0.75, 1)   # proportion nonzero for each peptide
  prop.nonzero<-sort(prop.nonzero)
  means.peptide<-2+prop.nonzero*2    # we create consonant peptides only
  X <- matrix(rnorm(n*p), n, p, byrow=TRUE) + matrix(means.peptide, n, p, byrow=TRUE)
  X.zi <- peptide * rbinom(n*p, size=1, prob=matrix(prop.nonzero, n, p, byrow=TRUE))

  # Predictive X
  predictiveX <- sample(1:p, size=10)
  b0 <- 1
  bX <- 2
  bsex <- 1.5
  bage <- -2
  b.true <- c(b0, bsex, bage, bX)
  
  # Outcome
  eps <- rnorm(n, mean=0, sd=5)
  y <- b0 + rowSums(2* X[,predictiveX] + bsex*sex + bage*age/100) + eps
  
  out <- data.frame("y"=y, "sex"=sex, "age"=age, "x"=X, "xz"=X.zi)
  return(out)
}


# ============================ Data analysis =================================
data_analysis <- function(df){
  
  # Prepare data
  
  
  # ridge-garrote
  
  
  # ridge-lasso
  
  
  # lasso-ridge
  
  
  # network-based group penalization
  
  
  # multi univariable
  
  
  # Merge results
  
  return()
}