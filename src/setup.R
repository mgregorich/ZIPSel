#' ==============================================================================
#' Author: MG
#' Date: 09/05/2023
#' Info: Setup of parameters, packages, paths, etc.
#' ==============================================================================


# Parameter
set.seed(666)
iter = 2
n <- c(100,200,400)         # sample size
p <- c(200,400)             # number of candidate predictors
rhomat <- list(rbind(c(.8,.2), c(.8,.2), c(.8,.2), c(.4,.2)))   # correlation ranges in each group
beta_max <- 5               # maximum coefficient
a <- 0.5                    # controls balance of U and D influence on y
epsstd <- 2.5
prop.nonzero <- 0.5         
sampthresh <- 0.05

# Scenario matrix
scenarios <- expand.grid(
  iter = iter,
  n = n,
  p = p,
  rhomat = rhomat,
  beta_max = beta_max,
  a = a,
  epsstd = epsstd,
  prop.nonzero = prop.nonzero,
  sampthresh = sampthresh) 



# Data simulation design
xmean = 0
xstd = 0.5
data_design <- scenarios %>% 
  data.frame() %>%
  dplyr::select(p) %>%
  filter(!duplicated(p))
list_design <- list()
data_design$dsgn <- NA
for(i in 1: nrow(data_design)){
  print(paste0("Sim design: ",i, "/", nrow(data_design)))
  # Groupwise Hub correlation design filled with toeplitz
  hub_cormat <- simcor.H(k = 4, size = rep(data_design[i,]$p/4,4),
                         rho = rhomat[[1]], power = 1, epsilon = 0.075, eidim = 2)
  if(!matrixcalc::is.positive.definite(hub_cormat)) hub_cormat <- nearPD(hub_cormat, base.matrix = TRUE, keepDiag = TRUE)$mat
  
  # Data design
  distlist <- rep(list(partial(function(x, meanlog, sdlog) qlnorm(x, meanlog = meanlog, sdlog = sdlog),
                               meanlog = xmean, sdlog = xstd)), nrow(hub_cormat)) 
  dsgn = simdata::simdesign_norta(cor_target_final = hub_cormat, dist = distlist, 
                                  transform_initial = data.frame,
                                  names_final = paste0("V",1:nrow(hub_cormat)), seed_initial = 1) 
  list_design[[i]] <- dsgn
  data_design[i,]$dsgn <- i
}

scenarios <- merge(scenarios, data_design, by=c("p")) %>%
  relocate(p, .after=n)
setup <- list("scenarios"=scenarios, "list_design"=list_design)
saveRDS(setup, here::here("src", "scenario_setup.rds"))



