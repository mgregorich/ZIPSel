#' ==============================================================================
#' Author: MG
#' Date: 09/05/2023
#' Info: Setup of parameters, packages, paths, etc.
#' ==============================================================================


# --- Parameter
set.seed(666)
iter <- 5
regenerate_simdata <- FALSE

n <- c(100,200,400)         # sample size
p <- c(200,400)             # number of candidate predictors
rhomat <- list(rbind(c(.8,.2), c(.8,.2), c(.8,.2), c(.4,.2)))   # correlation ranges in each group
beta_max <- 5                              # maximum coefficient
a <- c(0.2165, 0.28, 1)                    # medium balance of U and D influence on y
epsstd <- c(1, 2)
prop.nonzero <- 0.5         
sampthresh <- 0.05
scenario <- c("A", "B")

tbl_simdata <- data.frame("p" = p, "dsgns" = paste0("dsgn_", 1:length(p)))

# --- Scenario matrix
scenarios <- expand.grid(
  iter = iter,
  scenario = scenario,
  n = n,
  p = p,
  rhomat = rhomat,
  beta_max = beta_max,
  a = a,
  epsstd = epsstd,
  prop.nonzero = prop.nonzero,
  sampthresh = sampthresh) %>%
  merge(., tbl_simdata, by=c("p")) %>%
  relocate(p, .after=n)  %>% 
  mutate(epsstd = ifelse(scenario == "A" & p == 200 & a == .2165, epsstd * 5.5, epsstd),
         epsstd = ifelse(scenario == "A" & p == 200 & a == .28, epsstd * 6.5, epsstd),
         epsstd = ifelse(scenario == "A" & p == 200 & a == 1, epsstd * 20.25, epsstd),
         epsstd = ifelse(scenario == "A" & p == 400 & a == .2165, epsstd * 10, epsstd),
         epsstd = ifelse(scenario == "A" & p == 400 & a == .28, epsstd * 12, epsstd),
         epsstd = ifelse(scenario == "A" & p == 400 & a == 1, epsstd * 40, epsstd),
         epsstd = ifelse(scenario == "B" & p == 200 & a == .2165, epsstd * 5.525, epsstd),
         epsstd = ifelse(scenario == "B" & p == 200 & a == .28, epsstd * 6.1, epsstd),
         epsstd = ifelse(scenario == "B" & p == 200 & a == 1, epsstd * 15, epsstd),
         epsstd = ifelse(scenario == "B" & p == 400 & a == .2165, epsstd * 4.725, epsstd),
         epsstd = ifelse(scenario == "B" & p == 400 & a == .28, epsstd * 5.35, epsstd),
         epsstd = ifelse(scenario == "B" & p == 400 & a == 1, epsstd * 15, epsstd),
         beta_max = ifelse(scenario %in% "A", beta_max*0.4, beta_max)) %>% 
  arrange(p, a, epsstd)


# --- Data simulation design
if(regenerate_simdata){
  plan(multisession, workers = length(p))
  list_design <- future_lapply(1:length(p), function(x) generate_simdesign(p = p[x]), future.seed = TRUE)
  names(list_design) <- paste0("dsgn_", 1:length(p))
  plan(sequential)  
}else{
  list_design <- readRDS(here::here("src", "scenario_setup.rds"))$list_design
}

setup <- list("scenarios"=scenarios, "list_design"=list_design)
saveRDS(setup, here::here("src", "scenario_setup.rds"))



