#' ==============================================================================
#' Author: MG
#' Date: 09/05/2023
#' Info: Setup of parameters, packages, paths, etc.
#' ==============================================================================


# --- Parameter ----
print("1/5 - Start setup configuration:")
set.seed(666)
iter = 500
regenerate_simdata = FALSE
dsgn = "dsgn_1"
xmean = 10
xstd = 1
beta_max = 2                              
p = 200                                  

n <- c(100, 200, 400)           
UDdepen <- c("U", "U=2D", "U=D")             
epslvl <- c("none", "moderate", "high")
propzi <- c(0.25, 0.5, 0.75) 
revzi <- c(FALSE)
struczero <- c(0.33, 0.66)
OGM <- c("A", "B", "C")

# --- General scenario matrix ----
scenarios <- expand.grid(
  iter = iter,
  OGM = OGM,
  n = n,
  p = p,
  dsgn = dsgn,
  beta_max = beta_max,
  UDdepen = UDdepen,
  epslvl = epslvl,
  propzi = propzi,
  revzi = revzi, 
  xmean = xmean,
  xstd = xstd,
  struczero = struczero) %>%
  mutate(beta_max = ifelse(OGM %in% c("A", "C"), 1, beta_max),
         beta_max = ifelse(OGM %in% "B", 2, beta_max)) %>% 
  arrange(OGM, n, p, epslvl) %>%
  filter(!(OGM %in% c("B", "C") & revzi == TRUE))


# --- Generate sim data generator ----
if(regenerate_simdata){
  plan(multisession, workers=2)
  list_design <- future_lapply(1:length(unique(scenarios$xmean)), function(i){
    generate_simdesign(p = scenarios$p[i], xmean = unique(scenarios$xmean)[i], xstd = scenarios$xstd[i])}, future.seed = TRUE)
  names(list_design) <- paste0("dsgn_", 1:length(unique(scenarios$dsgn)))
  saveRDS(list_design, here::here("src", "scenario_setup.rds"))
  plan(sequential)
}else{
  list_design <- readRDS(here::here("src", "scenario_setup.rds"))$list_design
  
}
print("2/5 - Simdata generator finished.")

scenarios_UDdepen <- scenarios %>% 
  mutate(a=NA, check_varXD=NA) %>% 
  relocate(a, .after=UDdepen) %>%
  filter(n==100 & epslvl == "none")
for(i in 1:nrow(scenarios_UDdepen)){
  print(paste0(i, "/", nrow(scenarios_UDdepen)))
  scn = scenarios_UDdepen[i,]
  dsgn = list_design[[scenarios_UDdepen[i,]$dsgn]]
  if(scn$UDdepen=="U"){
    scenarios_UDdepen[i,]$a <- 1
    scenarios_UDdepen[i,]$check_varXD <- NA   
  }else{
    # Compute a
    data.val <- data_generation(dsgn = dsgn, OGM = scn$OGM, n = 100000, p=scn$p, beta_max = scn$beta_max, 
                                a = scn$a, epsstd = 0, xmean = scn$xmean, xstd = scn$xstd,
                                propzi = scn$propzi, revzi = scn$revzi, struczero = scn$struczero)
    
    Xb <- data.val$data_gen$X_struc %*% data.val$data_coef$beta_X
    Db <- data.val$data_gen$D_struc %*% data.val$data_coef$beta_D
    
    a_U2D <- (4 - sqrt(8*(var(Xb) / var(Db)))) / (2*(2 - var(Xb) / var(Db)))
    a_UD <- (2 - sqrt(4*(var(Xb) / var(Db)))) / (2*(1 - var(Xb) / var(Db)))
    scenarios_UDdepen[i,]$a <- ifelse(scn$UDdepen == "U=2D", a_U2D, a_UD)
    scenarios_UDdepen[i,]$check_varXD <- var(scenarios_UDdepen[i,]$a*Xb)/var((1-scenarios_UDdepen[i,]$a)*Db)
  }
}
scenarios <- scenarios %>% 
  full_join(., scenarios_UDdepen[,c("OGM", "p","UDdepen", "propzi", "struczero", "revzi", "a", "check_varXD")], 
            by=c("OGM", "p", "UDdepen", "propzi", "struczero", "revzi"))
print("3/5 - Determine UDdependence finished.")


# --- Determine residual variance ----
scenarios_eps <- scenarios %>% 
  mutate(epsstd=NA, check_R2=NA) %>% 
  relocate(epsstd, .after=epslvl) %>% 
  filter(n==100)
for(i in 1:nrow(scenarios_eps)){
  print(paste0(i, "/", nrow(scenarios_eps)))
  scn = scenarios_eps[i,]
  dsgn = list_design[[scenarios_eps[i,]$dsgn]]
  if(scn$epslvl == "none"){
    scenarios_eps[i,]$epsstd <- 0
    scenarios_eps[i,]$check_R2 <- 1
  }else{
    data.val <- data_generation(dsgn = dsgn, OGM = scn$OGM, n = 100000, p=scn$p, beta_max = scn$beta_max, 
                                a = scn$a, epsstd = 0, xmean = scn$xmean, xstd = scn$xstd,
                                propzi = scn$propzi, revzi = scn$revzi, struczero = scn$struczero)
    
    y_true <- data.val$data_gen$y-data.val$data_gen$eps
    r2 <- ifelse(scn$epslvl == "moderate", 0.6, 0.3)
    var.eps <- var(y_true) *((1 - r2)/r2)
    sd.eps <- sqrt(var.eps)
    
    # Check R2
    data.val <- data_generation(dsgn = dsgn, OGM = scn$OGM, n = 100000, p=scn$p, beta_max = scn$beta_max, 
                                a = scn$a, epsstd = sd.eps, xmean = scn$xmean, xstd = scn$xstd,
                                propzi = scn$propzi, revzi = scn$revzi, struczero = scn$struczero)
    checkR2 <- 1-(var(data.val$data_gen$eps)/var(data.val$data_gen$y))
    scenarios_eps[i,]$epsstd <- sd.eps
    scenarios_eps[i,]$check_R2 <- checkR2
  }
}
scenarios <- scenarios %>% 
  full_join(., scenarios_eps[,c("OGM", "p", "epslvl", "epsstd","UDdepen", "a","propzi", "struczero", "revzi", "a", "check_varXD", "check_R2")], 
            by=c("OGM", "p", "epslvl", "UDdepen", "a","propzi", "struczero", "revzi", "check_varXD")) %>%
  dplyr::select(-a.1) %>%
  relocate(a, .after = UDdepen) %>%
  relocate(epsstd, .after=epslvl)
print("4/5 - Determine residual variance finished.")

# ---- Save ----
setup <- list("scenarios"=scenarios, "list_design"=list_design)
saveRDS(setup, here::here("src", "scenario_setup.rds"))
print("5/5 - Setup saved.")



