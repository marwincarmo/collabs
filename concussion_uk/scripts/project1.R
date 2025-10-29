
# 0 Load required packages ------------------------------------------------

library(mlVAR)
library(imputeTS)
library(dplyr)
library(ggplot2)
library(stringr)
library(rstan)

# 1 Load in data ----------------------------------------------------------

d <- readRDS("data/SynthData.rds")
head(d$syn)

# 2 Select the most central symptoms --------------------------------------


sxd <- d$syn |> 
  dplyr::select(Main.SubjectIDNum, obs_num, Headache, PressureInHead, FeelingSlowedDown, FeelingLikeInAFog,
                DontFeelRight, DifficultyConcentrating, FatigueOrLowEnergy)


# 3 Impute missing data ---------------------------------------------------

## 3.1 Check missing pattern

sxd |> 
  group_by(Main.SubjectIDNum) |> 
  summarise(across(everything(),
                   ~sum(is.na(.x))/n())) |> 
  summarise(across(everything(),
                   ~mean(.x, na.rm = TRUE))) |> 
  # take average across all variables
  select(!Main.SubjectIDNum) |> 
  pull() |>
  mean()

## 3.2 Imput data with inputeTS()

sxd_imp <- sxd |> 
  group_by(Main.SubjectIDNum) |> 
  mutate(across(Headache:FatigueOrLowEnergy,
                ~imputeTS::na_kalman(.x, optim.control = 5000),
                .names = "{.col}_imp")) |> 
  ungroup() |> 
  select(Main.SubjectIDNum, contains("_imp")) |>
  rename_with(~str_remove(.x, "_imp"), contains("_imp")) #|>
  #mutate(id = as.numeric(as.factor(Main.SubjectIDNum)))


sxd_imp_std <- sxd_imp |> 
  dplyr::mutate(
           dplyr::across(Headache:DontFeelRight, ~scale(.x))
         ) |>
  dplyr::ungroup()

## sxd_imp %>%
##   group_by(Main.SubjectIDNum) %>%
##   summarize(has_duplicates = any(duplicated(obs_num))) |>
##   dplyr::pull(has_duplicates) |>
##   sum()

## Look for zero variance
sxd_imp |> 
  group_by(Main.SubjectIDNum) |> 
  summarise(across(everything(), sd)) |> 
  filter(if_any(everything(), ~ . == 0 | is.na(.)))
  print(n = Inf)

  sxd_imp <- sxd_imp |> 
    dplyr::filter(!Main.SubjectIDNum %in% c(576, 853, 2586))

# 4 mlVAR -----------------------------------------------------------------

vars <- colnames(sxd_imp[,-c(1)])

fit_mlvar <- mlVAR(data=sxd_imp, vars = vars, idvar = "id",
                   contemporaneous = "orthogonal", temporal = "orthogonal", estimator = "lmer")

summary(fit_mlvar)

# 5 Obtain centrality measures --------------------------------------------

# Obtain centrality of fitted mlVAR object
centrality_mlvar <- function(fit,
                             ignore_ar = TRUE){
  
  #--- Prepare
  n_var <- length(fit$fit$var)
  n_id <- length(fit$IDs)
  n_var_temp <- n_var
  
  #--- Obtain networks
  l_beta <- lapply(1:n_id, function(i){
    mlVAR::getNet(fit,
                  subject = i,
                  type = "temporal",
                  nonsig = "show")
  })
  l_pcor <- lapply(1:n_id, function(i){
    mlVAR::getNet(fit,
                  subject = i,
                  type = "contemporaneous",
                  nonsig = "show")
  })
  
  # Obtain overall adjacency matrix
  adj_beta <- ifelse(
    mlVAR::getNet(fit,
                  type = "temporal",
                  nonsig = "hide") != 0, 1, 0)
  adj_pcor <- ifelse(
    mlVAR::getNet(fit,
                  type = "contemporaneous",
                  nonsig = "hide") != 0, 1, 0) 
  
  # Set zero fixed effects to zero
  l_beta <- lapply(l_beta, function(x){
    x * adj_beta
  })
  l_pcor <- lapply(l_pcor, function(x){
    x * adj_pcor
  })
  
  
  #--- Density
  dens_temp <- lapply(l_beta, function(x){
    sum(abs(x))/(n_var^2)
  })
  dens_cont <- lapply(l_pcor, function(x){
    x <- x
    diag(x) <- 0
    x[lower.tri(x)] <- 0L
    sum(abs(x))/(n_var * (n_var-1)/2)
  })
  
  
  #--- Centrality
  # Important: in mlVAR, the lagged vars are rows, not columns
  # but only for the getNet function!
  outstrength <- lapply(l_beta, function(x){
    if(isTRUE(ignore_ar)){
      diag(x) <- 0
      n_var_temp <- n_var - 1
    }
    rowSums(abs(x))/n_var_temp
  })
  instrength <- lapply(l_beta, function(x){
    if(isTRUE(ignore_ar)){
      diag(x) <- 0
      n_var_temp <- n_var - 1
    }
    colSums(abs(x))/n_var_temp
  })
  strength <- lapply(l_pcor, function(x){
    # because diagonal is always zero, we subtract 1
    colSums(abs(x))/(n_var - 1)
  })
  
  # return list
  return(list(outstrength = outstrength,
              instrength = instrength,
              strength = strength,
              dens_temp = dens_temp,
              dens_cont = dens_cont))
  
}

centrality_res <- centrality_mlvar(fit_mlvar)


# 6 BmlVAR --------------------------------------------

vars <- colnames(sxd_imp[,-c(1)])

stan_model <- rstan::stan_model(file = "scripts/mlvar_lkj.stan",
                                model_name = "mlvar_lkj")

n_var <- length(vars)

# indicators for partial correlations
idx_rho <- upper.tri(matrix(1, n_var, n_var, byrow = F)) |>
      c() |>
      which()

Y_scaled <- sxd_imp_std |>
    dplyr::select(-c("Main.SubjectIDNum")) |>
    as.matrix()

Y <- sxd_imp_std |>
    dplyr::select(-c("Main.SubjectIDNum")) |>
    as.matrix()

# obtain number of time points per person
n_tp <- sxd_imp_std |>
  group_by(Main.SubjectIDNum) |>
  count() |>
  pull(n)

n_id <- length(unique(sxd_imp_std$Main.SubjectIDNum))

# Prepare stan data
stan_data_scaled <- list(
    K = n_var,
    I = n_id,
    N_total = nrow(sxd_imp_std),
    n_t = n_tp,
    n_pc = n_var * (n_var - 1) / 2,
    idx_rho = array(idx_rho, dim = length(idx_rho)),
    Y = Y_scaled,
    sparsity = 2,
    mean_center = 0
)

# perform mean-centering within Stan
stan_data <- list(
    K = n_var,
    I = n_id,
    N_total = nrow(sxd_imp_std),
    n_t = n_tp,
    n_pc = n_var * (n_var - 1) / 2,
    idx_rho = array(idx_rho, dim = length(idx_rho)),
    Y = Y/20,
    sparsity = 2,
    mean_center = 1
)

fit_bmlvar_centered <- rstan::sampling(
    object = stan_model,
    pars = c("Beta", "Sigma", "Rho", "Beta_in_strength", "Beta_out_strength", "Rho_strength", "Beta_density", "Rho_density"), # Explicitly include desired outputs
    include = TRUE, # Now we are including *only* these
    data = stan_data,
    chains = 4,
    cores = 8,
    warmup =500,
    iter = 2500,
    init = 0,
    control = list(adapt_delta = 0.99),
    verbose = FALSE
)
