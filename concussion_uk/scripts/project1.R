
# 0 Load required packages ------------------------------------------------

library(mlVAR)
library(imputeTS)
library(dplyr)
library(ggplot2)
library(stringr)

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
  select(Main.SubjectIDNum, obs_num, contains("_imp")) |>
  rename_with(~str_remove(.x, "_imp"), contains("_imp")) |> 
  mutate(id = as.numeric(as.factor(Main.SubjectIDNum)))






sxd_imp_std <- sxd_imp |> 
  dplyr::mutate(
    dplyr::across(Headache:DontFeelRight, scale)
  )

sxd_imp %>%
  group_by(Main.SubjectIDNum) %>%
  summarize(has_duplicates = any(duplicated(obs_num))) |> 
  dplyr::pull(has_duplicates) |> 
  sum()

## Look for zero variance
sxd_imp |> 
  group_by(Main.SubjectIDNum) |> 
  summarise(across(everything(), sd)) |> 
  filter(if_any(everything(), ~ . == 0 | is.na(.)))
  print(n = Inf)

  sxd_imp <- sxd_imp |> 
    dplyr::filter(!Main.SubjectIDNum %in% c(576, 853, 2586))

# 4 mlVAR -----------------------------------------------------------------

vars <- colnames(sxd_imp[,-c(1,2)])
sxd_imp$obs_num <- 

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



# Diagnostics -------------------------------------------------------------

library(dplyr)

# Define the variables in the model
vars_in_model <- vars

# Manually create a lagged dataset to see how many complete cases exist for a model
usable_rows_count <- sxd_imp %>%
  select(Main.SubjectIDNum, obs_num, all_of(vars_in_model)) %>%
  arrange(Main.SubjectIDNum, obs_num) %>%
  group_by(Main.SubjectIDNum) %>%
  # Create lagged versions of all variables to simulate the predictors
  mutate(across(all_of(vars_in_model), ~lag(.), .names = "{.col}_lag")) %>%
  ungroup() %>%
  # Drop the ID and time variables for counting
  select(-Main.SubjectIDNum, -obs_num) %>%
  # How many rows are complete across all original AND lagged variables?
  filter(complete.cases(.)) %>%
  nrow()

cat("Total number of usable rows for the temporal models:", usable_rows_count, "\n")

library(dplyr)

# Define your study's structure
beeps_per_day <- 3
total_days <- 7
required_obs <- beeps_per_day * total_days # Should be 21

# Create day/beep columns and filter for complete subjects
sxd_filtered_complete <- sxd %>%
  mutate(
    day  = ((obs_num - 1) %/% beeps_per_day) + 1,
    beep = ((obs_num - 1) %% beeps_per_day) + 1
  ) %>%
  group_by(Main.SubjectIDNum) %>%
  filter(n() == required_obs) %>%
  ungroup()

# --- Verification ---
# You can check the results to see how many subjects were kept

# Original number of subjects
cat("Original number of subjects:", n_distinct(sxd$Main.SubjectIDNum), "\n")

# Number of subjects with complete data (21 observations)
cat("Subjects with complete data:", n_distinct(sxd_filtered_complete$Main.SubjectIDNum), "\n")

plot(fit_mlvar, type = "temporal")
