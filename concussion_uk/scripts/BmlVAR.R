library(dplyr)
library(imputeTS)

# Load in data ---------------------------------------------

raw_data <- haven::read_sav("data/PCSNetworkModelAnalysis.sav") |>
  janitor::clean_names()

# 01 Create a new dataset with NAs ----------------------------------------

data <- raw_data |>
  dplyr::filter(survey_day != 0)

all_combinations <- expand.grid(
  survey_day = 1:30,
  survey_time = c(10, 3, 8)
)

symptoms <- names(raw_data)[19:length(names(raw_data))]

# Ensure all combinations for each ID are present, filling missing values with NA
df_complete <- data |>
  dplyr::group_by(subject_id) |>
  tidyr::complete(survey_day = all_combinations$survey_day, survey_time = all_combinations$survey_time) |>
  dplyr::ungroup()


# 02 Impute data ---------------------------------------------------------
new_param <- list(maxit = 5000)
data_imputed  <- df_complete |>
  group_by(subject_id) |>
  mutate(na_count = sum(!is.na(head))) |>
  dplyr::filter(na_count >= 3, subject_id != 262) |>
  mutate(na_count = sum(!is.na(head))) |>
  mutate(across(symptoms,
                ~imputeTS::na_kalman(.x, optim.control = new_param),
                .names = "{.col}_imp")) |>
  ungroup()

convergence_warnings <- dplyr::last_dplyr_warnings(n = 78)

data_imputed2 <- data_imputed |>
  dplyr::select(subject_id, dplyr::contains("_imp")) |>
  rename_with(~stringr::str_remove(.x, "_imp"), contains("_imp"))

# 03 Scale data --------------------------------------------------------

df_data_scaled <- data_imputed2 |>
  group_by(subject_id) |>
  mutate(across(everything(), ~scale(.x))) |>
  ungroup() |> # transform NaN into zeros
  dplyr::mutate(
         dplyr::across(dplyr::everything(),
                       ~tidyr::replace_na(.x,  0)))

## Add small noise

set.seed(999)
for (j in 2:ncol(df_data_scaled)) {
  for (i in 1:nrow(df_data_scaled)){
    noise <- rnorm(1, 0, 0.001)
    df_data_scaled[i, j] <- df_data_scaled[i, j] + noise
  }
}

data_list_scaled <- split(df_data_scaled, ~subject_id)

df_data_centered <- data_imputed2 |>
  group_by(subject_id) |>
  mutate(across(everything(), ~scale(.x, scale = FALSE))) |>
  ungroup()

data_list_centered <- split(df_data_centered, ~subject_id)

# BmlVAR ---------------------------------------------------------------

n_var <- 22

# indicators for partial correlations
idx_rho <- upper.tri(matrix(1, n_var, n_var, byrow = F)) |>
      c() |>
      which()

# repeat the outcome because we have two single-predictor regressions
## reg_data <- cbind(
##     neuroticism_vec_scaled,
##     neuroticism_vec_scaled
##   )

Y_scaled <- df_data_scaled |>
    dplyr::select(-c("subject_id")) |>
    as.matrix()

Y <- data_imputed2 |>
    dplyr::select(-c("subject_id")) |>
    as.matrix()

# obtain number of time points per person
n_tp <- data_imputed2 |>
  group_by(subject_id) |>
  count() |>
  pull(n)

n_id <- length(unique(data_imputed2$subject_id))

# Prepare stan data
stan_data_scaled <- list(
    K = n_var,
    I = n_id,
    P = 2,   # number of regression outcomes
    N_total = nrow(data_bringmann_imputed),
    n_t = n_tp,
    n_pc = n_var * (n_var - 1) / 2,
    idx_rho = array(idx_rho, dim = length(idx_rho)),
    Y = Y_scaled,
    reg_covariate = reg_data,
    sparsity = 2,
    mean_center = 0
  )

# perform mean-centering within Stan
stan_data <- list(
    K = n_var,
    I = n_id,
    P = 2,   # number of regression outcomes
    N_total = nrow(data_bringmann_imputed),
    n_t = n_tp,
    n_pc = n_var * (n_var - 1) / 2,
    idx_rho = array(idx_rho, dim = length(idx_rho)),
    Y = Y/20,
    reg_covariate = reg_data,
    sparsity = 2,
    mean_center = 1
  )

## fit-bmlvar-----------------------------------------------------------

model_name <- "mlvar_lkj"

stan_model <-   rstan::stan_model(
    file = here::here("scripts", paste0(model_name, ".stan")),
    model_name = model_name
  )

fit_bmlvar_centered <- rstan::sampling(
  object = stan_model,
  pars = c("Beta_raw", "Intercepts_raw", "L_Theta", "Rho", "lp__"),
  include = FALSE,
  data = stan_data,
  chains = 4,
  cores = 4,
  warmup =500,
  iter = 2500,
  init = 0,
  control = list(adapt_delta = 0.99),
  verbose = FALSE
  )
