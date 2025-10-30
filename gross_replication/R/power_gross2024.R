## 0 Load packages ---------------------------------------

library(lme4)
library(dplyr)
library(readr)
library(mediation)

## 1 Load and prepare data ------------------------------

esm <- readr::read_csv("ESM.Clean.Nov2022.csv")

## person-mean centering
esm_pc <- esm |>
  ## extracted from the original analysis
  ## create inner speach variables
  dplyr::mutate(clarity_all = coalesce(clarity, no_inner_speech_clarity),
                valence_all = coalesce(valence, no_inner_speech_valence),
                reaction_all = coalesce(reaction, no_inner_speech_reaction),
                interesting_all = coalesce(interesting, no_inner_speech_interesting),
                ## center variables
                wb = wb - 4,
                valence = valence - 4,
                valence_all = valence_all - 4,
                no_inner_speech_valence = no_inner_speech_valence - 4,
                ## recode activity
                activity = dplyr::case_when(
                                    activity == 1 ~ "Social",
                                    activity == 2 ~ "Physical",
                                    activity == 3 ~ "Restful",
                                    activity == 4 ~ "Household",
                                    activity == 5 ~ "Cognitive",
                                    activity == 6 ~ "Other",
                                    .default = NA
                                  ),
                attention  = ifelse(attention == 1, 1, -1)
                #attention = ifelse(attention == 1, "Present", "Not Present")
                ) |>
  ## Person-mean centering (from chunk 'Person mean centering')
  dplyr::group_by(PID) |>
  dplyr::mutate(
           wb_mean2 = mean(wb, na.rm = T),
           clarity_all_mean2 = mean(clarity_all, na.rm = T),
           valence_all_mean2 = mean(valence_all, na.rm = T),
           reaction_all_mean2 = mean(reaction_all, na.rm = T),
           interesting_all_mean2 = mean(interesting_all, na.rm = T)
         ) |>
  dplyr::ungroup() |>
  ## Within-person centered variables
  dplyr::mutate(
           wb_cw2 = wb - wb_mean2,
           clarity_all_cw2 = clarity_all - clarity_all_mean2,
           valence_all_cw2 = valence_all - valence_all_mean2,
           reaction_all_cw2 = reaction_all - reaction_all_mean2,
           interesting_all_cw2 = interesting_all - interesting_all_mean2)

## Create the "Inner Speech" only dataset (from chunk 'Inner speech datasets')
esm_inner <- esm_pc |>
  dplyr::filter(inner_speech == 0) # Inner Speech as they have on the paper

## Convert attention to a factor for the mediation package (from chunk 'With covariates')
## Note: The paper and script chunk 'As individual MLMs' use -1/1 coding.
## But the final mediation chunk 'With covariates' uses as.factor(). We MUST follow
## the final mediation chunk for our power analysis.
esm_inner$attention <- as.factor(esm_inner$attention)

## 3 Fit the two key models from the script

## Main hypothesis: Is Thought Valence more negative when people
## are mind wandering versus present?

## Table 3 - left panel

## These are from the chunk 'With covariates' under '##Model 2: Mediation of Valeince'

model1a_fit <- lmer(formula = wb_cw2 ~ activity + clarity_all_cw2+ #as a covariate
                    interesting_all_cw2 + attention+ (1|PID),
                    data = esm_inner,
                   REML = FALSE) #must set this to F to return loglik and AIC
summary(model1a_fit)
## Model a: Mediator model (Valence ~ Attention + Covariates)
## (Corresponds to 'fit.mediator1' in the script)
## Does Attention State predict  Thought Valence (Mediator)?

fit_mediator <- lme4::lmer(
  valence_all_cw2 ~ attention + activity + clarity_all_cw2 + interesting_all_cw2 + (1 | PID),
  data = esm_inner
)
summary(fit_mediator)
## Model b: Outcome model (Mood ~ Attention + Valence + Covariates)
## (Corresponds to 'fit.dv1' in the script)
fit_outcome <- lme4::lmer(
                       wb ~ attention + valence_all_cw2 +
                         activity + clarity_all_cw2 + interesting_all_cw2 + (1 | PID),
  data = esm_inner
)

# --- 0. Load required libraries ---
library(lme4)
library(mediation)
library(dplyr)
library(purrr)     # For iteration
library(ggplot2)   # For plotting

# --- !!! IMPORTANT !!! ---
# You MUST have the 'esm_inner' dataframe from our Step 1 script
# in your R environment for this to work. The simulation
# samples from it to create realistic covariates.
#
# If you don't have it, please run this part of your script again:
#
# esm_pc <- ... (all your dplyr data prep)
# esm_inner <- esm_pc |> dplyr::filter(inner_speech == 0)
# esm_inner$attention <- as.factor(esm_inner$attention)
#
# ---

# --- 1. Define Ground Truth Parameters ---
# (Extracted from your summary() outputs)

# Model A (Mediator) Fixed Effects
b_med <- list(
  "(Intercept)" = -0.28529,
  attention1 = 0.11700,
  activityHousehold = 0.09521,
  activityOther = -0.14934,
  activityPhysical = 0.21010,
  activityRestful = 0.09941,
  activitySocial = 0.21218,
  clarity_all_cw2 = 0.08972,
  interesting_all_cw2 = 0.15971
)
# Model A Residual SD
sd_med_resid <- 1.064

# Model B (Outcome) Fixed Effects
b_out <- list(
  "(Intercept)" = -0.227307,
  attention1 = 0.117671,
  valence_all_cw2 = 0.597218,
  activityHousehold = 0.059135,
  activityOther = 0.004509,
  activityPhysical = 0.102689,
  activityRestful = 0.116552,
  activitySocial = 0.179922,
  clarity_all_cw2 = 0.025502,
  interesting_all_cw2 = 0.071591
)
# Model B Residual SD
sd_out_resid <- 0.9444

# Sample structure from original data
j_prompts <- 13 # Avg prompts (13.2 in original, we'll use 13)


# --- 2. Define the Simulation Function ---

#' Runs a single iteration of the mediation simulation
#'
#' @param N Number of participants to simulate
#' @param j Number of prompts per participant
#' @param esm_inner_orig The original `esm_inner` dataframe to sample covariates from
#'
#' @return A data.frame with the p-value for the ACME (indirect effect)
run_sim_iteration <- function(N, j, esm_inner_orig) {

  # 1. Create participant and prompt structure
  sim_data <- data.frame(
    PID = rep(1:N, each = j)
  )
  n_total_obs <- N * j

  # 2. Simulate covariates by resampling from the original data
  # This is the easiest way to preserve their distributions and correlations
  covariate_indices <- sample(1:nrow(esm_inner_orig), size = n_total_obs, replace = TRUE)

  sim_data$attention <- esm_inner_orig$attention[covariate_indices]
  sim_data$activity <- esm_inner_orig$activity[covariate_indices]
  sim_data$clarity_all_cw2 <- esm_inner_orig$clarity_all_cw2[covariate_indices]
  sim_data$interesting_all_cw2 <- esm_inner_orig$interesting_all_cw2[covariate_indices]

  # We must re-level the factor to match the model's coefficients
  sim_data$activity <- factor(sim_data$activity,
                              levels = c("Cognitive", "Household", "Other",
                                         "Physical", "Restful", "Social"))

  # 3. Simulate Mediator (Valence)
  # Create a model matrix to easily calculate predicted values
  X_med <- model.matrix(~ attention + activity + clarity_all_cw2 + interesting_all_cw2, data = sim_data)

  # Get the beta coefficients in the correct order
  b_med_vec <- unlist(b_med[colnames(X_med)])

  # Calculate predicted mediator value (mu)
  mu_med <- X_med %*% b_med_vec

  # Add residual error
  sim_data$sim_med <- mu_med + rnorm(n_total_obs, mean = 0, sd = sd_med_resid)

  # 4. Simulate Outcome (Mood)
  # Create model matrix for the outcome model
  X_out <- model.matrix(~ attention + sim_med + activity + clarity_all_cw2 + interesting_all_cw2, data = sim_data)

  # Get betas in the correct order (note 'sim_med' replaces 'valence_all_cw2')
  b_out_vec <- unlist(b_out[colnames(X_out)])
  names(b_out_vec)[names(b_out_vec) == "valence_all_cw2"] <- "sim_med" # align names

  # Calculate predicted outcome value (mu)
  mu_out <- X_out %*% b_out_vec

  # Add residual error
  sim_data$sim_out <- mu_out + rnorm(n_total_obs, mean = 0, sd = sd_out_resid)

  # 5. Run the authors' exact (singular) analysis
  # We use lmer() even though Var=0 to perfectly replicate their workflow
  suppressMessages({ # Suppress singular fit warnings
    fit_med_sim <- lmer(sim_med ~ attention + activity + clarity_all_cw2 + interesting_all_cw2 + (1 | PID),
                        data = sim_data)

    fit_out_sim <- lmer(sim_out ~ attention + sim_med + activity + clarity_all_cw2 + interesting_all_cw2 + (1 | PID),
                        data = sim_data)
  })

  # 6. Run mediation test
  # Using sims=100 for speed. Increase to 1000 for final, more stable results.
  med_results <- suppressMessages(
    mediation::mediate(fit_med_sim, fit_out_sim,
                       treat = 'attention', mediator = 'sim_med',
                       sims = 100)
  )

  # 7. Return the p-value for the indirect effect (ACME)
  return(data.frame(
    N = N,
    p.value = summary(med_results)$ACME.avg.p
  ))
}

# --- 3. Run the Simulation (Step 3) ---

# Define your simulation parameters
N_to_test <- c(100, 150, 200, 250, 300, 350, 400) # Sample sizes to check
n_sims <- 500 # Number of iterations per sample size
              # START WITH 100-500. Run 1000+ for your final paper.

# This may take several minutes!
set.seed(123) # For reproducible results
simulation_results <- purrr::map_dfr(N_to_test, ~{
  cat("Simulating N =", .x, "\n")
  purrr::rerun(n_sims, run_sim_iteration(N = .x, j = j_prompts, esm_inner_orig = esm_inner)) %>%
    purrr::list_rbind()
})

# --- 4. Calculate and Plot Power Curve (Step 4) ---

# Calculate power (proportion of p < .05) for each sample size
power_curve <- simulation_results %>%
  group_by(N) %>%
  summarise(
    power = mean(p.value < 0.05, na.rm = TRUE)
  )

print(power_curve)

# Plot the power curve
ggplot(power_curve, aes(x = N, y = power)) +
  geom_line(linewidth = 1, color = "blue") +
  geom_point(size = 3, color = "blue") +
  geom_hline(yintercept = 0.80, linetype = "dashed", color = "red") +
  scale_y_continuous(labels = scales::percent, limits = c(0, 1)) +
  labs(
    title = "Power Analysis for Mediation Effect",
    subtitle = "Replicating Gross et al. (2025) Analysis",
    x = "Number of Participants (N)",
    y = "Power to Detect Indirect Effect (p < .05)",
    caption = paste(n_sims, "simulations per point")
  ) +
  theme_minimal()
