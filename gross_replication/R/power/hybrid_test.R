## Step 1: Fit Bayesian model to original data once
library(brms)
library(progressr)
source("R/power/hybrid_sim.R")

## Fit model with weakly informative priors
bayes_model1a <- brm(
  wb ~ attention + clarity_all_cw2 + interesting_all_cw2 + activity + (1|PID),
  data = esm_inner,
  prior = c(
    prior(normal(0, 1), class = "b"),
    prior(normal(0, 1), class = "sd"),
    prior(normal(0, 2), class = "sigma")
  ),
  cores = 4,
  chains = 4,
  iter = 2000,
  backend = "cmdstanr"  # Faster sampling
)

bayes_model1b <- brm(
  wb ~ attention + clarity_all_cw2 + interesting_all_cw2 + valence_all_cw2 + activity + (1|PID),
  data = esm_inner,
  prior = c(
    prior(normal(0, 1), class = "b"),
    prior(normal(0, 1), class = "sd"),
    prior(normal(0, 2), class = "sigma")
  ),
  cores = 4,
  chains = 4,
  iter = 2000,
  backend = "cmdstanr"  # Faster sampling
)

## Extract posterior draws
posteriorA <- as.data.frame(bayes_model1a)
posteriorB <- as.data.frame(bayes_model1b)


# 3. DEFINE data generation (YOUR LOGIC, REFACTORED)
single_draw <- posteriorB[sample(1:nrow(posteriorB), 1), ]
param_draw = extract_params_from_draw(single_draw)
N <- 100

## 5. RUN power analysis
plan(multisession, workers = 12)

param_grid <- expand.grid(
  sample_size = c(100, 150, 200, 250, 300, 350),
  attenuation_factor = c(1, 0.7),
  completion_adjust = c(1, 0.7),
  use_variable_prompts = TRUE
) |>
  dplyr::mutate(scenario_id = dplyr::row_number())

nsims <- 1000
total_sims <- nrow(param_grid) * nsims

power_results <- with_progress({
  p <- progressr::progressor(total_sims)
  furrr::future_map_dfr(
           1:nrow(param_grid),
           function(row_idx) {
             params <- param_grid[row_idx, ]
             pvals <- furrr::future_map(
                               1:nsims,
                               function(x) {
                                        # Run simulation
                                 result <- run_single_bayes_sim(N = params$sample_size,
                                                                posterior_draws = posterior_draws,
                                                                use_variable_prompts = params$use_variable_prompts,
                                                                attenuation_factor = params$attenuation_factor,
                                                                completion_adjust = params$completion_adjust)
                                 p()
                                 return(result)
                               },
                               .options = furrr_options(seed = TRUE)
                             )
             data.frame(scenario_id = params$scenario_id,
                        sample_size = params$sample_size,
                        power = mean(unlist(pvals) < 0.05, na.rm = TRUE),
                        n_converged = sum(!is.na(pvals)),
                        use_variable_prompts = params$use_variable_prompts,
                        attenuation_factor = params$attenuation_factor,
                        completion_adjust = params$completion_adjust,
                        total_simulations = nsims
                        )
           },
           .options = furrr_options(seed = TRUE)
         )
})


saveRDS(power_results, "output/power_analysis_hybrid.rds")

power_results$scenario <- ifelse(power_results$attenuation_factor == 1 & power_results$completion_adjust == 1, "Ideal",
                          ifelse(power_results$attenuation_factor == 0.7 & power_results$completion_adjust == 1, "70% effect",
                          ifelse(power_results$attenuation_factor == 1 & power_results$completion_adjust == 0.7, "70% completion",
                          ifelse(power_results$attenuation_factor == 0.7 & power_results$completion_adjust == 0.7, "Challenging", NA)
                          )))

plan(sequential)

# Plot
ggplot(power_results, aes(x = sample_size, y = power, color = scenario, group = scenario)) +
  geom_line(linewidth = 1) +
  geom_point(size = 2.5) +
  geom_hline(yintercept = 0.80, linetype = "dashed", color = "red") +
  scale_y_continuous(labels = scales::percent, limits = c(0, 1)) +
  labs(
    title = "Power Analysis for Main Attention Effects",
    x = "Number of Participants Recruited (N)",
    y = "Power (p < .05)",
    color = "Analysis Scenario"
  ) +
  theme_minimal()

run_single_bayes_sim(N = 100, posterior_draws = posteriorB, model_type = "B")
