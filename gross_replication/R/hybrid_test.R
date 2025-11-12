## Step 1: Fit Bayesian model to original data once
library(brms)

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

## Extract posterior draws
posterior_draws <- as.data.frame(bayes_model1a)


# 2. DEFINE parameter extraction function
extract_params_from_draw <- function(draw) {
  list(
    b_Intercept = draw$b_Intercept,
    b_attention = draw$b_attention1,
    b_clarity = draw$b_clarity_all_cw2,
    b_interesting = draw$b_interesting_all_cw2,
    b_activities = draw[, grep("^b_activity", colnames(draw))],
    sd_u0 = draw$sd_PID__Intercept,
    sd_e = draw$sigma
  )
}

# 3. DEFINE data generation (YOUR LOGIC, REFACTORED)
single_draw <- posterior_draws[sample(1:nrow(posterior_draws), 1), ]
param_draw = extract_params_from_draw(single_draw)
N <- 100

#' Generate simulated ESM data with preserved correlation structure
#'
#' @param N Number of participants
#' @param param_draw Single posterior draw or parameter list
#' @param prompt_counts_orig Vector of empirical prompt counts
#' @param esm_inner_orig Original data frame
#'
#' @return Simulated data frame
generate_data_joint <- function(N,
                                param_draw,
                                prompt_counts_orig = prompt_counts,
                                esm_inner_orig = esm_inner,
                                use_variable_prompts = TRUE,
                                attenuation_factor = 1,
                                completion_average = 13,  ## average prompt completion in the original research
                                completion_adjust = 0.8) {

  ## 1. Sample participants WITH their full covariance structure
  n_original_participants <- length(unique(esm_inner_orig$PID))

  # Create a "template participant" for each simulated person
  # This preserves the attention-activity-covariate correlations
  sampled_pids_original <- sample(unique(esm_inner_orig$PID), size = N, replace = TRUE)

  ## 2. Build simulated dataset by copying and perturbing participant histories
  sim_data_list <- lapply(seq_along(sampled_pids_original), function(i) {

    # Get all prompts from this original participant
    pid_data <- esm_inner_orig |> dplyr::filter(PID == sampled_pids_original[i])
    n_available_prompts <- nrow(pid_data)

    if (use_variable_prompts) {
      # Sample a subset of prompts to match desired prompt count
      n_desired_prompts <- sample(prompt_counts_orig, size = 1, replace = TRUE)
    } else {
      n_desired_prompts <- min(rpois(1, completion_average), 42)
    }

    ## adjust by desired completion rate
    ## ceiling() ensure values > 0
    n_desired_prompts <- ceiling(n_desired_prompts * completion_adjust)

    # If we need fewer prompts than available, sample randomly
    # If we need more, sample with replacement (duplicate prompts)
    pid_data <- pid_data |> dplyr::select(,
      attention, activity, clarity_all_cw2, interesting_all_cw2, valence_all_cw2
      )

    if (n_desired_prompts == n_available_prompts) {
      # Perfect match - use original data
      person_sim <- pid_data

    } else if (n_desired_prompts < n_available_prompts) {
      # Sample subset (no duplication needed)
      person_sim <- pid_data |> dplyr::slice_sample(n = n_desired_prompts, replace = FALSE)

    } else {
      # CASE C: Need more prompts than available
      # Sample with replacement, then add synthetic variation

      # Base sample (with duplicates)
      person_sim <- pid_data |> dplyr::slice_sample(n = n_desired_prompts, replace = TRUE)

      # ADD SYNTHETIC WITHIN-PERSON VARIATION
      # Extract the duplicated values
      duplicated_indices <- which(duplicated(person_sim |> dplyr::select(attention, activity)))

      if (length(duplicated_indices) > 0) {
        # For duplicated rows, add person-specific noise
        # This preserves the mean but adds realistic variation

        # Person-level random effect for each covariate (constant for this person)
        person_sim$clarity_all_cw2[duplicated_indices] <-
          person_sim$clarity_all_cw2[duplicated_indices] + rnorm(length(duplicated_indices), 0, .1)

        person_sim$interesting_all_cw2[duplicated_indices] <-
          person_sim$interesting_all_cw2[duplicated_indices] + rnorm(length(duplicated_indices), 0, .1)

        person_sim$valence_all_cw2[duplicated_indices] <-
          person_sim$valence_all_cw2[duplicated_indices] + rnorm(length(duplicated_indices), 0, .1)
      }
    }

    # Assign new participant ID
    person_sim$PID <- i

    return(person_sim)
  })

  ## 3. Combine participants
  sim_data <- dplyr::bind_rows(sim_data_list)

  ## 4. Add person-specific random intercept
  u_0j <- rnorm(N, mean = 0, sd = param_draw$sd_u0)
  sim_data$u_0j_mapped <- u_0j[sim_data$PID]

  ## 5. Generate outcome with residual error
  n_total_obs <- nrow(sim_data)

  # Build model matrix (ensures correct contrasts)
  X <- model.matrix(~ attention + activity + clarity_all_cw2 + interesting_all_cw2,
                    data = sim_data)


  b_vector <- unlist(param_draw[-grep("^sd", names(param_draw))]) * attenuation_factor
  mu_fixed <- X %*% b_vector

  sim_data$sim_wb <- as.numeric(mu_fixed) +
                     sim_data$u_0j_mapped +
                     rnorm(n_total_obs, sd = param_draw$sd_e)

  return(sim_data)
}

# 4. DEFINE simulation wrapper
run_single_bayes_sim <- function(N, posterior_draws, ...) {

  # Sample one posterior draw
  draw <- posterior_draws[sample(1:nrow(posterior_draws), 1), ]
  params <- extract_params_from_draw(draw)

  # Generate data
  sim_data <- generate_data_joint(N, params, prompt_counts, ...)
  contrasts(sim_data$attention) <- contr.sum(2)

  # Fit frequentist model
  model <- lmerTest::lmer(sim_wb ~ activity + clarity_all_cw2 + interesting_all_cw2 + attention + (1 | PID),
                          data = sim_data, REML = FALSE)

  # Return p-value
  summary(model)$coefficients["attention1", "Pr(>|t|)"]
}

sample_size <- c(100, 125, 150, 200, 225, 250, 300, 350)
## 5. RUN power analysis
plan(multisession, workers = availableCores() - 1)

power_results <- furrr::future_map_dfr(sample_size, ~{
  pvals <- furrr::future_map(1:1000,
                             function(x) run_single_bayes_sim(
                                           N = .x,
                                           posterior_draws = posterior_draws),
                             .options = furrr_options(seed = TRUE) )
  data.frame(N = .x, power = mean(unlist(pvals) < 0.05))
},
.options = furrr_options(seed = TRUE) )

plan(sequential)

ggplot(power_results, aes(x = N, y = power)) +
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
