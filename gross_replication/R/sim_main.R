#' Runs a simulation for Model 1a OR Model 1b
#'
#' @param N Initial number of participants.
#' @param model_type A string: "model1a" or "model1b".
#' @param esm_inner_orig The original `esm_inner` dataframe.
#' @param prompt_counts_orig A vector of prompt counts.
#' @param beta_attention_ratio Ratio to multiply the 'attention' beta by (for sensitivity).
#' @param use_variable_prompts (Boolean) Sample from prompt_counts?
#' @param attrition_rate (Numeric) Proportion of participants to drop.
#'
#' @return A data.frame with N, p.value, and model_type.
run_sim_main <- function(
                         N,
                         model_type = "model1a",
                         esm_inner_orig = esm_inner,
                         prompt_counts_orig = prompt_counts,
                         beta_attention_ratio = 1,
                         use_variable_prompts = TRUE,
                         attrition_rate = 0.0) {

  ## 1. Select parameters based on model_type
  if (model_type == "model1a") {
    b_to_use <- b_main_1a
    sd_intercept <- sd_u0_1a
    sd_residual <- sd_e_1a
    formula_matrix <- ~ activity + clarity_all_cw2 + interesting_all_cw2 + attention
    formula_lmer <- sim_wb ~ activity + clarity_all_cw2 + interesting_all_cw2 + attention + (1 | PID)
  } else if (model_type == "model1b") {
    b_to_use <- b_main_1b
    sd_intercept <- sd_u0_1b
    sd_residual <- sd_e_1b
    formula_matrix <- ~ activity + clarity_all_cw2 + interesting_all_cw2 + attention + valence_all_cw2
    formula_lmer <- sim_wb ~ activity + clarity_all_cw2 + interesting_all_cw2 + attention + valence_all_cw2 + (1 | PID)
  } else {
    stop("model_type must be 'model1a' or 'model1b'")
  }

  ## 2. Apply Attrition
  N_final <- floor(N * (1 - attrition_rate))
  if (N_final <= 0) return(data.frame(N = N, p.value = NA, model_type = model_type))

  ## 3. Create Bootstrap Sample (Hierarchical Resampling)
  # Sample N_final participant IDs *with replacement*
  sampled_pids <- sample(unique(esm_inner_orig$PID), size = N_final, replace = TRUE)

  # Build the new dataset by grabbing all rows for each sampled PID
  # We create new PIDs (1 to N_final) and map them to the original data
  sim_data_list <- purrr::map(1:N_final, ~{
    orig_pid <- sampled_pids[.x]
    esm_inner_orig %>%
      filter(PID == orig_pid) %>%
      mutate(PID = .x) # Assign new PID
  })

  sim_data <- bind_rows(sim_data_list)
  n_total_obs <- nrow(sim_data)

  ## 3. Define prompt counts
  if (use_variable_prompts) {
    prompts_per_pid <- sample(prompt_counts_orig, size = N_final, replace = TRUE)
  } else {
    prompts_per_pid <- rep(13, N_final) # Use fixed 13
  }

  sim_data <- data.frame(PID = rep(1:N_final, times = prompts_per_pid))
  n_total_obs <- nrow(sim_data)

  ## 4. Simulate Random Intercepts
  u_0j <- rnorm(N_final, mean = 0, sd = sd_intercept)
  sim_data$u_0j_mapped <- u_0j[sim_data$PID]

  ## 5. Simulate covariates
  covariate_indices <- sample(1:nrow(esm_inner_orig), size = n_total_obs, replace = TRUE)
  sim_data$attention <- esm_inner_orig$attention[covariate_indices]
  sim_data$activity <- esm_inner_orig$activity[covariate_indices]
  sim_data$clarity_all_cw2 <- esm_inner_orig$clarity_all_cw2[covariate_indices]
  sim_data$interesting_all_cw2 <- esm_inner_orig$interesting_all_cw2[covariate_indices]
  sim_data$valence_all_cw2 <- esm_inner_orig$valence_all_cw2[covariate_indices] # (Resampled even if not used)

  sim_data$activity <- factor(sim_data$activity,
                             levels = c("Cognitive", "Household", "Other",
                                        "Physical", "Restful", "Social"))

  ## 6. Simulate Fixed Effect component
  X_main <- model.matrix(formula_matrix, data = sim_data)

  b_main_vec <- b_to_use[colnames(X_main)]
  b_main_vec["attention1"] = b_to_use["attention1"] * beta_attention_ratio

  mu_fixed <- X_main %*% b_main_vec

  ## 7. Simulate Residual Error
  e_ij <- rnorm(n_total_obs, mean = 0, sd = sd_residual)

  ## 8. Create the simulated Outcome
  sim_data$sim_wb <- mu_fixed + sim_data$u_0j_mapped + e_ij

  ## 9. Run the analysis
  sim_fit <- tryCatch({
    lmerTest::lmer(formula_lmer, data = sim_data, REML = FALSE)
  }, error = function(e) { NULL })

  ## 10. Extract the p-value for 'attention'
  p.value <- NA
  if (!is.null(sim_fit)) {
    coeffs <- summary(sim_fit)$coefficients
    if ("attention1" %in% rownames(coeffs)) {
      p.value <- coeffs["attention1", "Pr(>|t|)"]
    }
  }

  return(data.frame(N = N, p.value = p.value, model_type = model_type))
}
