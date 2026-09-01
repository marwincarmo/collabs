#' Extract Parameters from a Single Posterior Draw
#'
#' Formats a raw row from a brms posterior dataframe into a named list
#' usable for data generation.
#'
#' @param draw A single row from as_draws_df(brms_model)
#' @return A named list containing fixed effects (betas) and variance components.
#'
extract_params_from_draw <- function(draw) {
  as.list(c(
    ## extract all fixed effects
    draw[grep("^b", names(draw))],
    ## extract random intercept sd
    "sd_u0" = draw$sd_PID__Intercept,
    ## extract residual sd
    "sd_e" = draw$sigma))
}

#' Generate Simulated ESM Data via Hierarchical Resampling
#'
#' This function creates a synthetic dataset by resampling participants from the
#' original data. This preserves the joint distribution and correlation structure
#' of the covariates (e.g., age, clarity, activity) without needing to model them explicitly.
#'
#' @param N Integer. Number of participants to simulate.
#' @param param_draw List. Parameters extracted from `extract_params_from_draw`.
#' @param prompt_counts_orig Vector. Empirical distribution of prompts per person.
#' @param esm_inner_orig Dataframe. The original "Inner Speech" dataset.
#' @param use_variable_prompts Logical. If TRUE, samples prompt counts from real data.
#' @param model_type String. "A" (Attention + covariates) or "B" (Attention + Covariates + Thought valence).
#' @param attenuation_factor Numeric. Scalar to shrink effect sizes.
#' @param completion_average Integer. Mean for Poisson sampling if variable prompts are FALSE.
#' @param completion_adjust Numeric. Multiplier to simulate attrition/lower completion rates.
#'
#' @return A dataframe.
generate_data_joint <- function(N,
                                param_draw,
                                prompt_counts_orig = prompt_counts,
                                esm_inner_orig = esm_inner,
                                use_variable_prompts = TRUE,
                                model_type,
                                attenuation_factor = 1,
                                completion_average = 13,
                                completion_adjust = 1) {

  ## Define predictors based on model choice
  if (model_type == "A") {
    formula_matrix <- ~ attention + activity + clarity_all_cw2 + interesting_all_cw2
  } else if (model_type == "B") {
    formula_matrix <- ~ activity + clarity_all_cw2 + interesting_all_cw2 + attention + valence_all_cw2
  } else {
    stop("model must be 'A' or 'B'")
  }

  ## 1. Sample participants WITH their full covariance structure
  ## Instead of generating covariates from scratch, we sample real participants
  ## This preserves the within-person dependencies and Level 2 covariate structures
  n_original_participants <- length(unique(esm_inner_orig$PID))

  # Create a "template participant" for each simulated person
  sampled_pids_original <- sample(unique(esm_inner_orig$PID), size = N, replace = TRUE)

  ## 2. Build simulated dataset by copying and perturbing participant histories
  sim_data_list <- lapply(seq_along(sampled_pids_original), function(i) {

    ## Get all prompts from this original participant
    pid_data <- esm_inner_orig |> dplyr::filter(PID == sampled_pids_original[i])
    n_available_prompts <- nrow(pid_data)

    if (use_variable_prompts) {
      ## Sample a subset of prompts to match desired prompt count
      n_desired_prompts <- sample(prompt_counts_orig, size = 1, replace = TRUE)
    } else {
      ## Theoretical approach: Poisson distribution capped at 42 (max theoretical prompts)
      n_desired_prompts <- min(rpois(1, completion_average), 42)
    }

    ## Adjust by desired completion rate
    ## ceiling() ensure values > 0
    n_desired_prompts <- ceiling(n_desired_prompts * completion_adjust)

    ## If we need fewer prompts than available, sample randomly
    ## If we need more, sample with replacement (duplicate prompts)
    pid_data <- pid_data[, all.vars(formula_matrix)]

    if (n_desired_prompts == n_available_prompts) {
      ## Perfect match - use original data
      person_sim <- pid_data

    } else if (n_desired_prompts < n_available_prompts) {
      ## Sample subset (no duplication needed)
      person_sim <- pid_data |> dplyr::slice_sample(n = n_desired_prompts, replace = FALSE)

    } else {
      ## CASE C: Need more prompts than available
      ## Sample with replacement, then add synthetic variation

      ## Base sample (with duplicates)
      person_sim <- pid_data |> dplyr::slice_sample(n = n_desired_prompts, replace = TRUE)

      ## ADD SYNTHETIC WITHIN-PERSON VARIATION
      ## Extract the duplicated values
      duplicated_indices <- which(duplicated(person_sim |> dplyr::select(attention, activity)))

      if (length(duplicated_indices) > 0) {
        ## For duplicated rows, add person-specific noise
        ## This preserves the mean but adds realistic variation

        ## Person-level random effect for each covariate (constant for this person)
        person_sim$clarity_all_cw2[duplicated_indices] <-
          person_sim$clarity_all_cw2[duplicated_indices] + rnorm(length(duplicated_indices), 0, .1)

        person_sim$interesting_all_cw2[duplicated_indices] <-
          person_sim$interesting_all_cw2[duplicated_indices] + rnorm(length(duplicated_indices), 0, .1)

        person_sim$valence_all_cw2[duplicated_indices] <-
          person_sim$valence_all_cw2[duplicated_indices] + rnorm(length(duplicated_indices), 0, .1)
      }
    }

    ## Assign new participant ID
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

  ## Build model matrix
  X <- model.matrix(formula_matrix,
                    data = sim_data)

  ## Extract betas
  b_vector <- unlist(param_draw[-grep("^sd", names(param_draw))]) * attenuation_factor
  ## Remove the "b_" prefix
  names(b_vector) <- sub("^b_", "", names(b_vector))
  ## If necessary, fix naming mismatches
  names(b_vector)[names(b_vector) == "attention1"] <- "attentionPresent"
  names(b_vector)[names(b_vector) == "Intercept"] <- "(Intercept)"
  ## Sort coefficients to match the order of columns in X
  b_vector <- b_vector[colnames(X)]
  ## Calculate fixed effects
  mu_fixed <- X %*% b_vector

  sim_data$sim_wb <- as.numeric(mu_fixed) +
                     sim_data$u_0j_mapped +
                     rnorm(n_total_obs, sd = param_draw$sd_e)

  return(sim_data)
}

## 4. Defome simulation wrapper
run_single_bayes_sim <- function(N, posterior_draws, model_type, prompt_counts, ...) {

  ## select predictors based on model
  if (model_type == "A") {
    formula_lmer <- sim_wb ~ activity + clarity_all_cw2 + interesting_all_cw2 + attention + (1 | PID)
  } else if (model_type == "B") {
    formula_lmer <- sim_wb ~ activity + clarity_all_cw2 + interesting_all_cw2 + attention + valence_all_cw2 + (1 | PID)
  } else {
    stop("model must be 'A' or 'B'")
  }

  ## Sample one posterior draw (the "true" parameters)
  draw <- posterior_draws[sample(1:nrow(posterior_draws), 1), ]
  params <- extract_params_from_draw(draw)

  ## Generate data
  sim_data <- generate_data_joint(N, params, prompt_counts, model_type = model_type, ...)
  contrasts(sim_data$attention) <- contr.sum(2)

  # Fit frequentist model
  model <- lmerTest::lmer(formula_lmer,
                          data = sim_data, REML = FALSE)

  ## Return p-value
  summary(model)$coefficients["attention1", "Pr(>|t|)"]
}
