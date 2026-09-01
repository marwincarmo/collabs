## -----------------------------------------------------------------------------
## 09  Sanity check: make up data with known effects and see if the model
##     gets them back.
##
## This is not part of the analysis. It is here because the centring, the
## contrast coding and the model formula have to agree with each other, and the
## only way to be sure of that before the real data arrive is to run the same
## steps on data where the right answer is known.
##
## Run it on its own:  Rscript R/09_check_recovery.R
## -----------------------------------------------------------------------------

source("R/00_setup.R")
set.seed(seed)

## -- Make up a complete study at the preregistered target N ------------------
n_people  <- 300
n_prompts <- 13          # Gross's participants averaged 13.2 usable prompts

true_b <- c(attention = 0.25, clarity = 0.10, interesting = 0.08, valence = 0.45)
true_intercept <- 0.5
true_sd_person <- 0.8    # between-person SD of mood
true_sd_error  <- 1.0

sim <- tibble(
  PID         = rep(sprintf("P%03d", 1:n_people), each = n_prompts),
  attention   = factor(sample(c("Present", "MW"), n_people * n_prompts, TRUE, c(.75, .25)),
                       levels = attention_levels),
  activity    = factor(sample(activity_levels, n_people * n_prompts, TRUE),
                       levels = activity_levels),
  clarity     = sample(1:7,   n_people * n_prompts, TRUE),
  valence     = sample(-3:3,  n_people * n_prompts, TRUE),
  interesting = sample(1:7,   n_people * n_prompts, TRUE)
)

## Centre exactly the way 04 does, so the generating model and the fitted model
## are on the same scale.
sim <- sim |>
  group_by(PID) |>
  mutate(clarity_cw2     = clarity     - mean(clarity),
         valence_cw2     = valence     - mean(valence),
         interesting_cw2 = interesting - mean(interesting)) |>
  ungroup()

## contr.sum with Present first means Present = +1 and MW = -1.
attention_effect <- ifelse(sim$attention == "Present", 1, -1)
person_intercept <- rnorm(n_people, 0, true_sd_person)[match(sim$PID, unique(sim$PID))]

sim$wb <- true_intercept +
  true_b["attention"]   * attention_effect +
  true_b["clarity"]     * sim$clarity_cw2 +
  true_b["interesting"] * sim$interesting_cw2 +
  true_b["valence"]     * sim$valence_cw2 +
  person_intercept + rnorm(nrow(sim), 0, true_sd_error)

contrasts(sim$attention) <- contr.sum(2)

## -- Fit the same model 05 fits ----------------------------------------------
m <- lmer(wb ~ attention + activity + clarity_cw2 + interesting_cw2 + valence_cw2 +
            (1 | PID), data = sim, REML = use_reml)
co <- summary(m)$coefficients

comparison <- tibble(
  parameter = c("attention (half the Present-MW gap)", "clarity_cw2",
                "interesting_cw2", "valence_cw2", "person SD", "residual SD"),
  true      = c(true_b[["attention"]], true_b[["clarity"]], true_b[["interesting"]],
                true_b[["valence"]], true_sd_person, true_sd_error),
  estimated = c(co["attention1", "Estimate"], co["clarity_cw2", "Estimate"],
                co["interesting_cw2", "Estimate"], co["valence_cw2", "Estimate"],
                as.data.frame(VarCorr(m))$sdcor[1], sigma(m)),
  se        = c(co["attention1", "Std. Error"], co["clarity_cw2", "Std. Error"],
                co["interesting_cw2", "Std. Error"], co["valence_cw2", "Std. Error"],
                NA, NA)
) |>
  mutate(difference = estimated - true,
         within_2se = ifelse(is.na(se), NA, abs(difference) < 2 * se))

cat("\n---- Recovering known values from made-up data -----------------------\n")
cat(n_people, "participants x", n_prompts, "prompts\n\n")
print(as.data.frame(comparison), digits = 3, row.names = FALSE)

cat("\nSign check: attention1 is",
    ifelse(co["attention1", "Estimate"] > 0, "POSITIVE, as it should be when mood",
           "NEGATIVE -- something is wrong with the contrast coding"),
    "\nis higher while present.\n")

cat("\nReported as the full difference: Present - MW =",
    sprintf("%+.3f (true %+.3f)\n", 2 * co["attention1", "Estimate"],
            2 * true_b[["attention"]]))
