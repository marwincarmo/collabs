## -----------------------------------------------------------------------------
## 09  Sanity check: make up data with known effects and see if the models get
##     them back.
##
## Not part of the analysis. It is here because the centring, the coding and the
## model formula have to agree with each other, and the only way to be sure of
## that before the real data arrive is to run the same steps on data where the
## right answer is known.
##
## Two things are checked.
##
##  A. The primary specification (att_cw2 + att_cb2) recovers the true
##     within-person effect when a between-person confound is present, and the
##     two codings this analysis rejects do not. This is the simulation evidence
##     behind docs/deviations.md 13.
##  B. What sparse data does to each coding -- half the participants with a
##     single prompt, as in the current pilot.
##
## Run it on its own:  Rscript R/09_check_recovery.R
## -----------------------------------------------------------------------------

source("R/00_setup.R")
set.seed(seed)

n_people  <- 300
n_prompts <- 13          # Gross's participants averaged 13.2 usable prompts

true_within    <- 0.50   # the full Present - MW difference WITHIN a person
true_b         <- c(clarity = 0.10, interesting = 0.08, valence = 0.45)
true_intercept <- 0.5
true_sd_person <- 0.8    # between-person SD of mood
true_sd_error  <- 1.0
confound       <- 1.2    # how strongly being Present-prone tracks a high intercept
n_reps         <- 25     # a single draw is too noisy to read a bias off


## -- A study where attention is confounded between people ---------------------
## Each person gets an intercept u_i and a propensity to be Present that rises
## with it. So people who are Present more often really do have higher mood, on
## top of the within-person effect. That is the situation the split exists for:
## a random intercept alone will NOT separate the two.

make_study <- function(prompts_per_person) {
  u <- rnorm(n_people, 0, true_sd_person)
  p <- plogis(0.8 + confound * u / true_sd_person * 0.6)

  sim <- tibble(
    PID     = rep(sprintf("P%03d", 1:n_people), times = prompts_per_person),
    u       = rep(u, times = prompts_per_person),
    att_num = rbinom(sum(prompts_per_person), 1,
                     rep(p, times = prompts_per_person))
  ) |>
    mutate(
      attention   = factor(ifelse(att_num == 1, "Present", "MW"),
                           levels = attention_levels),
      activity    = factor(sample(activity_levels, n(), TRUE), levels = activity_levels),
      clarity     = sample(1:7,  n(), TRUE),
      valence     = sample(-3:3, n(), TRUE),
      interesting = sample(1:7,  n(), TRUE)
    )

  ## Centre exactly the way 04 does, so the generating and fitted models agree.
  sim <- sim |>
    group_by(PID) |>
    mutate(clarity_cw2     = clarity     - mean(clarity),
           valence_cw2     = valence     - mean(valence),
           interesting_cw2 = interesting - mean(interesting),
           att_mean2       = mean(att_num)) |>
    ungroup() |>
    mutate(att_cw2 = att_num - att_mean2,
           att_c2  = att_num - mean(att_num),
           att_cb2 = att_c2 - att_cw2)

  ## att_num enters at its true within-person weight; the confound arrives only
  ## through u, which is shared by the person's prompts.
  sim$wb <- true_intercept +
    true_within           * sim$att_num +
    true_b["clarity"]     * sim$clarity_cw2 +
    true_b["interesting"] * sim$interesting_cw2 +
    true_b["valence"]     * sim$valence_cw2 +
    sim$u + rnorm(nrow(sim), 0, true_sd_error)

  sim <- sim |> group_by(PID) |> mutate(wb_cw2 = wb - mean(wb)) |> ungroup()
  contrasts(sim$attention) <- contr.sum(2)
  sim
}

covars <- "+ activity + clarity_cw2 + interesting_cw2 + valence_cw2 + (1 | PID)"

## The three codings, each reported as the full Present - MW difference.
codings <- list(
  "Split (att_cw2 + att_cb2)" = list(f = paste("wb ~ att_cw2 + att_cb2", covars),
                                     term = "att_cw2",    k = 1),
  "Blended (attention raw)"   = list(f = paste("wb ~ attention", covars),
                                     term = "attention1", k = 2),
  "Centred outcome (wb_cw2)"  = list(f = paste("wb_cw2 ~ attention", covars),
                                     term = "attention1", k = 2)
)

fit_one <- function(sim, x, reml) {
  m  <- lmer(as.formula(x$f), data = sim, REML = reml)
  co <- summary(m)$coefficients[x$term, ]
  ## unname: co["Estimate"] carries its own name, which would corrupt the
  ## row labels the caller indexes on.
  c(est  = unname(x$k * co["Estimate"]),
    se   = unname(x$k * co["Std. Error"]),
    sing = as.numeric(isSingular(m)))
}

## One draw is nowhere near enough to read a bias off -- the Monte Carlo error on
## a single fit is about the size of the bias being looked for. Averaged instead.
three_ways <- function(prompts_per_person, reps = n_reps, reml = use_reml) {
  draws <- replicate(reps, {
    sim <- make_study(prompts_per_person)
    vapply(codings, fit_one, numeric(3), sim = sim, reml = reml)
  }, simplify = "array")

  tibble(
    coding           = names(codings),
    ## draws is [statistic, coding, replication]; average over replications.
    present_minus_mw = apply(draws["est", , ], 1, mean),
    bias             = apply(draws["est", , ], 1, mean) - true_within,
    mean_se          = apply(draws["se", , ], 1, mean),
    singular         = paste0(round(100 * apply(draws["sing", , ], 1, mean)), "%")
  )
}


## == A. Complete data, 13 prompts each ========================================

cat("\n---- A. Recovering a known within-person effect ----------------------\n")
cat(n_people, "participants x", n_prompts, "prompts, with a between-person",
    "confound built in;\n", n_reps, "replications.\n")
cat("True within-person Present - MW =", sprintf("%+.3f", true_within), "\n\n")
print(as.data.frame(three_ways(rep(n_prompts, n_people))), digits = 3,
      row.names = FALSE)
cat("\nOnly the split recovers the truth. The blended coding runs HIGH, having\n",
    "absorbed the between-person confound. The centred outcome runs LOW, and is\n",
    "singular in every replication -- by construction, not by bad luck.\n")

sim <- make_study(rep(n_prompts, n_people))

## The rest of Model B, from the primary specification.
m <- lmer(as.formula(paste("wb ~ att_cw2 + att_cb2", covars)), data = sim,
          REML = use_reml)
co <- summary(m)$coefficients

comparison <- tibble(
  parameter = c("att_cw2 (Present - MW, within)", "clarity_cw2",
                "interesting_cw2", "valence_cw2", "person SD", "residual SD"),
  true      = c(true_within, true_b[["clarity"]], true_b[["interesting"]],
                true_b[["valence"]], true_sd_person, true_sd_error),
  estimated = c(co["att_cw2", "Estimate"], co["clarity_cw2", "Estimate"],
                co["interesting_cw2", "Estimate"], co["valence_cw2", "Estimate"],
                as.data.frame(VarCorr(m))$sdcor[1], sigma(m)),
  se        = c(co["att_cw2", "Std. Error"], co["clarity_cw2", "Std. Error"],
                co["interesting_cw2", "Std. Error"], co["valence_cw2", "Std. Error"],
                NA, NA)
) |>
  mutate(difference = estimated - true,
         within_2se = ifelse(is.na(se), NA, abs(difference) < 2 * se))

cat("\n-- the rest of the primary model -------------------------------------\n")
print(as.data.frame(comparison), digits = 3, row.names = FALSE)

cat("\nSign check: att_cw2 is",
    ifelse(co["att_cw2", "Estimate"] > 0,
           "POSITIVE, as it should be when mood is higher while present.",
           "NEGATIVE -- something is wrong with the coding."), "\n")

cat("\nPerson SD comes in BELOW the", true_sd_person, "used to generate it, and",
    "should: att_cb2 =", sprintf("%+.3f", co["att_cb2", "Estimate"]),
    "\nis a person-level term, so it explains part of the between-person",
    "variance that\nwould otherwise land in the random intercept. The two are",
    "not comparable\nagainst the generating value in isolation.\n")


## == B. What sparse data does ==================================================
## Half the participants contribute one prompt, as in the current pilot. The
## within-person effect is unchanged; only the design is thinner.

sparse_counts <- rep(c(1L, n_prompts), each = n_people / 2)

cat("\n---- B. The same study with half the sample at one prompt ------------\n")
cat(sum(sparse_counts), "prompts;", sum(sparse_counts == 1),
    "participants with exactly one;", n_reps, "replications.\n")
cat("True within-person Present - MW =", sprintf("%+.3f", true_within), "\n\n")
print(as.data.frame(three_ways(sparse_counts)), digits = 3, row.names = FALSE)

cat("\nThe split loses precision but stays unbiased: a single-prompt participant\n",
    "has att_cw2 exactly 0 and simply sits out the within slope.\n")
cat("Both rejected codings get WORSE as the data thin. The centred outcome\n",
    "attenuates further because those participants have wb_cw2 exactly 0 while\n",
    "attention still varies, which drags the slope toward zero.\n")
