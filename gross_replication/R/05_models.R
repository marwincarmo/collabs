## -----------------------------------------------------------------------------
## 05  The two preregistered hypotheses.
##
##  H1  Mood is higher when Present than when Mind-Wandering, controlling for
##      current activity, thought clarity and thought interestingness.
##  H2  That relationship is mediated by thought valence.
##
## REML, as the prereg specifies. (R/power/ used ML for simulation speed;
## 07 refits under ML so the difference is visible.)
##
## Reading the attention coefficient: contr.sum codes Present +1 and MW -1, so
## `attention1` is HALF the Present-minus-MW difference. Doubling is done
## explicitly below, never silently.
## -----------------------------------------------------------------------------

if (!exists("read_expiwell")) source("R/00_setup.R")
d <- readRDS(file.path(dir_derived, "04_analysis.rds"))
cat("\n", nrow(d), "prompts from", nlevels(d$PID), "participants\n")


## == H1 =======================================================================

## Model A -- attention plus the covariates Gross used in Table 3.
m_a <- lmer(wb ~ attention + activity + clarity_cw2 + interesting_cw2 + (1 | PID),
            data = d, REML = use_reml)

cat("\n======== Model A: attention only ======================================\n")
print(summary(m_a))

## Model B -- the same, plus thought valence.
m_b <- lmer(wb ~ attention + activity + clarity_cw2 + interesting_cw2 + valence_cw2 +
              (1 | PID), data = d, REML = use_reml)

cat("\n======== Model B: attention, given valence ============================\n")
print(summary(m_b))

cat("\n-- the attention effect ----------------------------------------------\n")
for (nm in c("Model A", "Model B")) {
  m  <- if (nm == "Model A") m_a else m_b
  co <- summary(m)$coefficients["attention1", ]
  ci <- confint(m, "attention1", method = "Wald")
  cat(sprintf("%s  b = %+.3f (SE %.3f), p = %.4f\n", nm, co["Estimate"], co["Std. Error"],
              co["Pr(>|t|)"]))
  cat(sprintf("          Present - MW = %+.3f, 95%% CI [%+.3f, %+.3f]\n",
              2 * co["Estimate"], 2 * ci[1], 2 * ci[2]))
}

## Type III F tests, for the activity omnibus among others.
cat("\n-- Model A, type III tests -------------------------------------------\n")
print(anova(m_a, type = 3))


## == Variance explained by attention state ====================================
## pseudo-R2 in the sense of Blanke et al. (2018): how much of each variance
## component disappears when attention state is added. Both models are fit the
## same way on the same rows, so this is a variance comparison, not a test.

varcomp <- function(m) {
  v <- as.data.frame(VarCorr(m))
  c(level1 = v$vcov[v$grp == "Residual"],
    level2 = v$vcov[v$grp == "PID" & is.na(v$var2)])
}

m_a_noatt <- lmer(wb ~ activity + clarity_cw2 + interesting_cw2 + (1 | PID),
                  data = d, REML = use_reml)
m_b_noatt <- lmer(wb ~ activity + clarity_cw2 + interesting_cw2 + valence_cw2 +
                    (1 | PID), data = d, REML = use_reml)

pseudo_r2 <- tibble(
  model  = c("Model A", "Model B"),
  var_l1_without = c(varcomp(m_a_noatt)["level1"], varcomp(m_b_noatt)["level1"]),
  var_l1_with    = c(varcomp(m_a)["level1"],       varcomp(m_b)["level1"]),
  var_l2_without = c(varcomp(m_a_noatt)["level2"], varcomp(m_b_noatt)["level2"]),
  var_l2_with    = c(varcomp(m_a)["level2"],       varcomp(m_b)["level2"])
) |>
  mutate(pseudo_r2_within  = (var_l1_without - var_l1_with) / var_l1_without,
         pseudo_r2_between = (var_l2_without - var_l2_with) / var_l2_without)

cat("\n-- pseudo-R2 for attention state -------------------------------------\n")
print(as.data.frame(pseudo_r2), digits = 3, row.names = FALSE)


## == H2: mediation by thought valence =========================================
## The three-model chain from Gross's own script, on within-person centred mood.
## Fit with lme4 rather than lmerTest, and with treatment contrasts, because
## mediate() wants an explicit control and treatment value.

fit_total    <- lme4::lmer(wb_cw2      ~ attention_tx + activity + clarity_cw2 +
                             interesting_cw2 + (1 | PID), data = d, REML = use_reml)
fit_mediator <- lme4::lmer(valence_cw2 ~ attention_tx + activity + clarity_cw2 +
                             interesting_cw2 + (1 | PID), data = d, REML = use_reml)
fit_dv       <- lme4::lmer(wb_cw2      ~ attention_tx + valence_cw2 + activity +
                             clarity_cw2 + interesting_cw2 + (1 | PID),
                           data = d, REML = use_reml)

cat("\n-- total effect: attention -> mood -----------------------------------\n")
print(round(summary(fit_total)$coefficients, 4))
cat("\n-- a path: attention -> valence --------------------------------------\n")
print(round(summary(fit_mediator)$coefficients, 4))
cat("\n-- b path: valence -> mood, attention still in ------------------------\n")
print(round(summary(fit_dv)$coefficients, 4))

cat("\nSingular fits (a warning sign at small N):",
    paste(c("total", "mediator", "outcome")[c(isSingular(fit_total),
      isSingular(fit_mediator), isSingular(fit_dv))], collapse = ", "), "\n")

set.seed(seed)
cat("\nRunning mediation with", mediation_sims, "simulations...\n")
med <- mediation::mediate(fit_mediator, fit_dv,
                          treat = "attention_tx", mediator = "valence_cw2",
                          control.value = "MW", treat.value = "Present",
                          sims = mediation_sims)

cat("\n======== H2: mediation ================================================\n")
print(summary(med))

saveRDS(list(m_a = m_a, m_b = m_b, pseudo_r2 = pseudo_r2, med = med,
             fit_total = fit_total, fit_mediator = fit_mediator, fit_dv = fit_dv),
        file.path(dir_derived, "05_models.rds"))
cat("\nSaved:", file.path(dir_derived, "05_models.rds"), "\n")
