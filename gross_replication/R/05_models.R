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
## HOW ATTENTION ENTERS -- read this before reading any coefficient.
##
## Both hypotheses are within-person claims: the SAME person feels better in a
## Present moment than in a Mind-Wandering one. Three codings are possible and
## they do not estimate the same thing.
##
##   att_cw2 + att_cb2   the split used here. att_cw2 is the within-person
##                       Present-minus-MW difference, att_cb2 the between-person
##                       association. att_num is 0/1, so att_cw2 is the FULL
##                       difference -- nothing to double.
##
##   attention (raw)     what Gross's Table 3 models do, and what this script
##                       used to do. A random intercept does NOT absorb the
##                       between-person part of a level-1 predictor, so the
##                       coefficient is a blend of the two above. Fit below as
##                       m_a_blend / m_b_blend so the size of the blending is
##                       visible rather than assumed away.
##
##   wb_cw2 outcome      Gross's mediation chain. Person-mean centring the
##                       OUTCOME does isolate the within effect, but crudely: it
##                       forces every person's mean to exactly zero, which
##                       leaves (1 | PID) with no variance to estimate (hence
##                       the singular fits) and, where a participant has a single
##                       prompt, contributes a zero outcome against a varying
##                       predictor. Kept in 07 as a sensitivity variant.
##
## See docs/deviations.md 13. This is a departure from both Gross et al. and the
## registered plan, and is reported as one.
## -----------------------------------------------------------------------------

if (!exists("read_expiwell")) source("R/00_setup.R")
d <- readRDS(file.path(dir_derived, "04_analysis.rds"))
cat("\n", nrow(d), "prompts from", nlevels(d$PID), "participants\n")


## == H1 =======================================================================

## Model A -- attention plus the covariates Gross used in Table 3.
m_a <- lmer(wb ~ att_cw2 + att_cb2 + activity + clarity_cw2 + interesting_cw2 +
              (1 | PID), data = d, REML = use_reml)

cat("\n======== Model A: attention only ======================================\n")
print(summary(m_a))

## Model B -- the same, plus thought valence.
m_b <- lmer(wb ~ att_cw2 + att_cb2 + activity + clarity_cw2 + interesting_cw2 +
              valence_cw2 + (1 | PID), data = d, REML = use_reml)

cat("\n======== Model B: attention, given valence ============================\n")
print(summary(m_b))

cat("\n-- the attention effect ----------------------------------------------\n")
cat("att_cw2 = within-person Present - MW (this is H1)\n")
cat("att_cb2 = between-person: do people who are Present more often report",
    "higher mood overall\n\n")
for (nm in c("Model A", "Model B")) {
  m <- if (nm == "Model A") m_a else m_b
  for (term in c("att_cw2", "att_cb2")) {
    co <- summary(m)$coefficients[term, ]
    ci <- confint(m, term, method = "Wald")
    cat(sprintf("%-8s %-8s b = %+.3f (SE %.3f), 95%% CI [%+.3f, %+.3f], p = %.4f\n",
                nm, term, co["Estimate"], co["Std. Error"], ci[1], ci[2],
                co["Pr(>|t|)"]))
  }
}

## Type III F tests, for the activity omnibus among others.
cat("\n-- Model A, type III tests -------------------------------------------\n")
print(anova(m_a, type = 3))


## -- What the old coding gave -------------------------------------------------
## Attention entered raw, contr.sum coded, exactly as in Gross's Table 3. Fit
## only so the blend can be compared with the split above; not used for H1.

m_a_blend <- lmer(wb ~ attention + activity + clarity_cw2 + interesting_cw2 +
                    (1 | PID), data = d, REML = use_reml)
m_b_blend <- lmer(wb ~ attention + activity + clarity_cw2 + interesting_cw2 +
                    valence_cw2 + (1 | PID), data = d, REML = use_reml)

cat("\n-- blended coding, for comparison only -------------------------------\n")
cat("contr.sum puts +1 on Present, so attention1 is HALF the difference.\n\n")
blend <- bind_rows(lapply(
  list("Model A" = m_a_blend, "Model B" = m_b_blend),
  function(m) {
    co <- summary(m)$coefficients["attention1", ]
    tibble(present_minus_mw = 2 * co["Estimate"], se = 2 * co["Std. Error"],
           p = co["Pr(>|t|)"])
  }), .id = "model")
blend$within_from_split <- c(summary(m_a)$coefficients["att_cw2", "Estimate"],
                             summary(m_b)$coefficients["att_cw2", "Estimate"])
print(as.data.frame(blend), digits = 3, row.names = FALSE)
cat("\nThe gap between the two columns is the between-person contamination.\n")


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
cat("Both att_cw2 and att_cb2 are dropped, so the level-2 column now moves too.\n")
print(as.data.frame(pseudo_r2), digits = 3, row.names = FALSE)


## == H2: mediation by thought valence =========================================
## The three-model chain, on RAW mood and RAW valence with the attention split
## carrying the within-person contrast. Fit with lme4 rather than lmerTest.
##
## Treatment is att_cw2, a continuous variable, going from 0 to 1: within a
## person, holding their overall rate of being Present fixed, that is exactly
## the move from a Mind-Wandering prompt to a Present one.
##
## One compromise. mediate() needs the mediator to be the same variable in both
## models, and the a-path's outcome cannot be valence_cw2 without reintroducing
## the centred-outcome problem. So the mediator is raw `valence`, and the b path
## is a within/between blend. The b path is printed both ways below; in the
## pilot the two are within 0.01 of each other, but check it again on the full
## data rather than trusting that.

fit_total    <- lme4::lmer(wb      ~ att_cw2 + att_cb2 + activity + clarity_cw2 +
                             interesting_cw2 + (1 | PID), data = d, REML = use_reml)
fit_mediator <- lme4::lmer(valence ~ att_cw2 + att_cb2 + activity + clarity_cw2 +
                             interesting_cw2 + (1 | PID), data = d, REML = use_reml)
fit_dv       <- lme4::lmer(wb      ~ att_cw2 + att_cb2 + valence + activity +
                             clarity_cw2 + interesting_cw2 + (1 | PID),
                           data = d, REML = use_reml)

cat("\n-- total effect: attention -> mood -----------------------------------\n")
print(round(summary(fit_total)$coefficients, 4))
cat("\n-- a path: attention -> valence --------------------------------------\n")
print(round(summary(fit_mediator)$coefficients, 4))
cat("\n-- b path: valence -> mood, attention still in ------------------------\n")
print(round(summary(fit_dv)$coefficients, 4))

## The same b path with valence split within/between, to size the compromise.
fit_dv_split <- lme4::lmer(wb ~ att_cw2 + att_cb2 + valence_cw2 + valence_cb2 +
                             activity + clarity_cw2 + interesting_cw2 + (1 | PID),
                           data = d, REML = use_reml)
cat("\n   b path, valence blended  :",
    sprintf("%+.3f", summary(fit_dv)$coefficients["valence", "Estimate"]), "\n")
cat("   b path, valence within   :",
    sprintf("%+.3f", summary(fit_dv_split)$coefficients["valence_cw2", "Estimate"]), "\n")

sing <- c("total", "mediator", "outcome")[c(isSingular(fit_total),
  isSingular(fit_mediator), isSingular(fit_dv))]
cat("\nSingular fits:", if (length(sing)) paste(sing, collapse = ", ") else "none",
    "\n")
cat("(Under Gross's centred outcome all three are singular by construction --",
    "see 07.)\n")

set.seed(seed)
cat("\nRunning mediation with", mediation_sims, "simulations...\n")
med <- mediation::mediate(fit_mediator, fit_dv,
                          treat = "att_cw2", mediator = "valence",
                          control.value = 0, treat.value = 1,
                          sims = mediation_sims)

cat("\n======== H2: mediation ================================================\n")
print(summary(med))

saveRDS(list(m_a = m_a, m_b = m_b, m_a_blend = m_a_blend, m_b_blend = m_b_blend,
             blend = blend, pseudo_r2 = pseudo_r2, med = med,
             fit_total = fit_total, fit_mediator = fit_mediator, fit_dv = fit_dv,
             fit_dv_split = fit_dv_split),
        file.path(dir_derived, "05_models.rds"))
cat("\nSaved:", file.path(dir_derived, "05_models.rds"), "\n")
