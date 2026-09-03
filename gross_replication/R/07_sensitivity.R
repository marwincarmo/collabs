## -----------------------------------------------------------------------------
## 07  Does the H1 result depend on the choices we had to make?
##
## Two kinds of choice are varied here, and they are independent of each other.
##
## 1. DATA CONSTRUCTION -- one setting from 00_setup.R per row of `variants`.
##    03 and 04 are re-run with it before refitting.
##
## 2. HOW ATTENTION ENTERS -- four specifications, fitted for every variant.
##    This is the axis that matters most, and the one where this analysis
##    departs from Gross et al. See docs/deviations.md 13 and the header of 05.
##
##      Model A / Model B    within/between split. att_cw2 is the within-person
##                           Present-minus-MW difference. PRIMARY.
##      Blended              attention entered raw with contr.sum, as in Gross's
##                           Table 3. The random intercept does not absorb the
##                           between-person part, so this is a blend of the
##                           within and between effects.
##      Centred outcome      wb_cw2 as the outcome, as in Gross's mediation
##                           chain. Isolates the within effect but zeroes the
##                           random-intercept variance, so it is singular by
##                           construction, and single-prompt participants
##                           contribute a zero outcome against a varying
##                           predictor. Reported, not used.
##
## The primary specification is fixed in 00_setup.R and 05, and does not change
## on the basis of what this table says -- the point is to report the
## alternatives, not to pick among them.
## -----------------------------------------------------------------------------

source("R/00_setup.R")

## The defaults, then one change per row.
variants <- tribble(
  ~label,                             ~effort, ~reml, ~centre,         ~age, ~drop_prompts,
  "Primary (preregistered)",                2L,  TRUE, "inner_speech",  18L, FALSE,
  "ML instead of REML",                     2L, FALSE, "inner_speech",  18L, FALSE,
  "Mood centred on all prompts",            2L,  TRUE, "all_moments",   18L, FALSE,
  "Strict effort screen (option 1 only)",   1L,  TRUE, "inner_speech",  18L, FALSE,
  "Prompt-level effort screen",             2L,  TRUE, "inner_speech",  18L, TRUE,
  "No age screen",                          2L,  TRUE, "inner_speech",   0L, FALSE
)

## label, formula, the coefficient to read, and whether it needs doubling.
## contr.sum puts +1 on Present, so `attention1` is HALF the difference;
## att_num is 0/1, so `att_cw2` is already the whole of it.
specs <- tribble(
  ~model,            ~rhs,                                             ~term,        ~double,
  "Model A",         "att_cw2 + att_cb2",                              "att_cw2",      FALSE,
  "Model B",         "att_cw2 + att_cb2 + valence_cw2",                "att_cw2",      FALSE,
  "Blended",         "attention",                                      "attention1",   TRUE,
  "Centred outcome", "attention",                                      "attention1",   TRUE
)
specs$lhs <- c("wb", "wb", "wb", "wb_cw2")

results <- list()

for (i in seq_len(nrow(variants))) {
  v <- variants[i, ]
  cat("\n----", v$label, "----------------------------------------\n")

  ## Apply this variant's settings, then re-run the two scripts that use them.
  ## They are noisy, so their printing is captured rather than shown.
  effort_pass_max        <<- v$effort
  use_reml               <<- v$reml
  centering_base         <<- v$centre
  min_age                <<- v$age
  drop_prompts_on_effort <<- v$drop_prompts

  invisible(capture.output({
    source("R/03_exclusions.R")
    source("R/04_derive.R")
  }))

  dv <- readRDS(file.path(dir_derived, "04_analysis.rds"))

  for (j in seq_len(nrow(specs))) {
    sp <- specs[j, ]
    f  <- as.formula(paste(sp$lhs, "~", sp$rhs,
                           "+ activity + clarity_cw2 + interesting_cw2 + (1 | PID)"))
    m  <- lmer(f, data = dv, REML = v$reml)
    co <- summary(m)$coefficients[sp$term, ]
    k  <- if (sp$double) 2 else 1

    results[[length(results) + 1]] <- tibble(
      variant      = v$label,
      model        = sp$model,
      participants = nlevels(dv$PID),
      prompts      = nrow(dv),
      ## always reported as the full Present - MW difference
      present_minus_mw = k * co["Estimate"],
      ci_low           = k * (co["Estimate"] - 1.96 * co["Std. Error"]),
      ci_high          = k * (co["Estimate"] + 1.96 * co["Std. Error"]),
      p                = co["Pr(>|t|)"],
      singular         = isSingular(m)
    )
    cat(sprintf("  %-16s Present - MW = %+.3f [%+.3f, %+.3f]  p = %.4f%s\n",
                sp$model, k * co["Estimate"],
                k * (co["Estimate"] - 1.96 * co["Std. Error"]),
                k * (co["Estimate"] + 1.96 * co["Std. Error"]), co["Pr(>|t|)"],
                if (isSingular(m)) "  SINGULAR" else ""))
  }
}

sensitivity <- bind_rows(results)

cat("\n---- All specifications ----------------------------------------------\n")
print(as.data.frame(sensitivity), digits = 3, row.names = FALSE)

cat("\nEvery 'Centred outcome' row is singular. That is the specification, not",
    "the data:\nperson-mean centring the outcome leaves (1 | PID) nothing to",
    "estimate.\n")

## Put the settings back, and rebuild the primary analysis data so the files on
## disk match 00_setup.R again rather than the last variant run above.
source("R/00_setup.R")
drop_prompts_on_effort <- FALSE
invisible(capture.output({ source("R/03_exclusions.R"); source("R/04_derive.R") }))
cat("\nRestored the primary specification and rebuilt 03/04 outputs.\n")

saveRDS(sensitivity, file.path(dir_derived, "07_sensitivity.rds"))
cat("Saved:", file.path(dir_derived, "07_sensitivity.rds"), "\n")
