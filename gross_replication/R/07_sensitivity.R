## -----------------------------------------------------------------------------
## 07  Does the H1 result depend on the choices we had to make?
##
## Each variant changes ONE setting from 00_setup.R, then re-runs 03 and 04 with
## it and refits. The primary specification is fixed in 00_setup.R and does not
## change on the basis of what this table says -- the point is to report the
## alternatives, not to pick among them.
##
## The three models per variant:
##   Model A / Model B  -- outcome is raw mood
##   Total effect       -- outcome is wb_cw2, the model the mediation chain uses
## The centring variant only moves the third one, because it is the only model
## whose outcome is person-mean centred.
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

  fits <- list(
    "Model A"      = lmer(wb ~ attention + activity + clarity_cw2 + interesting_cw2 +
                            (1 | PID), data = dv, REML = v$reml),
    "Model B"      = lmer(wb ~ attention + activity + clarity_cw2 + interesting_cw2 +
                            valence_cw2 + (1 | PID), data = dv, REML = v$reml),
    "Total effect" = lmer(wb_cw2 ~ attention + activity + clarity_cw2 +
                            interesting_cw2 + (1 | PID), data = dv, REML = v$reml)
  )

  for (nm in names(fits)) {
    co <- summary(fits[[nm]])$coefficients["attention1", ]
    results[[length(results) + 1]] <- tibble(
      variant      = v$label,
      model        = nm,
      participants = nlevels(dv$PID),
      prompts      = nrow(dv),
      ## reported as the full Present - MW difference, i.e. twice the coefficient
      present_minus_mw = 2 * co["Estimate"],
      ci_low           = 2 * (co["Estimate"] - 1.96 * co["Std. Error"]),
      ci_high          = 2 * (co["Estimate"] + 1.96 * co["Std. Error"]),
      p                = co["Pr(>|t|)"],
      singular         = isSingular(fits[[nm]])
    )
    cat(sprintf("  %-13s Present - MW = %+.3f [%+.3f, %+.3f]  p = %.4f\n",
                nm, 2 * co["Estimate"],
                2 * (co["Estimate"] - 1.96 * co["Std. Error"]),
                2 * (co["Estimate"] + 1.96 * co["Std. Error"]), co["Pr(>|t|)"]))
  }
}

sensitivity <- bind_rows(results)

cat("\n---- All specifications ----------------------------------------------\n")
print(as.data.frame(sensitivity), digits = 3, row.names = FALSE)

## Put the settings back, and rebuild the primary analysis data so the files on
## disk match 00_setup.R again rather than the last variant run above.
source("R/00_setup.R")
drop_prompts_on_effort <- FALSE
invisible(capture.output({ source("R/03_exclusions.R"); source("R/04_derive.R") }))
cat("\nRestored the primary specification and rebuilt 03/04 outputs.\n")

saveRDS(sensitivity, file.path(dir_derived, "07_sensitivity.rds"))
cat("Saved:", file.path(dir_derived, "07_sensitivity.rds"), "\n")
