## -----------------------------------------------------------------------------
## 04  Build the analysis data: keep inner-speech prompts, centre the
##     predictors, set the contrasts.
##
## Centering follows Gross et al.'s 'Person mean centering' chunk:
##   _cw2   within-person  (prompt minus that person's mean)
##   _c2    grand-mean centred
##   _cb2   between-person (= _c2 - _cw2)
##
## What differs is the BASE. Gross took person means over every prompt, using
## items that were asked of everyone. Here the thought items are only asked on
## inner-speech prompts, so the means are taken over those. That means each
## person's mean rests on roughly 13 prompts rather than 42, so the _cw2
## predictors are noisier than in the original. See docs/deviations.md.
## -----------------------------------------------------------------------------

if (!exists("read_expiwell")) source("R/00_setup.R")
inc <- readRDS(file.path(dir_derived, "03_included.rds"))

## Only inner-speech prompts. This also removes the handful of prompts flagged
## in 02 that said "no inner speech" but carried thought answers.
d <- inc$esm |> filter(inner_speech == "Yes")
cat("\nInner-speech prompts:", nrow(d), "from", n_distinct(d$PID), "participants\n")
cat("Prompts per person: ")
print(summary(as.integer(table(d$PID))))


## -- Person means and centring -----------------------------------------------

if (centering_base == "all_moments") {
  ## Mood is asked at every prompt, so its person mean can use them all.
  ## The thought items cannot -- they only exist here.
  cat("\nCentring mood on ALL prompts, thought items on inner-speech prompts.\n")
  wb_means <- inc$esm |> group_by(PID) |> summarise(wb_mean2 = mean(wb, na.rm = TRUE))
  d <- d |> left_join(wb_means, by = "PID")
} else {
  cat("\nCentring everything on inner-speech prompts.\n")
  d <- d |> group_by(PID) |> mutate(wb_mean2 = mean(wb, na.rm = TRUE)) |> ungroup()
}

d <- d |>
  group_by(PID) |>
  mutate(
    clarity_mean2     = mean(clarity, na.rm = TRUE),
    valence_mean2     = mean(valence, na.rm = TRUE),
    interesting_mean2 = mean(interesting, na.rm = TRUE)
  ) |>
  ungroup() |>
  mutate(
    ## within-person
    wb_cw2          = wb          - wb_mean2,
    clarity_cw2     = clarity     - clarity_mean2,
    valence_cw2     = valence     - valence_mean2,
    interesting_cw2 = interesting - interesting_mean2,
    ## grand-mean centred
    wb_c2           = wb          - mean(wb, na.rm = TRUE),
    clarity_c2      = clarity     - mean(clarity, na.rm = TRUE),
    valence_c2      = valence     - mean(valence, na.rm = TRUE),
    interesting_c2  = interesting - mean(interesting, na.rm = TRUE),
    ## between-person
    wb_cb2          = wb_c2          - wb_cw2,
    clarity_cb2     = clarity_c2     - clarity_cw2,
    valence_cb2     = valence_c2     - valence_cw2,
    interesting_cb2 = interesting_c2 - interesting_cw2
  )

## -- Attention state, decomposed the same way --------------------------------
## `att_num` is 1 at a Present prompt and 0 at a Mind-Wandering one, so `att_cw2`
## carries the FULL Present-minus-MW difference within a person. No doubling,
## unlike the contr.sum coding set up below.
##
## Why bother, when there is already a random intercept: `(1 | PID)` does not
## remove the between-person part of a level-1 predictor. Attention entered raw
## gives a precision-weighted blend of the within- and between-person effects
## (Neuhaus & Kalbfleisch, 1998; Enders & Tofighi, 2007; Curran & Bauer, 2011).
## Splitting it is what makes the H1 coefficient a within-person quantity, which
## is what the hypothesis is about. See docs/deviations.md 13.

d <- d |>
  mutate(att_num = as.numeric(attention == "Present")) |>
  group_by(PID) |>
  mutate(att_mean2 = mean(att_num, na.rm = TRUE)) |>
  ungroup() |>
  mutate(
    att_cw2 = att_num - att_mean2,                        # within-person
    att_c2  = att_num - mean(att_num, na.rm = TRUE),      # grand-mean centred
    att_cb2 = att_c2 - att_cw2                            # between-person
  )


cat("\n-- centring checks ---------------------------------------------------\n")
cat("largest person-mean of wb_cw2 (should be ~0):",
    max(abs(tapply(d$wb_cw2, d$PID, mean, na.rm = TRUE)), na.rm = TRUE), "\n")
cat("_c2 equals _cw2 + _cb2:", isTRUE(all.equal(d$wb_c2, d$wb_cw2 + d$wb_cb2)), "\n\n")

## A participant with a single prompt has att_cw2 exactly 0: their one prompt IS
## their mean. They therefore contribute nothing to the within-person slope while
## still informing the intercept and the variance components -- which is the
## right treatment. Under a person-mean-centred OUTCOME they are not neutral but
## actively harmful, contributing a zero outcome against a varying predictor.
n_single <- sum(table(d$PID)[as.character(d$PID)] == 1)
cat("single-prompt rows:", n_single, "of", nrow(d),
    "-- att_cw2 exactly 0 for all of them:",
    all(d$att_cw2[table(d$PID)[as.character(d$PID)] == 1] == 0), "\n\n")
print(summary(d[c("wb", "wb_cw2", "clarity_cw2", "valence_cw2", "interesting_cw2")]))


## -- Contrasts ---------------------------------------------------------------
## contr.sum(2) puts +1 on the first level. With Present first, `attention1` is
## positive when mood is higher while present -- the direction H1 predicts, and
## HALF the Present-minus-MW difference.
##
## Careful: R/power/power_gross2024.R let the levels sort alphabetically
## ("MW", "Present"), which gives the coefficient the opposite sign.

d$attention <- factor(as.character(d$attention), levels = attention_levels)
contrasts(d$attention) <- contr.sum(2)
cat("\ncontrasts(attention):\n"); print(contrasts(d$attention))

## A plain treatment-coded copy for the mediation model in 05, which needs an
## explicit control and treatment value rather than sum-to-zero contrasts.
d$attention_tx <- factor(as.character(d$attention), levels = c("MW", "Present"))

d$activity <- factor(as.character(d$activity), levels = activity_levels)
cat("\nactivity (reference =", activity_levels[1], "):\n")
print(table(d$activity, useNA = "ifany"))


## -- Add the person-level variables ------------------------------------------
d <- d |>
  left_join(
    inc$person |> select(PID, age, gender, gender3, ladder, n_inner,
                         starts_with("ffmq_")),
    by = "PID"
  ) |>
  mutate(
    ## centred so the attention main effect stays interpretable in 06
    age_c        = age        - mean(age, na.rm = TRUE),
    ladder_c     = ladder     - mean(ladder, na.rm = TRUE),
    ffmq_total_c = ffmq_total - mean(ffmq_total, na.rm = TRUE),
    PID          = factor(PID)
  )

cat("\n-- analysis data -----------------------------------------------------\n")
cat(nrow(d), "prompts,", nlevels(d$PID), "participants\n\n")
print(d |> select(PID, day, attention, att_cw2, att_cb2, wb, clarity_cw2,
                  valence_cw2, interesting_cw2, activity), n = 6)

saveRDS(d, file.path(dir_derived, "04_analysis.rds"))
cat("\nSaved:", file.path(dir_derived, "04_analysis.rds"), "\n")
