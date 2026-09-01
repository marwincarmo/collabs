## -----------------------------------------------------------------------------
## 02  Turn the label text into numbers.
##
## The exports store answers as label text ("+3 (muito bem)", "4 (moderadamente
## claro)"), not codes. Every recode below works by matching the answer against
## the option list the survey itself defines (printed in 01), then converting
## the option number to the analysis value.
##
## Each recode prints a table of raw label against recoded value. Read them --
## that table is the check. If the survey ever changes, a row will show NA.
## -----------------------------------------------------------------------------

source("R/00_setup.R")
raw <- readRDS(file.path(dir_derived, "01_raw.rds"))

e  <- raw$esm$data       # ESM prompts, one row per prompt
ch <- raw$esm$choices    # the ESM answer options, in survey order
p  <- raw$person$data    # day-1 survey, one row per participant
pc <- raw$person$choices


## A label that is not in the survey's option list would silently become NA.
## Each recode is checked for that immediately, before any rows are reordered.
check_na <- function(new, old, what) {
  bad <- sum(is.na(new) & !is.na(old))
  cat(sprintf("   [check] %-14s %s\n", what,
              if (bad == 0) "every answer recoded" else paste(bad, "FAILED TO RECODE")))
  stopifnot(bad == 0)
}


## == ESM prompts ==============================================================

esm <- tibble(
  PID       = e$`Participant ID`,
  day       = as.integer(e$`Day of Survey`),
  start_at  = as.POSIXct(e$`Start Date`, format = "%m/%d/%Y %I:%M%p", tz = "UTC"),
  finished  = e$Finished
)

## Attention state. Option 1 is "PRESENTE...", option 2 is "DIVAGANDO...".
esm$attention <- factor(c("Present", "MW")[match(e$Q1, ch$Q1)], levels = attention_levels)
cat("\n-- attention ---------------------------------------------------------\n")
print(table(str_trunc(e$Q1, 40), esm$attention, useNA = "ifany"))
check_na(esm$attention, e$Q1, "attention")

## Mood. Presented as options 1-7 for a -3..+3 scale, so subtract 4. This is the
## same shift Gross et al. applied (`wb - 4`).
esm$wb <- match(e$Q2, ch$Q2) - 4
cat("\n-- wb (mood) ---------------------------------------------------------\n")
print(table(e$Q2, esm$wb, useNA = "ifany"))
check_na(esm$wb, e$Q2, "wb")

## Inner speech. Option 1 = experiencing it, option 2 = not.
## NOTE: Gross recoded this to 0 = inner speech, 1 = none, which is why their
## script filters on `inner_speech == 0`. Here it stays a plain Yes/No.
esm$inner_speech <- factor(c("Yes", "No")[match(e$Q3, ch$Q3)], levels = c("Yes", "No"))
cat("\n-- inner_speech ------------------------------------------------------\n")
print(table(str_trunc(e$Q3, 40), esm$inner_speech, useNA = "ifany"))
check_na(esm$inner_speech, e$Q3, "inner_speech")

## Thought-content items. These are asked ONLY when the participant reports
## inner speech, so they are missing by design on the other prompts.
esm$clarity     <- match(e$Q4, ch$Q4)          # 1-7
esm$valence     <- match(e$Q5, ch$Q5) - 4      # -3..+3
esm$interesting <- match(e$Q6, ch$Q6)          # 1-7
cat("\n-- clarity / valence / interesting -----------------------------------\n")
print(table(e$Q4, esm$clarity, useNA = "ifany"))
print(table(e$Q5, esm$valence, useNA = "ifany"))
check_na(esm$clarity, e$Q4, "clarity")
check_na(esm$valence, e$Q5, "valence")
check_na(esm$interesting, e$Q6, "interesting")

## Current activity. The option order is the same as Gross's numeric coding:
## 1 social, 2 physical, 3 restful, 4 household, 5 cognitive, 6 other.
activity_by_option <- c("Social", "Physical", "Restful", "Household", "Cognitive", "Other")
esm$activity <- factor(activity_by_option[match(e$Q7, ch$Q7)], levels = activity_levels)
cat("\n-- activity ----------------------------------------------------------\n")
print(table(str_trunc(e$Q7, 35), esm$activity, useNA = "ifany"))
check_na(esm$activity, e$Q7, "activity")

## Per-prompt effort check, 1 = read carefully ... 4 = rushed.
esm$effort_prompt <- match(e$Q8, ch$Q8)
cat("\n-- effort_prompt -----------------------------------------------------\n")
print(table(str_trunc(e$Q8, 45), esm$effort_prompt, useNA = "ifany"))
check_na(esm$effort_prompt, e$Q8, "effort_prompt")


## -- Prompt order --------------------------------------------------------
## `Occasion within Day` is constant in the pilot export, so the order within
## each day comes from the timestamp instead. If a later export fills the
## column in properly, compare the two before trusting this.
esm <- esm |>
  arrange(PID, start_at) |>
  group_by(PID, day) |>
  mutate(occasion = row_number()) |>
  group_by(PID) |>
  mutate(prompt_index = row_number(), n_prompts = n()) |>
  ungroup()

cat("\n-- prompts per participant -------------------------------------------\n")
print(summary(esm$n_prompts[!duplicated(esm$PID)]))

## -- Branch check --------------------------------------------------------
## The thought items should be blank whenever inner_speech is not "Yes", but
## the platform does not enforce it: someone who goes back and changes their
## inner-speech answer can leave answers behind. Flag those and treat the
## inner-speech answer as the true one -- 04 filters them out.
esm <- esm |>
  mutate(branch_odd = (is.na(inner_speech) | inner_speech == "No") &
           (!is.na(clarity) | !is.na(valence) | !is.na(interesting)))
cat("\n-- prompts saying 'no inner speech' but carrying thought answers ------\n")
cat(sum(esm$branch_odd), "of", nrow(esm), "prompts\n")
if (any(esm$branch_odd)) {
  print(esm |> filter(branch_odd) |>
          select(PID, day, inner_speech, wb, clarity, valence, interesting))
}


## == Day-1 survey =============================================================

person <- tibble(
  PID      = p$`Participant ID`,
  start_at = as.POSIXct(p$`Start Date`, format = "%m/%d/%Y %I:%M%p", tz = "UTC")
)

## FFMQ. Q5-Q18 are the fourteen items, all on the same 1-5 scale.
## NOTE: the prereg names the FFMQ-15. This survey has 14 items -- three each
## for four facets but only two for non-reactivity. See docs/deviations.md.
ffmq_q <- paste0("Q", 5:18)
ffmq   <- sapply(ffmq_q, function(q) match(p[[q]], pc[[q]]))
colnames(ffmq) <- sprintf("ffmq_%02d", 1:14)

cat("\n-- FFMQ item responses (1-5) -----------------------------------------\n")
print(apply(ffmq, 2, table, useNA = "ifany"))

## Facet of each item, and whether it is reverse scored.
ffmq_facet <- c("Observe", "Describe", "ActAware", "NonJudge",
                "Observe", "Describe", "ActAware", "NonJudge",
                "NonReact", "Observe", "Describe", "ActAware",
                "NonJudge", "NonReact")
ffmq_rev   <- c(FALSE, FALSE, TRUE,  TRUE,  FALSE, TRUE,  TRUE,
                TRUE,  FALSE, FALSE, FALSE, TRUE,  TRUE,  FALSE)

ffmq_scored <- ffmq
ffmq_scored[, ffmq_rev] <- 6 - ffmq_scored[, ffmq_rev]     # 1-5 scale, so 6 - x

cat("\n-- FFMQ facets -------------------------------------------------------\n")
print(data.frame(item = colnames(ffmq), facet = ffmq_facet, reversed = ffmq_rev))

for (f in unique(ffmq_facet)) {
  person[[paste0("ffmq_", tolower(f))]] <- rowMeans(ffmq_scored[, ffmq_facet == f], na.rm = TRUE)
}
## Mean of the 14 administered items. Not comparable to a published FFMQ-15.
person$ffmq_total <- rowMeans(ffmq_scored, na.rm = TRUE)

## Demographics.
person$age <- suppressWarnings(as.integer(p$Q19))
cat("\n-- age ---------------------------------------------------------------\n")
print(summary(person$age))
cat("under 18:", sum(person$age < min_age, na.rm = TRUE), "\n")

person$gender <- factor(p$Q20, levels = pc$Q20)
cat("\n-- gender ------------------------------------------------------------\n")
print(table(person$gender, useNA = "ifany"))

## Collapsed for the moderation models: the non-cis and undisclosed cells are
## too small to estimate separately.
person$gender3 <- factor(case_when(
  p$Q20 == pc$Q20[1] ~ "Man",
  p$Q20 == pc$Q20[2] ~ "Woman",
  TRUE               ~ "Other / undisclosed"
), levels = c("Woman", "Man", "Other / undisclosed"))
cat("\ncollapsed for moderation models:\n")
print(table(person$gender3, useNA = "ifany"))

person$ladder <- match(p$Q21, pc$Q21)                # MacArthur ladder, 1-10
person$effort_person <- match(p$Q22, pc$Q22)         # 1 = careful ... 4 = rushed
cat("\n-- effort_person -----------------------------------------------------\n")
print(table(str_trunc(p$Q22, 45), person$effort_person, useNA = "ifany"))
check_na(person$effort_person, p$Q22, "effort_person")
check_na(person$gender, p$Q20, "gender")
check_na(person$ladder, p$Q21, "ladder")

## A final-day export, if one ever arrives, contributes a second effort check.
person$effort_post <- NA_integer_
if (!is.na(raw_post_file)) {
  post <- read_expiwell(raw_post_file)
  idx  <- match(person$PID, post$data$`Participant ID`)
  person$effort_post <- match(post$data$Q1, post$choices$Q1)[idx]
  cat("final-day effort check joined for", sum(!is.na(person$effort_post)), "participants\n")
} else {
  cat("\nNo final-day export configured; exclusions in 03 use the day-1 check only.\n")
}

## One row per participant. One duplicate in the pilot -- keep the fuller record.
cat("\n-- duplicate day-1 records -------------------------------------------\n")
cat(sum(duplicated(person$PID)), "duplicate(s)\n")
person <- person |>
  mutate(n_answered = rowSums(!is.na(across(starts_with("ffmq_"))))) |>
  arrange(PID, desc(n_answered), start_at) |>
  distinct(PID, .keep_all = TRUE) |>
  select(-n_answered)


## == Final look ===============================================================

cat("\n-- value ranges ------------------------------------------------------\n")
print(sapply(esm[c("wb", "clarity", "valence", "interesting")], range, na.rm = TRUE))

saveRDS(list(esm = esm, person = person), file.path(dir_derived, "02_recoded.rds"))
cat("\nSaved:", file.path(dir_derived, "02_recoded.rds"), "\n")
