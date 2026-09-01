## -----------------------------------------------------------------------------
## 03  Apply the preregistered exclusions, one at a time.
##
## From the prereg: exclude people who said they did not take the study
## seriously, people who failed the attention checks, and people with no
## inner-speech prompt at all. The Brazilian survey has a single effort item per
## survey rather than Gross's separate validity and attention checks, so the
## first two collapse into one screen.
##
## Each step prints how many participants and prompts are left. The table at the
## end is the participant flow figure for the paper.
## -----------------------------------------------------------------------------

## Re-sourcing setup would undo the settings 07 changes when it re-runs this
## script, so only load it when it has not been loaded already.
if (!exists("read_expiwell")) source("R/00_setup.R")
rec <- readRDS(file.path(dir_derived, "02_recoded.rds"))

person <- rec$person
esm    <- rec$esm

## Small helper so each step below is one readable line.
count_now <- function(person, esm) c(
  participants = n_distinct(person$PID),
  prompts      = sum(esm$PID %in% person$PID),
  inner        = sum(esm$PID %in% person$PID & esm$inner_speech == "Yes", na.rm = TRUE)
)
steps <- list(`Enrolled (did the day-1 survey)` = count_now(person, esm))


## -- Did they respond at all? ------------------------------------------------
## Not an exclusion, but the paper needs the number.
never_responded <- setdiff(person$PID, esm$PID)
cat("\nEnrolled but never answered a prompt:", length(never_responded), "\n")
person <- person |> filter(PID %in% esm$PID)
steps[["Answered >= 1 prompt"]] <- count_now(person, esm)


## -- 1. Effort / seriousness -------------------------------------------------
## The item runs 1 = read everything carefully, 2 = went quickly but answered
## correctly, 3 = clicked some answers without reading, 4 = rushed and did not
## answer properly. Options 3 and 4 say the person did not take it seriously;
## option 2 says the opposite. So the default cutoff keeps 1 and 2.
##
## Applied at participant level. The per-prompt effort item is asked only on the
## inner-speech branch, so screening prompts on it would drop prompts in a way
## that is tied to the outcome. 07 runs that stricter version.

cat("\nEffort check (cutoff: keep options 1-", effort_pass_max, ")\n", sep = "")
print(table(day1 = person$effort_person, useNA = "ifany"))
if (any(!is.na(person$effort_post))) {
  cat("A final-day check is also available and is applied here.\n")
  print(table(final_day = person$effort_post, useNA = "ifany"))
} else {
  cat("No final-day check available; using the day-1 answer only.\n")
}

failed_effort <- person$PID[
  (!is.na(person$effort_person) & person$effort_person > effort_pass_max) |
  (!is.na(person$effort_post)   & person$effort_post   > effort_pass_max)
]
cat("Excluded:", length(failed_effort), "\n")
person <- person |> filter(!PID %in% failed_effort)
steps[[paste0("Passed effort check (1-", effort_pass_max, ")")]] <- count_now(person, esm)


## -- 2. Age ------------------------------------------------------------------
## Not in the prereg, but consent restricts the study to adults and Gross
## applied the same screen. 07 reports the analysis without it.
too_young <- person$PID[!is.na(person$age) & person$age < min_age]
cat("\nUnder", min_age, ":", length(too_young), "\n")
person <- person |> filter(!PID %in% too_young)
steps[[paste0("Aged >= ", min_age)]] <- count_now(person, esm)


## -- 3. At least one inner-speech prompt -------------------------------------
## The analysis only uses inner-speech prompts, so someone with none of them
## contributes nothing.
inner_counts <- esm |>
  filter(inner_speech == "Yes") |>
  count(PID, name = "n_inner")

person <- person |>
  left_join(inner_counts, by = "PID") |>
  mutate(n_inner = replace_na(n_inner, 0L))

cat("\nInner-speech prompts per participant:\n")
print(table(person$n_inner))
cat("With none at all:", sum(person$n_inner == 0), "\n")
person <- person |> filter(n_inner > 0)
steps[["Has >= 1 inner-speech prompt"]] <- count_now(person, esm)


## -- Keep only the prompts belonging to retained participants ----------------
esm <- esm |> filter(PID %in% person$PID)

## Optional stricter screen, off by default. Used by 07.
if (exists("drop_prompts_on_effort") && isTRUE(drop_prompts_on_effort)) {
  before <- nrow(esm)
  esm <- esm |> filter(is.na(effort_prompt) | effort_prompt <= effort_pass_max)
  cat("\nPrompt-level effort screen dropped", before - nrow(esm), "prompts\n")
  steps[["Prompt-level effort screen"]] <- count_now(person, esm)
}


## -- Flow table --------------------------------------------------------------
flow <- as_tibble(do.call(rbind, steps)) |>
  mutate(step = names(steps), .before = 1) |>
  mutate(dropped = c(NA, -diff(participants)))

cat("\n---- Participant flow ------------------------------------------------\n")
print(as.data.frame(flow), row.names = FALSE)

saveRDS(list(person = person, esm = esm, flow = flow,
             never_responded = never_responded),
        file.path(dir_derived, "03_included.rds"))
cat("\nSaved:", file.path(dir_derived, "03_included.rds"), "\n")
