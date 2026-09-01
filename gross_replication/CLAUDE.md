# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project

A Brazilian replication of Gross et al. (2024/2025) — an experience-sampling (ESM)
study of inner speech, attention state (Present vs. Mind-Wandering) and momentary
well-being. Preregistered at OSF `txfsv` (registered 2026-05-04); PUC-SP; target
N = 200 (Model A) / 300 (Model B), 7 ESM days at 6 prompts/day after a day-1 survey.

`gross_replication/` is one project inside the `collabs/` git repository (siblings:
`afghanistan/`, `concussion_uk/`). The RStudio project file (`collabs.Rproj`) is at
the repo root, one level up. **Run everything with `gross_replication/` as the
working directory.**

## How the analysis is organised

Numbered scripts, run top to bottom, in order. Each one prints what it did, saves
its result to `data/derived/`, and is picked up by the next. Run them one at a
time and read the output — that is the point. There is no wrapper script and no
single entry point.

```r
source("R/01_read.R")        # what is in the export: questions, options, counts
source("R/02_recode.R")      # label text -> numbers, with a check table per item
source("R/03_exclusions.R")  # prereg exclusions, one at a time + flow table
source("R/04_derive.R")      # inner-speech prompts, centring, contrasts
source("R/05_models.R")      # H1 and H2
source("R/06_secondary.R")   # covariates, activity, moderators
source("R/07_sensitivity.R") # H1 under each contested choice
source("R/08_tables_figures.R")
```

`R/00_setup.R` holds packages, file paths and the analysis settings, and is
sourced by each script. `R/09_check_recovery.R` is a standalone sanity check:
it simulates data with known effects and confirms the model recovers them.

To analyse the final export, change the paths at the top of `R/00_setup.R`.

Scripts 03 and 04 begin with `if (!exists("read_expiwell")) source("R/00_setup.R")`
so that `07_sensitivity.R` can change a setting and re-run them without the setup
file resetting it. That is the only bit of indirection in the pipeline.

No build or lint harness, no `renv`, no `testthat`. `gt` is not installed.

## Things about the data that will bite

**The exports are not plain CSVs.** Four header rows sit above the data: question
ids, question text, response type, then the real column names followed by each
question's answer options as `[ 1 = 'label' 2 = 'label' ]`. `read_expiwell()` in
`00_setup.R` parses this. Answers are stored as **label text**, not codes, and
inconsistently: scale endpoints carry words (`"+3 (muito bem)"`) while midpoints
are bare (`"5"`).

**The option code is not the analysis value.** Mood and valence are options 1–7
for a −3..+3 scale, so `02_recode.R` subtracts 4 — the same shift Gross applied.
Recoding matches against the survey's own option list rather than typed-out
strings, and each recode prints a raw-label × value table plus a check that
nothing became `NA`.

**`attention1` is half the Present − MW difference.** `contr.sum(2)` puts +1 on
the first level, and `00_setup.R` sets `c("Present", "MW")` so the coefficient is
positive in the direction H1 predicts. `R/power/` let the levels sort
alphabetically, which gives the opposite sign — the two are not comparable.
Scripts 05 and 07 report the doubled difference explicitly.

**The instrument differs from Gross's in ways that matter.** There is no
no-inner-speech branch (so person-mean centring uses a different base), no
`reaction` item, and the FFMQ has 14 items where the prereg names the FFMQ-15.
`docs/deviations.md` lists all twelve departures, each tied to a setting in
`00_setup.R`. Read it before changing a model or an exclusion rule.

**`Occasion within Day` is constant** in the pilot export, so prompt order comes
from `Start Date`. **Three prompts** say "no inner speech" yet carry thought
answers; `02_recode.R` flags them and `04_derive.R` drops them with the rest of
the non-inner-speech prompts.

## Data

- `data/raw/` — the exports as received. **Gitignored**: they contain
  `ip address`, `Location - Lat` and `Location - Long`. `02_recode.R` keeps only
  the columns it needs, so nothing in `output/` has them, but the raw files do.
  Strip them before any public deposit.
- `data/derived/` — one `.rds` per script (gitignored).
- `osfstorage-archive/` — the original authors' OSF download. `ESM cleaning.Rmd`
  and `ESM analysis.Rmd` are their own scripts; when a decision here looks
  arbitrary, the chunk named in the comment is where it came from.

## `R/power/` — the power analysis (frozen)

The simulation that justified the target N, kept as preregistration evidence. Do
not change it to match the analysis scripts; they legitimately differ (ML vs
REML, contrast sign), and `docs/deviations.md` records why.

`power_gross2024.R` is meant to be stepped through, not sourced: section 3 refits
two `brms`/cmdstanr models and runs 28 scenarios × 1000 simulations (hours), with
its `saveRDS` commented out and the cached `output/power_analysis_hybrid.rds`
read back instead. Section 4 is leftover and errors if sourced.
`hybrid_test.R` and `sim_main.R` are superseded scratch files.

`power_simulation_gross_repl.R` (repo root of this project) is the earlier fully
parametric simulation (J. Bottesini, June 2025) — superseded, but it documents
the missingness and attention-state assumptions. It uses **repo-root-relative**
paths, unlike everything in `R/`.

## Status

Data collection is not finished. The current numbers come from 95 participants
with a median of one inner-speech prompt each, which is why the mediation chain
reports singular fits. The specifications are fixed in `00_setup.R` and the
sensitivity table reports alternatives rather than choosing among them.
