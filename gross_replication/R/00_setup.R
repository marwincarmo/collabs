## Packages, paths and analysis settings. Sourced at the top of every script.
## Nothing here does any work -- it just sets things up.

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(readr)
  library(stringr)
  library(lme4)
  library(lmerTest)
  library(ggplot2)
})

## -- Paths -------------------------------------------------------------------
raw_person_file <- "data/raw/Primeiro_dia.csv"
raw_esm_file    <- "data/raw/Questionário.csv"
raw_post_file   <- NA          # set this when/if a final-day export arrives
dir_derived     <- "data/derived"
dir_output      <- "output"

## -- Analysis settings -------------------------------------------------------
## Reasons for each are in docs/deviations.md.

effort_pass_max  <- 2          # effort item: 1-2 pass, 3-4 fail. 1 = strict
min_age          <- 18         # consent restricts the study to adults
use_reml         <- TRUE       # the prereg specifies REML

## Person means are taken over inner-speech prompts only, because the thought
## items are not asked on the other prompts. "all_moments" centres mood over
## every prompt instead, which is closer to Gross; 07 reports both.
centering_base   <- "inner_speech"
mediation_sims   <- 5000
seed             <- 20260504   # OSF registration date

## contr.sum(2) puts +1 on the FIRST level, so with Present first the
## `attention1` coefficient is positive when mood is higher while present.
attention_levels <- c("Present", "MW")

## Social is the reference level, so each activity coefficient reads as
## "this activity vs. social activity".
activity_levels  <- c("Social", "Cognitive", "Household", "Other", "Physical", "Restful")

set.seed(seed)
for (d in c(dir_derived, dir_output)) if (!dir.exists(d)) dir.create(d, recursive = TRUE)


## -- Reading the Expiwell exports --------------------------------------------
## The exports are not plain CSVs. There are four header rows above the data:
##
##   row 1  question ids           {{Q://Q7/ChoiceTextEntryValue}}
##   row 2  the question text
##   row 3  response type          "Single selection" / "Text"
##   row 4  the real column names, then the literal "Choices", then, for each
##          question, its answer options as  [ 1 = 'label' 2 = 'label' ]
##   row 5+ the data, stored as LABEL TEXT rather than numeric codes
##
## This is the one thing worth wrapping in a function, because it is fiddly and
## we do it twice. It returns the data plus a small codebook table.

read_expiwell <- function(path) {
  raw <- read_csv(path, col_names = FALSE, na = character(),
                  col_types = cols(.default = col_character()), progress = FALSE)

  ## Squish whitespace and drop stray byte-order marks; they turn up mid-file.
  tidy_txt <- function(x) str_squish(gsub("\ufeff", "", x, fixed = TRUE))

  qid   <- tidy_txt(as.character(unlist(raw[1, ])))
  qtext <- tidy_txt(as.character(unlist(raw[2, ])))
  qtype <- tidy_txt(as.character(unlist(raw[3, ])))
  hdr   <- tidy_txt(as.character(unlist(raw[4, ])))

  ## Everything left of the "Choices" cell is metadata; everything right of it
  ## is a question. That cell holds the participant label in the data rows.
  label_col <- which(hdr == "Choices")
  meta_cols <- seq_len(label_col - 1)
  q_cols    <- seq(label_col + 1, ncol(raw))

  dat <- raw[-(1:4), ]
  names(dat)[meta_cols] <- hdr[meta_cols]
  names(dat)[label_col] <- "participant_label"
  names(dat)[q_cols]    <- str_match(qid[q_cols], "Q\\d+")[, 1]
  dat <- dat |> mutate(across(everything(), ~ na_if(tidy_txt(.x), "")))

  ## Pull each question's answer options out of the row-4 choice map, so we
  ## recode against the survey's own definitions rather than typed-out strings.
  choices <- lapply(hdr[q_cols], function(s) {
    m <- str_match_all(s, "(\\d+)\\s*=\\s*'(.*?)'(?=\\s*\\d+\\s*=\\s*'|\\s*\\]\\s*$)")[[1]]
    if (nrow(m) == 0) NULL else tidy_txt(m[, 3])[order(as.integer(m[, 2]))]
  })
  names(choices) <- names(dat)[q_cols]

  codebook <- tibble(
    column    = names(dat)[q_cols],
    n_options = vapply(choices, function(x) if (is.null(x)) NA_integer_ else length(x), integer(1)),
    type      = qtype[q_cols],
    question  = str_trunc(qtext[q_cols], 70)
  )

  list(data = dat, choices = choices, codebook = codebook)
}
