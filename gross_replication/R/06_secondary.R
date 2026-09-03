## -----------------------------------------------------------------------------
## 06  The secondary hypotheses from the prereg.
##
##  S1  Thought valence, clarity, interestingness and activity predict mood.
##  S2  FFMQ and demographics do NOT moderate the attention-mood relationship.
##  S3  Social activity shows the highest mood.
##  S4  Age moderates the attention-mood relationship.
##
## S2 is a predicted null. A non-significant interaction is not evidence that
## there is no moderation, so read the intervals, not the stars.
## -----------------------------------------------------------------------------

if (!exists("read_expiwell")) source("R/00_setup.R")
d   <- readRDS(file.path(dir_derived, "04_analysis.rds"))
mod <- readRDS(file.path(dir_derived, "05_models.rds"))


## == S1: do the covariates predict mood? ======================================
## Straight off Model B -- no new model needed.
cat("\n---- S1: covariate effects (Model B) ---------------------------------\n")
print(round(summary(mod$m_b)$coefficients[c("clarity_cw2", "interesting_cw2",
                                            "valence_cw2"), ], 4))


## == S3: is social activity the happiest? =====================================
## Social is the reference level, so every activity coefficient below reads as
## "this activity compared with social activity". S3 predicts all are negative.
cat("\n---- S3: activity vs. social activity (Model A) ----------------------\n")
print(round(summary(mod$m_a)$coefficients[grep("^activity", rownames(
  summary(mod$m_a)$coefficients)), ], 4))

activity_means <- d |>
  group_by(activity) |>
  summarise(prompts = n(), participants = n_distinct(PID),
            mean_wb = mean(wb, na.rm = TRUE), sd_wb = sd(wb, na.rm = TRUE),
            .groups = "drop") |>
  arrange(desc(mean_wb))

cat("\nRaw mood by activity:\n")
print(as.data.frame(activity_means), digits = 3, row.names = FALSE)


## == S2 and S4: moderation ====================================================
## Each moderator is added to Model A along with its interaction with attention.
## Only the interaction row matters for the hypothesis.

moderators <- c(
  "S2 FFMQ total"        = "ffmq_total_c",
  "S2 FFMQ observe"      = "ffmq_observe",
  "S2 FFMQ describe"     = "ffmq_describe",
  "S2 FFMQ act-aware"    = "ffmq_actaware",
  "S2 FFMQ non-judging"  = "ffmq_nonjudge",
  "S2 FFMQ non-reactive" = "ffmq_nonreact",
  "S2 gender"            = "gender3",
  "S2 social status"     = "ladder_c",
  "S4 age"               = "age_c"
)

moderation <- list()

for (i in seq_along(moderators)) {
  v     <- moderators[[i]]
  label <- names(moderators)[i]

  ## The moderated term is the WITHIN-person attention effect, att_cw2, to match
  ## the H1 specification in 05. att_cb2 stays in as a main effect so the split
  ## is not undone. att_cw2 is 0/1 coded, so the interaction is the change in the
  ## full Present-minus-MW difference per unit of the moderator.
  f <- as.formula(paste(
    "wb ~ att_cw2 * ", v,
    "+ att_cb2 + activity + clarity_cw2 + interesting_cw2 + (1 | PID)"
  ))
  m <- lmer(f, data = d, REML = use_reml)

  co  <- summary(m)$coefficients
  rows <- grep("att_cw2:", rownames(co), fixed = TRUE)

  cat("\n----", label, "-----------------------------------------------------\n")
  cat("formula:", deparse1(formula(m)), "\n\n")
  print(round(co[rows, , drop = FALSE], 4))
  if (isSingular(m)) cat("NOTE: singular fit\n")

  moderation[[label]] <- tibble(
    hypothesis = label,
    term       = rownames(co)[rows],
    estimate   = co[rows, "Estimate"],
    se         = co[rows, "Std. Error"],
    df         = co[rows, "df"],
    t          = co[rows, "t value"],
    p          = co[rows, "Pr(>|t|)"]
  )
}

moderation <- bind_rows(moderation) |>
  mutate(ci_low = estimate - 1.96 * se, ci_high = estimate + 1.96 * se)

cat("\n---- All interaction terms together ----------------------------------\n")
print(as.data.frame(moderation |> select(hypothesis, term, estimate, ci_low, ci_high, p)),
      digits = 3, row.names = FALSE)

saveRDS(list(moderation = moderation, activity_means = activity_means),
        file.path(dir_derived, "06_secondary.rds"))
cat("\nSaved:", file.path(dir_derived, "06_secondary.rds"), "\n")
