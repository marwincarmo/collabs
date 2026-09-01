## -----------------------------------------------------------------------------
## 08  Tables and figures for the paper. Nothing new is estimated here.
## -----------------------------------------------------------------------------

if (!exists("read_expiwell")) source("R/00_setup.R")
inc  <- readRDS(file.path(dir_derived, "03_included.rds"))
d    <- readRDS(file.path(dir_derived, "04_analysis.rds"))
mod  <- readRDS(file.path(dir_derived, "05_models.rds"))
sec  <- readRDS(file.path(dir_derived, "06_secondary.rds"))
sens <- readRDS(file.path(dir_derived, "07_sensitivity.rds"))


## -- Sample descriptives ------------------------------------------------------
person <- inc$person
cat("\n---- Sample ----------------------------------------------------------\n")
cat("Participants:", nrow(person), "\n")
cat("Age:", sprintf("M = %.1f, SD = %.1f, range %d-%d", mean(person$age, na.rm = TRUE),
                    sd(person$age, na.rm = TRUE), min(person$age, na.rm = TRUE),
                    max(person$age, na.rm = TRUE)), "\n")
print(table(person$gender3, useNA = "ifany"))
cat("\nSocial status ladder: M =", round(mean(person$ladder, na.rm = TRUE), 2), "\n")
cat("FFMQ (14-item mean):  M =", round(mean(person$ffmq_total, na.rm = TRUE), 2), "\n")
cat("Inner-speech prompts per person: median", median(person$n_inner),
    "range", min(person$n_inner), "-", max(person$n_inner), "\n")

descriptives <- d |>
  summarise(across(c(wb, clarity, valence, interesting),
                   list(n = ~sum(!is.na(.x)), mean = ~mean(.x, na.rm = TRUE),
                        sd = ~sd(.x, na.rm = TRUE)))) |>
  pivot_longer(everything(), names_to = c("variable", "stat"), names_sep = "_(?=[^_]+$)") |>
  pivot_wider(names_from = stat, values_from = value)

cat("\n---- Prompt-level descriptives ---------------------------------------\n")
print(as.data.frame(descriptives), digits = 3, row.names = FALSE)
cat("Present:", sprintf("%.1f%%", 100 * mean(d$attention == "Present")), "of prompts\n")


## -- Table 3 analogue ---------------------------------------------------------
tidy_fixed <- function(m, label) {
  co <- summary(m)$coefficients
  tibble(model = label, term = rownames(co),
         b = co[, "Estimate"], se = co[, "Std. Error"],
         df = co[, "df"], t = co[, "t value"], p = co[, "Pr(>|t|)"])
}
table3 <- bind_rows(tidy_fixed(mod$m_a, "Model A"), tidy_fixed(mod$m_b, "Model B"))

cat("\n---- Table 3: fixed effects ------------------------------------------\n")
print(as.data.frame(table3), digits = 3, row.names = FALSE)


## -- Figures ------------------------------------------------------------------

## Mood by attention state, one grey line per participant.
person_means <- d |>
  group_by(PID, attention) |>
  summarise(wb = mean(wb, na.rm = TRUE), .groups = "drop")

p_attention <- ggplot(person_means, aes(attention, wb)) +
  geom_line(aes(group = PID), colour = "grey75", linewidth = .3, alpha = .7) +
  geom_point(colour = "grey55", size = 1, alpha = .6) +
  stat_summary(fun = mean, geom = "point", size = 3.5, colour = "#1b6ca8") +
  stat_summary(fun.data = mean_se, geom = "errorbar", width = .12,
               linewidth = .8, colour = "#1b6ca8") +
  labs(x = "Attention state", y = "Mood (person mean, -3 to +3)",
       title = "Mood by attention state",
       subtitle = "Grey lines join each participant's two means") +
  theme_minimal(base_size = 12)

p_activity <- ggplot(d, aes(reorder(activity, wb, FUN = mean), wb)) +
  geom_boxplot(outlier.alpha = .3, fill = "grey93") +
  stat_summary(fun = mean, geom = "point", size = 2.5, colour = "#c0392b") +
  coord_flip() +
  labs(x = NULL, y = "Mood (-3 to +3)", title = "Mood by current activity",
       subtitle = "Red points are means") +
  theme_minimal(base_size = 12)

p_sensitivity <- ggplot(sens, aes(present_minus_mw, reorder(variant, present_minus_mw),
                                  colour = model)) +
  geom_vline(xintercept = 0, linetype = "dashed", colour = "grey50") +
  geom_pointrange(aes(xmin = ci_low, xmax = ci_high),
                  position = position_dodge(width = .5), size = .4) +
  labs(x = "Present - MW difference in mood", y = NULL, colour = NULL,
       title = "H1 across specifications") +
  theme_minimal(base_size = 12) + theme(legend.position = "bottom")

ggsave(file.path(dir_output, "fig_attention.png"),   p_attention,   width = 6, height = 5, dpi = 150)
ggsave(file.path(dir_output, "fig_activity.png"),    p_activity,    width = 7, height = 4, dpi = 150)
ggsave(file.path(dir_output, "fig_sensitivity.png"), p_sensitivity, width = 8, height = 5, dpi = 150)


## -- Write the tables out -----------------------------------------------------
write.csv(inc$flow,            file.path(dir_output, "flow.csv"), row.names = FALSE)
write.csv(descriptives,        file.path(dir_output, "descriptives.csv"), row.names = FALSE)
write.csv(table3,              file.path(dir_output, "table3_fixed_effects.csv"), row.names = FALSE)
write.csv(mod$pseudo_r2,       file.path(dir_output, "pseudo_r2.csv"), row.names = FALSE)
write.csv(sec$moderation,      file.path(dir_output, "moderation.csv"), row.names = FALSE)
write.csv(sec$activity_means,  file.path(dir_output, "activity_means.csv"), row.names = FALSE)
write.csv(sens,                 file.path(dir_output, "sensitivity.csv"), row.names = FALSE)

cat("\n---- Written to", dir_output, "-------------------------------------\n")
print(list.files(dir_output, pattern = "csv|png"))

cat("\nNOTE: these numbers come from partial data --",
    nlevels(d$PID), "participants with a median of",
    median(as.integer(table(d$PID))), "inner-speech prompt(s) each.\n")
cat("The prereg targets 200 (Model A) / 300 (Model B).\n")
