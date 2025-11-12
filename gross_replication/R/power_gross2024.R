## 0 Load packages ---------------------------------------

library(lme4)
library(dplyr)
library(readr)
library(purrr)
library(future)
library(furrr)
library(ggplot2)
library(mediation)
source("R/sim_main.R") # load the simulation function

## 1 Load and prepare data -------------------------------

esm <- readr::read_csv("osfstorage-archive/ESM.Clean.Nov2022.csv")

## person-mean centering
esm_pc <- esm |>
  ## extracted from the original analysis
  ## create inner speach variables
  dplyr::mutate(clarity_all = coalesce(clarity, no_inner_speech_clarity),
                valence_all = coalesce(valence, no_inner_speech_valence),
                reaction_all = coalesce(reaction, no_inner_speech_reaction),
                interesting_all = coalesce(interesting, no_inner_speech_interesting),
                ## center variables
                wb = wb - 4,
                valence = valence - 4,
                valence_all = valence_all - 4,
                no_inner_speech_valence = no_inner_speech_valence - 4,
                ## recode activity
                activity = dplyr::case_when(
                                    activity == 1 ~ "Social",
                                    activity == 2 ~ "Physical",
                                    activity == 3 ~ "Restful",
                                    activity == 4 ~ "Household",
                                    activity == 5 ~ "Cognitive",
                                    activity == 6 ~ "Other",
                                    .default = NA
                                  ),
                attention  = factor(ifelse(attention == 1, "Present", "MW"))
                ) |>
  ## Person-mean centering (from chunk 'Person mean centering')
  dplyr::group_by(PID) |>
  dplyr::mutate(
           wb_mean2 = mean(wb, na.rm = T),
           clarity_all_mean2 = mean(clarity_all, na.rm = T),
           valence_all_mean2 = mean(valence_all, na.rm = T),
           reaction_all_mean2 = mean(reaction_all, na.rm = T),
           interesting_all_mean2 = mean(interesting_all, na.rm = T)
         ) |>
  dplyr::ungroup() |>
  ## Within-person centered variables
  dplyr::mutate(
           wb_cw2 = wb - wb_mean2,
           clarity_all_cw2 = clarity_all - clarity_all_mean2,
           valence_all_cw2 = valence_all - valence_all_mean2,
           reaction_all_cw2 = reaction_all - reaction_all_mean2,
           interesting_all_cw2 = interesting_all - interesting_all_mean2,
           ## Grand mean centered variables
           wb_c2 = scale(wb, center = T, scale = F),
           clarity_all_c2 = scale(clarity_all, center = T, scale = F),
           valence_all_c2 = scale(valence_all, center = T, scale = F),
           reaction_all_c2 = scale(reaction_all, center = T, scale = F),
           interesting_all_c2 = scale(interesting_all, center = T, scale = F),
           ## Betwee-person centered variables
           wb_cb2 = wb_c2-wb_cw2,
           clarity_all_cb2 = clarity_all_c2-clarity_all_cw2,
           valence_all_cb2 = valence_all_c2-valence_all_cw2,
           reaction_all_cb2 = reaction_all_c2-reaction_all_cw2,
           interesting_all_cb2 = interesting_all_c2-interesting_all_cw2
         )

## Create the "Inner Speech" only dataset (from chunk 'Inner speech datasets')
esm_inner <- esm_pc |>
  dplyr::filter(inner_speech == 0) # Inner Speech as they have on the paper

## Convert attention to a factor for the mediation package (from chunk 'With covariates')
## Note: The paper and script chunk 'As individual MLMs' use -1/1 coding.
## But the final mediation chunk 'With covariates' uses as.factor(). We MUST follow
## the final mediation chunk for our power analysis.
#esm_inner$attention <- as.factor(esm_inner$attention)
contrasts(esm_inner$attention) <- contr.sum(2)

## Get prompt counts for variable-prompt simulation
prompt_counts <- esm_inner |>
  dplyr::count(PID) |>
  dplyr::pull(n)


## 2 Fit the two key models from the script --------------

## Main hypothesis: Is Thought Valence more negative when people
## are mind wandering versus present?

## Table 3 - left panel ------------------------------------

## These are from the chunk 'With covariates' under '##Model 2: Mediation of Valeince'

model1a_fit <- lmerTest::lmer(formula = wb ~ activity + clarity_all_cw2 +
                    interesting_all_cw2 + attention + (1|PID),
                    data = esm_inner, REML = FALSE)
b_main_1a <- summary(model1a_fit)$coefficients[,"Estimate"]
vc_1a <- as.data.frame(VarCorr(model1a_fit))
sd_u0_1a <- vc_1a[1, "sdcor"]
sd_e_1a  <- vc_1a[2, "sdcor"]

## Model a: Mediator model (Valence ~ Attention + Covariates)
## (Corresponds to 'fit.mediator1' in the script)
## Does Attention State predict  Thought Valence (Mediator)?

## Table 3 reports  the coefficient for attention as 0.20. This value is only obtained
## when attention is not transformed to factor. However, according to their analysis plan,
## attention should be coded as Present(2) = 1, and Mind wandering(1) = -1

## Table 3 - right panel -----------------------------------

## This model is not in the original script, but corresponds to the results
## shown in the righ-hand panel of Table 3.
## The original analysis script shows a model fit.dv1 that contains
## within-person-centered valence, but the outcome is within-person
## centered well-beign, which deviates from the model shown in the paper

model1b_fit <- lmerTest::lmer(formula = wb ~ activity + clarity_all_cw2 +
                         interesting_all_cw2 + attention + valence_all_cw2 + (1|PID),
                       data = esm_inner, REML = FALSE)
b_main_1b <- summary(model1b_fit)$coefficients[,"Estimate"]
vc_1b <- as.data.frame(VarCorr(model1b_fit))
sd_u0_1b <- vc_1b[1, "sdcor"]
sd_e_1b  <- vc_1b[2, "sdcor"]

## Moderation model ----------------------------------------
## code extracted from original script

fit.totaleffect1 <- lme4::lmer(wb_cw2 ~ attention + (1|PID)
                               + activity+ clarity_all_cw2 + interesting_all_cw2 #covariate
                             , data=esm_inner)
summary(fit.totaleffect1) #estimate of attention = total effect; IV to DV line; .19, p<.001

Anova(fit.totaleffect1, type=3)

fit.mediator1<- lme4::lmer(valence_all_cw2 ~ attention + (1|PID)
                           +activity+ clarity_all_cw2+ interesting_all_cw2   #covariate
                         , data=esm_inner)
summary(fit.mediator1) #effect of IV onto mediator; IV to Mediator line; .12, p<.001
Anova(fit.mediator1, type=3)

fit.dv1 <- lme4::lmer(wb_cw2 ~ attention + valence_all_cw2 + (1|PID)
                      + activity+ clarity_all_cw2 + interesting_all_cw2 #covariate
                    , data=esm_inner)
summary(fit.dv1) #effect of mediator on DV; Mediator to DV line (ignore IV here); .60
                                        #since IV is still sig here, we say it's a partial meditation as opposed to a complete mediation
Anova(fit.dv1, type=3)


results1<-mediation::mediate(fit.mediator1, fit.dv1, treat='attention',
                            mediator='valence_all_cw2',
              covariates =  c("activity", "clarity_all_cw2", "interesting_all_cw2"),  #NEW
                            na.action="na.omit")
summary(results1)

## Data inspection

esm_sub <- esm_inner |>
  dplyr::select(wb, activity, clarity_all_cw2,
                interesting_all_cw2, attention, valence_all_cw2)

GGally::ggpairs(esm_sub)

## 3 Run the simulation --------------------------------------

## Define simulation parameters
N_to_test <- c(50, 75, 100, 125, 150, 175, 200)
n_sims <- 300
set.seed(1234)

plan(multisession, workers = availableCores() - 1)

## Model 1a ----
sim1a <- furrr::future_map_dfr(N_to_test, ~{
  n = .x
  purrr::map(1:n_sims, ~run_sim_main(
                         N = n,
                         esm_inner_orig = esm_inner,
                         prompt_counts_orig = prompt_counts,
                         beta_attention_ratio = 1,
                         use_variable_prompts = TRUE,
                         attrition_rate = 0.0,
                         model_type = "model1a"
  )) |>  purrr::list_rbind()
}, .options = furrr_options(seed = TRUE))

## Model 1b -----
sim1b <- furrr::future_map_dfr(N_to_test, ~{
  n = .x
  purrr::map(1:n_sims, ~run_sim_main(
                         N = n,
                         model_type = "model1b",
                         esm_inner_orig = esm_inner,
                         prompt_counts_orig = prompt_counts,
                         beta_attention_ratio = 0.75,
                         use_variable_prompts = TRUE,
                         attrition_rate = 0.3
  )) |>  purrr::list_rbind()
}, .options = furrr_options(seed = TRUE))

## 4 Calculate and Plot Power Curve ------------------------

power_1a <- sim1a |>
  dplyr::group_by(N) |>
  dplyr::summarise(power = mean(p.value < 0.05, na.rm = TRUE)) |>
  dplyr::mutate(Scenario = "Model 1a (Attention Only)")

power_1b <- sim1b |>
  dplyr::group_by(N) |>
  dplyr::summarise(power = mean(p.value < 0.05, na.rm = TRUE)) |>
  dplyr::mutate(Scenario = "Model 1b (Attention, given Valence)")


all_power_curves <- dplyr::bind_rows(power_1a, power_1b)

# Plot
ggplot(all_power_curves, aes(x = N, y = power, color = Scenario, group = Scenario)) +
  geom_line(linewidth = 1) +
  geom_point(size = 2.5) +
  geom_hline(yintercept = 0.80, linetype = "dashed", color = "red") +
  scale_y_continuous(labels = scales::percent, limits = c(0, 1)) +
  labs(
    title = "Power Analysis for Main Attention Effects",
    x = "Number of Participants Recruited (N)",
    y = "Power (p < .05)",
    color = "Analysis Scenario"
  ) +
  theme_minimal()
