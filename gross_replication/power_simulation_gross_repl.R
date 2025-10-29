# Gross et al 2024 Replication Power Simulation
# Julia Bottesini
# June 10, 2025

library(lme4)
library(lmerTest)
library(tidyverse)

# Statistical models -----
# Hypothesis 1: Wellbeing will be higher for 
# From Gross et al 2024 prereg: First, using data from moments in which participants report having experienced “inner speech”, 
# we predict higher Well-Being for moments when people report being Present versus Not-Present 
# (i.e., a main effect of Attention State), in line with the results of Killingsworth & Gilbert, 2010. 
# For this analysis we use a multilevel model where the unit of analysis is “moment”,
# Attention State is the predictor variable, Well-Being is the dependent variable, 
# and Participant ID is entered as a random effect. 

# NOTE: prereg says nothing about covariates, but analysis in paper includes covariates

# y = B0 + B1*AttentionState
# B0 = GM + PM * Participant


# Analysis of existing data from Gross et al. 2024 ----
# setwd("")
g24_names <- read_csv("gross_replication/Corpus 6.csv") %>% names()
g24_full <- read_csv("gross_replication/Corpus 6.csv", skip = 4, col_names = g24_names)

g24_full <- g24_full %>% 
  dplyr::rename(
    duration = ...4,
    finished = ...5,
    participantID = ...6,
    participantNUM = Survey,
    attention = Q2,
    wb = Q3,
    inner_speech = Q4,
    clarity = Q5, 
    valence = Q6,
    reaction = Q7,
    interesting = Q8,
    inner_experience = Q9...18,
    no_inner_speech_clarity = Q10,
    no_inner_speech_valence = Q11,
    no_inner_speech_reaction = Q12,
    no_inner_speech_interesting = Q13,
    activity = Q14...24) %>% 
  distinct()

g24_full <- g24_full %>% 
  mutate(inner_speech = if_else(inner_speech == 1, "Yes", "No"),
         attention = if_else(attention == 1, "Present", "Not present"))

g24 <- g24_full %>% 
  select(participantID, attention, wb, inner_speech, valence) %>% 
  rename(wellbeing = wb)

g24 <- g24 %>% 
  filter(inner_speech == 'Yes')

# Hypothesis 1
m1 <- lmer(wellbeing ~ attention + (1|participantID),
           data = g24)
summary(m1)

# Hypothesis 2 (? unclear if this is right)
m2 <- lmer(wellbeing ~ attention + valence + (1|participantID),
           data = g24)
summary(m2)

## Data Ns ----
g24_full %>% 
  filter(inner_speech == "Yes" | inner_speech == "No") %>%
  distinct(...3, participantID, wb, inner_speech, attention, valence) %>% # remove repeated observations
  filter(inner_speech == "Yes") %>% 
  group_by(participantID) %>% 
  summarize(n_obs = n()) %>% 
  ungroup() %>% 
  summarize(n_participants = n(),
            n_obs = sum(n_obs))

# n_participants = 562 
# n_obs = 6059

# Data missingness distribution ----

# Get the distribution of number of prompts answered that contained inner speech out of 42 
# this is the actual distribution of "useful" data
# but it makes lots of assumptions about our data being the same as
# theirs. It's the best I can do, but far from ideal.
g24_full %>% 
  filter(inner_speech == "Yes" | inner_speech == "No") %>%
  distinct(...3, participantID, wb, inner_speech, attention, valence) %>% # remove repeated observations
  group_by(participantID, inner_speech) %>% 
  count() %>% 
  ungroup(inner_speech) %>% 
  mutate(participant_n = sum(n)) %>%
  pivot_wider(id_cols = c(participantID, participant_n), 
              names_from = inner_speech, values_from = n, names_prefix = "IS_") %>% 
  replace_na(replace = list(IS_Yes = 0, IS_No = 0)) %>% # add 0s where there are NAs (no observations of that kind)
  ungroup() %>% 
  mutate(prop_responded = participant_n/42,
         prop_inner_speech = IS_Yes/42) -> prob_data

hist(prob_data$IS_Yes)

# Attention state distribution ----
g24_full %>% 
  filter(inner_speech == "Yes" | inner_speech == "No") %>%
  distinct(...3, participantID, wb, inner_speech, attention, valence) %>% # remove repeated observations
  group_by(participantID) %>% 
  summarize(attention_prop = sum(attention == "Present")/n()) -> prob_attention
  
prob_attention %>% 
  ggplot(aes(x = attention_prop)) +
  geom_histogram()
  


# SIMULATION -----

## Parameters ----
# Number of participants
Ns = c(180, 200, 220, 240, 260) # x-axis values for the power plot

# Participant params
P_SD = 0.71 # standard deviation of participants (extracted from Gross et al. 2024 data, participant variance = 0.4955, sd = .7039, rounded up)
INNERSPEECH_DIST = prob_data$prop_inner_speech # Distribution of inner speech probability
ATTENTION_DIST = prob_attention$attention_prop # Distribution of probability of attention present

# Simulation params
GRAND_MEAN = 4 # overall average on DV, expected to be in the middle of scale
B1 = .2 # raw SESOI (smallest effect size of interest) for main effect, 1/5th of a point on the 7pt scale (gross et al raw effect = .46)
ERR_SD = 1.21 # std deviation of the noise (mean = 0) (extracted from Gross et al 2024 data, residual variance 1.4537, sd = 1.206, rounded up)
ITER = 10000 # number of iterations for this simulation


# y = B0 + B1*AttentionState
# B0 = GM + PM * Participant

### Participant level ----
one_participant <- function(){ #function to simulate data for one person
  
  one_part = expand_grid(Day = seq(1, 7),
                         Time = seq(1, 6))
  
  
  one_part = one_part %>% # add participant parameters
    mutate(GM = GRAND_MEAN, # grand mean wellbeing
           PM = rnorm(n = 1, sd = P_SD), # participant intercept wellbeing
           ERROR = rnorm(n = 42, sd = ERR_SD)) # error / noise
  
  attention_state_prob = sample(x = ATTENTION_DIST, size = 1, replace = T) # sample one attention state proportion from real data 
  attention_state_prob = attention_state_prob + rnorm(n = 1, sd = .1) # add a tiny bit of noise
  if(attention_state_prob > 1) {attention_state_prob = 1} # cut off if > 1
  if(attention_state_prob < 0) {attention_state_prob = 0} # cut off if negative (< 0)
  
  attention_state = sample(x = c(1, 0), # sample 42 attention states (Present = 1; Not present = 0)
                           size = 42, 
                           prob = c(attention_state_prob, 1-attention_state_prob), # with probability sampled above from data dist + noise
                           replace = T) 
  
  # add the attention variable
  one_part <- one_part %>% 
    mutate(attention = attention_state) #add attention variable
  
  
  # Sample number of obs that are present AND innerspeech
  inner_speech_prob = sample(x = INNERSPEECH_DIST, size = 1, replace = T) # sample one inner speech proportion from real data 
  inner_speech_prob = inner_speech_prob + rnorm(n = 1, sd = .1) # add a tiny bit of noise
  if(inner_speech_prob > 1) {inner_speech_prob = 1} # cut off if > 1
  if(inner_speech_prob < 0) {inner_speech_prob = 0} # cut off if negative (< 0)
  # hist(inner_speech_prob) # looks okay
  
  inner_speech_dat = sample(x = c("Yes", "No"), # sample 42 inner speech yes/no
                           size = 42, 
                           prob = c(inner_speech_prob, 1-inner_speech_prob), # with probability sampled above from data dist + noise
                           replace = T) 

  one_part = one_part %>% 
    mutate(inner_speech = inner_speech_dat)
  
  return(one_part)
}

### Study level ----
one_study <- function(sampleN) { # run one study with sample N
  
  # initialize dataframe for this study
  df_study = tibble(Day = numeric(),
                    Time = numeric(),
                    GM = numeric(),
                    PM = numeric(),
                    ERROR = numeric(),
                    attention = numeric())
  
  
  for (n in 1:sampleN) { # number of participants in this run
    new_p = one_participant() %>% 
      mutate(participantID = paste0("part",n))
    df_study = bind_rows(df_study, new_p)
  }
  
  # filter only inner speech == yes, calculate effective n
  df_study = df_study %>% 
    filter(inner_speech == "Yes") 
  
  df_study %>% 
    group_by(participantID) %>% 
    count() %>% 
    ungroup() %>% 
    summarize(part_n = n(),
              obs_n = sum(n)) -> samplesize
  
  # calculate wellbeing
  df_study = df_study %>% 
    mutate(wellbeing = GM + PM + B1*attention + ERROR)
  
  
  # fit the model
  m1 <- lmer(wellbeing ~ attention + (1|participantID),
             data = df_study)
  
  pval = summary(m1)$coef[2,5] # extract p-value for attention effect
  
  return(tibble(pvalue = pval, 
                part_n = samplesize$part_n, 
                obs_n = samplesize$obs_n))
  
}


## Check data distributions against real data ----

### Distributions of real data for comparison ----

# Sample size to use for simulated study
## This will result in about the same sample the data has
## after removing the 25% no-inner-speech people.
## Each time a study is run it will look a bit different
## but as long as it's close enough it's fine
sampleN = 750

# initialize dataframe for a study
df_study = tibble(Day = numeric(),
                  Time = numeric(),
                  GM = numeric(),
                  PM = numeric(),
                  ERROR = numeric(),
                  attention = numeric())


for (n in 1:sampleN) { # number of participants in this run
  new_p = one_participant() %>% 
    mutate(participantID = paste0("part",n))
  df_study = bind_rows(df_study, new_p)
}

# calculate wellbeing
df_study = df_study %>% 
  mutate(wellbeing = GM + PM + B1*attention + ERROR)

### check the data characteristics vs the simulated data ----

# distribution of inner speech observations (data)
g24_full %>% 
  filter(inner_speech == "Yes") %>%
  distinct(...3, participantID, wb, inner_speech, attention, valence) %>% # remove repeated observations
  group_by(participantID) %>% 
  count() %>% 
  ggplot(aes(x = n)) +
  geom_histogram()


#distribution of inner speech (simulation)
df_study %>% 
  filter(inner_speech == "Yes") %>% 
  group_by(participantID) %>%
  count() %>% 
  ggplot(aes(x = n)) +
  geom_histogram()


# proportion of attention variable (data)
g24_full %>% 
  filter(inner_speech == "Yes" | inner_speech == "No") %>%
  distinct(...3, participantID, wb, inner_speech, attention, valence) %>% # remove repeated observations
  group_by(attention) %>% 
  count() %>% 
  ungroup() %>% 
  mutate(prop = n/sum(n)) 

# proportion of attention variable (simulation)
df_study %>% 
  group_by(attention) %>% 
  count() %>%
  ungroup() %>% 
  mutate(prop = n/sum(n))


# distribution of wellbeing (data)
g24_full %>% 
  filter(inner_speech == "Yes" | inner_speech == "No") %>%
  distinct(...3, participantID, wb, inner_speech, attention, valence) %>% # remove repeated observations
  ggplot(aes(x = wb)) +
  geom_histogram()

# summary stats
g24_full %>% 
  filter(inner_speech == "Yes" | inner_speech == "No") %>%
  distinct(...3, participantID, wb, inner_speech, attention, valence) %>%
  summarize(M = mean(wb),
            Med = median(wb),
            SD = sd(wb))

# distribution of wellbeing (simulation)
# Note: obviously it's not discrete like the real data but it's pretty close
# in terms of its summary stats which is what we want
df_study %>% 
  ggplot(aes(x = wellbeing)) +
  geom_histogram()

# summary stats
df_study %>% 
  summarize(M = mean(wellbeing),
            Med = median(wellbeing),
            SD = sd(wellbeing))


# RUN SIM: Sample Size level ----

set.seed(123)

# tibble to save results
RESULTS <- tibble(TOTAL_N = numeric(), Avg_Pwr = numeric(), 
                  Avg_PartN = numeric(), SD_PartN = numeric(), 
                  Avg_ObsN = numeric(), SD_ObsN = numeric())

for (N in Ns) { # run for every N in Ns
  print(paste0("Testing sample N = ",N))
  
  outcomes = tibble(pvalue = numeric(),
                    part_n = numeric(),
                    obs_n = numeric())
  
  for (i in 1:ITER) { # run ITER studies of size N
    print(paste("Iteration #",i))
    outcomes = bind_rows(outcomes, one_study(N))
  }
  
  outcomes %>% # summarize outcome data
    mutate(detected = pvalue < 0.05) %>% 
    summarize(Avg_Pwr = mean(detected), # how often do we detect the effect with this N
              Avg_PartN = mean(part_n), # average effective N (participants)
              SD_PartN = sd(part_n),
              Avg_ObsN = mean(obs_n), # average effective N (observations)
              SD_ObsN = sd(obs_n)) %>% 
    mutate(TOTAL_N = N) -> outcomes

  
  RESULTS = bind_rows(RESULTS, outcomes)
  
}

RESULTS

# > RESULTS
# # A tibble: 5 × 6
#   TOTAL_N Avg_Pwr Avg_PartN SD_PartN Avg_ObsN SD_ObsN
#     <dbl>   <dbl>     <dbl>    <dbl>    <dbl>   <dbl>
# 1     180   0.798      154.     4.69    1941.    130.
# 2     200   0.833      172.     4.95    2160.    137.
# 3     220   0.869      189.     5.17    2375.    143.
# 4     240   0.894      206.     5.41    2588.    149.
# 5     260   0.917      223.     5.67    2806.    157.

# save results
# write_rds(RESULTS, "./power_sim_results.rds")

# read results back in (so you don't have to run this again)
RESULTS <- read_rds("./power_sim_results.rds")

# Power plot
RESULTS %>%
  ggplot(aes(x = TOTAL_N, y = Avg_Pwr*100)) +
  geom_point() +
  geom_label(aes(label = sprintf("%0.2f", round(Avg_Pwr*100, digits = 2)), hjust=.53, vjust=1.7)) +
  geom_line() +
  theme_classic() +
  ylim(60,100) +
  labs(x = "Number of participants (N)", 
       y = "Power to detect main effect for attention on wellbeing (%)")

