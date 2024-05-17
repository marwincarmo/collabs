
# 00 Load packages --------------------------------------------------------

library(tidyr)
library(dplyr)

# 01 Load in data ---------------------------------------------------------

raw_data <- haven::read_sav("data/PCSNetworkModelAnalysis.sav") |> 
  janitor::clean_names()

dplyr::glimpse(raw_data)


# 02 Check for missing data -----------------------------------------------

symptoms <- names(raw_data)[19:length(names(raw_data))]

long_data <- raw_data |> 
  dplyr::select(c(subject_id, survey_time, survey_day, head:vis)) |> 
  tidyr::pivot_longer(cols =head:vis,
                      names_to = "symptom",
                      values_to = "score")
long_data |> 
  dplyr::with_groups(symptom,
                     dplyr::summarize,
                       n_missing = sum(is.na(score))
                     ) |> 
  print(n = Inf)

# are the missing observations already deleted?


# 03 Create a new dataset with NAs ----------------------------------------

data <- raw_data |> 
  dplyr::filter(survey_day != 0)

all_combinations <- expand.grid(
  survey_day = 1:30,
  survey_time = c(10, 3, 8)
)

# Ensure all combinations for each ID are present, filling missing values with NA
df_complete <- data %>%
  group_by(subject_id) %>%
  complete(survey_day = all_combinations$survey_day, survey_time = all_combinations$survey_time) %>%
  ungroup()

# Print the resulting dataframe
print(df_complete)


# 04 impute data with BGGM ------------------------------------------------

imp_data <- impute_data(df_complete[, symptoms])

saveRDS(df_complete, "data/PSCMissings.rds")
