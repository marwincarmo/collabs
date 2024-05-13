
# 00 Load packages --------------------------------------------------------

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

# 02 Create between persons data ------------------------------------------

# Obtain participant means of each symptom averaging across all survey times.

raw_data |> 
  group_by(SubjectID) |> 
  arrange(SubjectID)

count(raw_data, SurveyDay) |> print(n=Inf)

count(raw_data, SurveyTime) |> print(n=Inf)
