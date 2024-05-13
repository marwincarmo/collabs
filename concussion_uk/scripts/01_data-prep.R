
# 00 Load packages --------------------------------------------------------

library(dplyr)

# 01 Load in data ---------------------------------------------------------

raw_data <- haven::read_sav("data/PCSNetworkModelAnalysis.sav")

dplyr::glimpse(raw_data)


# 02 Create between persons data ------------------------------------------

# Obtain participant means of each symptom averaging across all survey times.

raw_data |> 
  group_by(SubjectID) |> 
  arrange(SubjectID)

count(raw_data, SurveyDay) |> print(n=Inf)

count(raw_data, SurveyTime) |> print(n=Inf)
