
# 00 Load packages --------------------------------------------------------

library(tidyr)
library(dplyr)
library(BGGM)
library(micemd)

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
  dplyr::group_by(subject_id) %>%
  tidyr::complete(survey_day = all_combinations$survey_day, survey_time = all_combinations$survey_time) %>%
  dplyr::ungroup()


# 04 Check missing percentage ---------------------------------------------

# by survey time
df_complete %>% 
  dplyr::group_by(survey_time) %>% 
  dplyr::summarise(pct_missing = sum(is.na(head))/n())

# by survey day
df_complete %>% 
  dplyr::group_by(survey_day) %>% 
  dplyr::summarise(pct_missing = sum(is.na(head))/n()) %>% 
  print(n = Inf)

# by participant
df_complete %>% 
  dplyr::group_by(subject_id) %>% 
  dplyr::summarise(pct_missing = sum(is.na(head))/n()) |> 
  dplyr::arrange(desc(pct_missing))

## Check if there are missing values for group and sex ----

baseline <- raw_data[raw_data$survey_day == 0, ]

table(baseline$sex, useNA = "always")

table(baseline$concussion, useNA = "always")


# 05 Imputation -----------------------------------------------------------

data_for_imp <- df_complete |> 
  dplyr::mutate(measurement = dplyr::case_when(
    survey_time == 10 ~ paste0(survey_day, "a"),
    survey_time == 3 ~ paste0(survey_day, "b"),
    survey_time == 8 ~ paste0(survey_day, "c")
  )) |> 
  dplyr::select(subject_id, sex, concussion, measurement, dplyr::all_of(symptoms)) |> 
  dplyr::mutate(sex = factor(sex),
                concussion = factor(concussion)) |> 
  tidyr::fill(sex, .direction = "down") |> 
  tidyr::fill(concussion, .direction = "down")

# We want to impute data for the symptoms, that are count variables
# so we use the Poisson distribution

impMethod <- character(ncol(data_for_imp[, symptoms]))
names(impMethod) <- colnames(data_for_imp[, symptoms])

# set up predictor matrix
predMatrix <- matrix(0, ncol(data_for_imp), ncol(data_for_imp)) # create empty predictor
rownames(predMatrix) <- colnames(predMatrix) <- colnames(data_for_imp) # matrix

### Define predictors for each variable

predictors <-  c("subject_id", "sex", "concussion", "measurement", "head", 
                 "nau", "vom", "bal","diz","fat","flslp","slmo","slls",
                 "drw","lsen","nsen","irr","sad","nerv",      
                 "emo","num","slw","fog","attn","mem", "vis")

for (i in symptoms) {
  predMatrix[ i , predictors ] <- c(-2,1,1,3, rep(3, 22))
}

## Multilevel imputation ----

imputed <- mice.par(data_for_imp, m=30, imputationMethod="2l.glm.pois",
                    predictorMatrix=predMatrix, nnodes = 16) 

mice::complete(imputed, "long")

saveRDS(imputed, "out/PCSimputed.rds")

### Centering variables in imputed data set

PCSimputed_long <- mice::complete(imputed, "long", include=TRUE)

#### Standardizing symptoms across measurements

PCSimputed_z <- PCSimputed_long |> 
  dplyr::with_groups(
    .imp,
    dplyr::mutate,
    dplyr::across(
      dplyr::all_of(symptoms),
      ~as.vector(scale(.x)) 
    )
  ) |> 
  dplyr::with_groups(
    c(.imp, subject_id), # cluster (subject) mean
    dplyr::mutate,
    dplyr::across(
      dplyr::all_of(symptoms),
      ~mean(.x, na.rm=TRUE), .names = "{col}_mean") 
    ) 
# creating person centered variables

for (i in symptoms) {
  new_name <- paste0(i, "_c")
  mean_name <- paste0(i, "_mean")
  PCSimputed_z[new_name] <- PCSimputed_z[i] - PCSimputed_z[mean_name]
}

# Turning back into mids file.
## Entire sample
PCSimputed_z_mids <- mice::as.mids(
  PCSimputed_z[, !(names(PCSimputed_z) %in% c("subject_id", "concussion", "measurement"))], 
                                   .imp = ".imp", .id = ".id")
saveRDS(PCSimputed_z_mids, "out/imputed_datasets/PCSimputed_z_mids.rds")

## Concussion
PCSimputed_z_mids1 <- mice::as.mids(
  PCSimputed_z[PCSimputed_z$concussion==1, !(names(PCSimputed_z) %in% c("subject_id", "concussion", "measurement"))], 
  .imp = ".imp", .id = ".id")
saveRDS(PCSimputed_z_mids1, "out/imputed_datasets/PCSimputed_z_mids1.rds")

## Control
PCSimputed_z_mids0 <- mice::as.mids(
  PCSimputed_z[PCSimputed_z$concussion==0, !(names(PCSimputed_z) %in% c("subject_id", "concussion", "measurement"))], 
  .imp = ".imp", .id = ".id")
saveRDS(PCSimputed_z_mids0, "out/imputed_datasets/PCSimputed_z_mids0.rds")

## Create WITHIN mids file
### Entire sample
PCSimputed_z_within <- PCSimputed_z |> 
  dplyr::select(.imp, .id, sex, concussion,
                dplyr::ends_with("_c"))
PCSimputed_z_within_mids <- mice::as.mids(
  PCSimputed_z_within[, -4], .imp = ".imp", .id = ".id")
saveRDS(PCSimputed_z_within_mids, "out/imputed_datasets/PCSimputed_z_within_mids.rds")

### Concussion
PCSimputed_z_within_mids1 <- mice::as.mids(
  PCSimputed_z_within[PCSimputed_z_within$concussion==1, -4], 
  .imp = ".imp", .id = ".id")
saveRDS(PCSimputed_z_within_mids1, "out/imputed_datasets/PCSimputed_z_within_mids1.rds")

### Control
PCSimputed_z_within_mids0 <- mice::as.mids(
  PCSimputed_z_within[PCSimputed_z_within$concussion==0, -4], 
  .imp = ".imp", .id = ".id")
saveRDS(PCSimputed_z_within_mids0, "out/imputed_datasets/PCSimputed_z_within_mids0.rds")

## Create BETWEEN mids file

PCSimputed_z_between <- PCSimputed_z |> 
  dplyr::select(.imp, .id, sex, concussion,
                dplyr::ends_with("_mean"))

### Entire sample
PCSimputed_z_between_mids <- mice::as.mids(PCSimputed_z_between[, -4], 
                                           .imp = ".imp", .id = ".id")
saveRDS(PCSimputed_z_between_mids, "out/imputed_datasets/PCSimputed_z_between_mids.rds")

### Concussion
PCSimputed_z_between_mids1 <- mice::as.mids(
  PCSimputed_z_between[PCSimputed_z_between$concussion==1 , -4], 
                                           .imp = ".imp", .id = ".id")
saveRDS(PCSimputed_z_between_mids1, "out/imputed_datasets/PCSimputed_z_between_mids1.rds")

### Control
PCSimputed_z_between_mids0 <- mice::as.mids(
  PCSimputed_z_between[PCSimputed_z_between$concussion==0 , -4], 
  .imp = ".imp", .id = ".id")
saveRDS(PCSimputed_z_between_mids0, "out/imputed_datasets/PCSimputed_z_between_mids0.rds")


# 06 Create a dataset with NAs to impute with BGGM ------------------------

PCS_z <- df_complete |> 
  #dplyr::group_by(survey_day, survey_time) |> 
  dplyr::mutate(
    dplyr::across(
      dplyr::all_of(symptoms),
      ~as.vector(scale(.x)) 
    )
  ) |>
  #dplyr::ungroup() |> 
  dplyr::with_groups(
    subject_id, # cluster (subject) mean
    dplyr::mutate,
    dplyr::across(
      dplyr::all_of(symptoms),
      ~mean(.x, na.rm=TRUE), .names = "{col}_mean") 
  ) 
# creating person centered variables

for (i in symptoms) {
  new_name <- paste0(i, "_c")
  mean_name <- paste0(i, "_mean")
  PCS_z[new_name] <- PCS_z[i] - PCS_z[mean_name]
}

saveRDS(PCS_z, "out/PCS_z.rds")
