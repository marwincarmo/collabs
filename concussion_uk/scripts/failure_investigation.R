
# Finding observations that failed to converge ----------------------------

raw_data <- haven::read_sav("data/PCSNetworkModelAnalysis.sav") |> 
  janitor::clean_names()
data_dictionary <- readxl::read_excel("Data Dictionary.xlsx")
symptoms <- names(raw_data)[19:length(names(raw_data))]
symptoms_names <- data_dictionary[data_dictionary$`Variable Name` %in% symptoms,]$`Variable Label`
long_data <- raw_data |> 
  tidyr::pivot_longer(
    cols = dplyr::all_of(symptoms),
    names_to = "symptom",
    values_to = "score"
  )

## Transform data to long 
# person centering
c_data <- long_data |> 
  dplyr::with_groups(
    symptom,
    dplyr::mutate,
    z_score = as.vector(scale(score)) # scale the symptoms across all individuals and all time points
  ) |> 
  dplyr::with_groups(
    c(symptom, subject_id),
    dplyr::mutate,
    mean_score = mean(z_score, na.rm = TRUE) # mean person score in standardized units
  ) |> 
  dplyr::with_groups(
    symptom,
    dplyr::mutate,
    pc_score = z_score - mean_score # person centered score is the person deviation 
    # at each measurement from their own average 
  )
# revert back to wide
pc_data <- c_data[,c("subject_id", "pc_score", "symptom")] |> 
  tidyr::pivot_wider(id_cols = subject_id,
                     names_from = symptom,
                     values_from = pc_score) |>
  tidyr::unnest()

converged <- as.numeric(gsub("[^0-9]", "", list.files("out/results_noise")))

failed <- unique(raw_data$subject_id)[!(unique(raw_data$subject_id) %in% converged)]


# Looking at each case ----------------------------------------------------

pc_data |> 
  dplyr::filter(subject_id %in% failed)

## First lets rule out how many had less data than nodes

n_nodes <- ncol(pc_data) -1 #excluding subject_id

low_n <- pc_data |> 
  dplyr::filter(subject_id %in% failed
  ) |> 
  dplyr::with_groups(subject_id, dplyr::summarise, n = dplyr::n()) |> 
  dplyr::filter(n < n_nodes) |> 
  dplyr::pull(subject_id)
# There were 58/113

## Now looking at the other cases

pc_data |> 
  dplyr::filter(subject_id %in% failed,
                !(subject_id %in% low_n)) |> 
  dplyr::with_groups(subject_id, dplyr::summarise, n = dplyr::n()) |> 
  dplyr::arrange(desc(n)) |> 
  print(n = Inf)

library(dplyr)

# Compute standard deviation for each symptom per subject (wide format)
sd_wide_df <- pc_data |> 
  dplyr::filter(subject_id %in% failed,
                !(subject_id %in% low_n)) |> 
  dplyr::group_by(subject_id) |> 
  dplyr::summarise(dplyr::across(
    dplyr::everything(), sd, na.rm = TRUE, .names = "sd_{col}"))

# View the result
sd_wide_df


pc_data |> 
  dplyr::filter(subject_id %in% converged) |> 
  dplyr::group_by(subject_id) |> 
  dplyr::summarise(dplyr::across(
    dplyr::everything(), sd, na.rm = TRUE, .names = "sd_{col}"))
