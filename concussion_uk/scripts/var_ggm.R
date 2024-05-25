
# 1 Load packages and data ------------------------------------------------

library(parallel)
library(BGGM)

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

## Transform data to long ----

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

# 2 Estimate var for each subject -----------------------------------------

# Get the unique ids
unique_ids <- unique(pc_data$subject_id)

# Create an empty list to store results
results_list <- list()

# Create an empty vector to store ids that caused errors
error_ids <- c()
  
# Loop through each unique id
for (id in unique_ids) {
  
  tryCatch({
    # Subset the data frame for the current id
    sub_df <- subset(pc_data, subject_id == id)
  
    ## Check for constant values in variables:
    if( sum(apply(sub_df, 2, FUN = function(x) sd(x)) == 0) > 0 ) {
      pos <- which(apply(sub_df, 2, FUN = function(x) sd(x)) == 0)
      zero_var <- names(sub_df)[pos]
      keep <- setdiff(names(sub_df), zero_var)
      sub_df <- as.data.frame(sub_df[, keep])
    }
  
    # Compute the function
    result <- var_estimate(sub_df, progress = FALSE, iter = 5000)
  
    # Store the result in the list with id as the name
    results_list[[as.character(id)]] <- result
  
    # Remove the intermediate object to free up memory
    #rm(sub_df)
    #gc() # Run garbage collection
  }, error = function(e) {
      error_ids <<- c(error_ids, id)
    })
  }

# Save the results list to an .rds file
saveRDS(results_list, file = "results.rds")