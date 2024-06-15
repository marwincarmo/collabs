
# 1 Load packages and data ------------------------------------------------



# 2 Estimate var for each subject -----------------------------------------

# Get the unique ids
unique_ids <- unique(pc_data$subject_id)[1:3]

# Create an empty list to store results
results_list <- list()

# Create an empty vector to store ids that caused errors
error_ids <- c()

# Loop through each unique id

n <- length(unique_ids)
  
for (id in unique_ids) {
  
  ii <- which(unique_ids == id)
  extra <- nchar('||100%')
  width <- options()$width
  step <- round(ii / n * (width - extra))
  text <- sprintf('|%s%s|% 3s%%', strrep('=', step),
                  strrep(' ', width - step - extra), round(ii / n * 100))
  cat(text)
  
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
    rm(sub_df)
    gc() # Run garbage collection
  }, error = function(e) {
      error_ids <<- c(error_ids, id)
    })
  
  cat(if (ii == n) '\n' else '\014')
  }

# Save the results list to an .rds file
# file_name <- paste0("out/results_noise/id_", names(results_list[i]), ".rds")
# saveRDS(results_list[i], file_name)
# saveRDS(results_list, file = "results.rds")

# Load the results
#results_noise <- readRDS("out/results_noise.rds")

# Iterate over the list and save each element as a separate .rds file
# for (i in seq_along(results_noise)) {
#   # Define the filename for each element
#   file_name <- paste0("out/results_noise/id_", names(results_noise[i]), ".rds")
#   # Save the element
#   saveRDS(results_noise[i], file_name)
#   }


## IDs that failed to converge ----

converged <- as.numeric(gsub("[^0-9]", "", list.files("out/results_noise")))

failed <- unique(raw_data$subject_id)[!(unique(raw_data$subject_id) %in% converged)]

converged[!converged %in% unique_ids]

## Re-run with failed ids

error_ids <- c()
estimate_var_ggm(failed)

# The convergence failure was due to some participants having a small number 
# of observations. Let's check how many observations each of these participants have.

raw_data |> 
  dplyr::filter(subject_id %in% failed,
                concussion == 1) |> 
  dplyr::with_groups(subject_id, dplyr::summarise, n = dplyr::n()) |> 
  dplyr::arrange(desc(n)) |> 
  print(n = Inf)

# A few had a many observations. However, these have very little to no variance in 
# a lot of symptoms.


# 3 Estimate bridge strength ----------------------------------------------
source("scripts/roll_your_own.R")
library(networktools)
bridgestrenght <- function(x, ...){
  bridge(x, ...)$`Bridge Strength`
}

for (j in 1:length(results_list)) {
  
  res <- roll_your_own2(results_list[[j]], iter = 2000, select = TRUE, 
                        FUN = bridgestrenght, communities = 1:22, progress = TRUE, 
                        cores = 8)
  saveRDS(res, paste0("out/results_noise/bs/bs_", names(results_list[j]), ".rds") )
}


