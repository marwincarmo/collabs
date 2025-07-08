
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

raw_data |> 
  dplyr::filter(subject_id %in% failed) |> 
  dplyr::with_groups(concussion, dplyr::summarise, n = dplyr::n()) |> 
  dplyr::arrange(desc(n)) 

## First lets rule out how many had less data than nodes

n_nodes <- ncol(pc_data) -1 #excluding subject_id

low_n <- pc_data |> 
  dplyr::filter(subject_id %in% failed
  ) |> 
  dplyr::with_groups(subject_id, dplyr::summarise, n = dplyr::n()) |> 
  dplyr::filter(n < n_nodes) |> 
  dplyr::pull(subject_id)
# There were 58/113

### Testing
library(BGGM)
sub_df <- subset(pc_data, subject_id == 45)
result <- var_estimate(sub_df[,-1], progress = TRUE, iter = 5000)

## Now looking at the other cases

multicol <- pc_data |> 
  dplyr::filter(subject_id %in% failed,
                !(subject_id %in% low_n)) 

multicol |> 
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
    dplyr::everything(), ~sd(., na.rm = TRUE), .names = "sd_{col}"))

# View the result
sd_wide_df


pc_data |> 
  dplyr::filter(subject_id %in% converged) |> 
  dplyr::group_by(subject_id) |> 
  dplyr::summarise(dplyr::across(
    dplyr::everything(), sd, na.rm = TRUE, .names = "sd_{col}"))


# Eric's questions --------------------------------------------------------

## 1. Those who didn't complete enough responses (n=58)
### i. What are the subsample sizes for the control and concussion groups?

raw_data |> 
  dplyr::filter(subject_id %in% low_n) |> 
  dplyr::with_groups(subject_id,
                     dplyr::slice_head) |> 
  dplyr::count(concussion) |> 
  dplyr::mutate(pct = n/sum(n))

# concussion     n   pct
# <dbl+lbl>  <int> <dbl>
# 1 0 [No]        38 0.655
# 2 1 [Yes]       20 0.345

### ii. For this group overall, what was the Mean, Standard Deviation, and range of # of responses completed

raw_data |> 
  dplyr::filter(subject_id %in% low_n) |> 
  dplyr::with_groups(subject_id,
                     dplyr::summarise, n = dplyr::n()) |> 
  dplyr::pull(n) |> 
  psych::describe()

#    vars  n mean  sd median trimmed  mad min max range skew kurtosis   se
# X1    1 58 8.59 6.4      7    8.15 7.41   1  21    20 0.51    -1.08 0.84

## 2. Those with multicollinearity (n=55)
### i. What are the subsample sizes for the control and concussion groups?

raw_data |> 
  dplyr::filter(subject_id %in% unique(multicol$subject_id)) |> 
  dplyr::with_groups(subject_id,
                     dplyr::slice_head) |> 
  dplyr::count(concussion) |> 
  dplyr::mutate(pct = n/sum(n))

# concussion     n   pct
# <dbl+lbl>  <int> <dbl>
# 1 0 [No]        42 0.764
# 2 1 [Yes]       13 0.236

### ii. Table summarizing the most common pairs of highly correlated symptoms and the magnitude of correlations
### Example: "For instance, for n=17 cases, vomiting was highly correlated with Sensitivity to noise and light (r>.xx)"

# add random noise
set.seed(999)
for (j in 2:ncol(pc_data)) {
  for (i in 1:nrow(pc_data)){
    noise <- rnorm(1, 0, 0.001)
    pc_data[i, j] <- pc_data[i, j] + noise 
  }
}


# Step 1: Compute personal correlation matrices for each subject
personal_cor_matrices <- pc_data |> 
  dplyr::filter(subject_id %in% 
                  c(unique(multicol$subject_id))) |> 
  dplyr::group_by(subject_id) |> 
  tidyr::nest() |> 
  dplyr::mutate(
    cor_matrix = purrr::map(
      data,
      ~cor(.[, -1], use = "complete.obs")
    )
  )

# Step 2: Average correlations across subjects

# Extract all correlation matrices into a list
cor_list <- personal_cor_matrices$cor_matrix

# Initialize a matrix to store the sum of correlations
sum_cor_matrix <- matrix(0, nrow = ncol(cor_list[[1]]), ncol = ncol(cor_list[[1]]))
colnames(sum_cor_matrix) <- colnames(cor_list[[1]])
rownames(sum_cor_matrix) <- colnames(cor_list[[1]])

# Sum all correlation matrices
for (cor_matrix in cor_list) {
  sum_cor_matrix <- sum_cor_matrix + cor_matrix
}

# Compute the average correlation matrix
average_cor_matrix <- sum_cor_matrix / length(cor_list)

# Step 3: Count how frequently each pair surpasses a threshold
threshold <- 0.7  # Set your correlation threshold here

# Initialize a matrix to store the count of high correlations
high_cor_count_matrix <- matrix(0, nrow = ncol(cor_list[[1]]), ncol = ncol(cor_list[[1]]))
colnames(high_cor_count_matrix) <- colnames(cor_list[[1]])
rownames(high_cor_count_matrix) <- colnames(cor_list[[1]])

# Count how often each pair exceeds the threshold
for (cor_matrix in cor_list) {
  # Create a logical matrix indicating whether correlations exceed the threshold
  # Use `is.na()` to ignore NA values
  exceeds_threshold <- cor_matrix > threshold
  exceeds_threshold[is.na(exceeds_threshold)] <- FALSE  # Treat NA as FALSE
  
  # Add to the count matrix
  high_cor_count_matrix <- high_cor_count_matrix + exceeds_threshold
}

# Convert the count matrix to a data frame for easier interpretation
high_cor_count_df <- as.data.frame(high_cor_count_matrix) |> 
  tibble::rownames_to_column(var = "Symptom1") |> 
  tidyr::pivot_longer(-Symptom1, names_to = "Symptom2", values_to = "Frequency") |> 
  dplyr::filter(Symptom1 < Symptom2) |> # Avoid duplicate pairs (e.g., A-B and B-A)
  dplyr::arrange(dplyr::desc(Frequency))
high_cor_count_df


## For the  sample that converged
# A tibble: 210 × 3
# Symptom1 Symptom2 Frequency
# <chr>    <chr>        <dbl>
# 1 emo      sad             57
# 2 flslp    slls            48
# 3 fog      slw             45
# 4 attn     fog             44
# 5 lsen     nsen            39
# 6 drw      fat             30
# 7 emo      nerv            30
# 8 attn     mem             28
# 9 fog      mem             27
# 10 attn     slw             26

## How many responses were obtained from control and concussion groups

raw_data |> 
  dplyr::filter(subject_id %in% low_n) |> 
  #dplyr::filter(subject_id %in% failed) |> 
  dplyr::with_groups(c(concussion, subject_id), 
                     dplyr::count) |> 
  #psych::describe()
  psych::describeBy( group = "concussion")
