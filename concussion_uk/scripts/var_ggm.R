
# 1 Load packages and data ------------------------------------------------

library(parallel)
library(BGGM)
library(network)
library(ggplot2)


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

# add random noise
set.seed(999)
for (j in 2:ncol(pc_data)) {
  for (i in 1:nrow(pc_data)){
    noise <- rnorm(1, 0, 0.001)
    pc_data[i, j] <- pc_data[i, j] + noise 
  }
}

# 2 Estimate var for each subject -----------------------------------------

# Get the unique ids
unique_ids <- unique(pc_data$subject_id)

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

## VAR for the whole sample

mle_net <- var_estimate(pc_data[-1], progress = TRUE, iter = 5000)

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

source("helpers/roll_your_own.R")
library(networktools)
bridgestrenght <- function(x, ...){
  bridge(x, ...)$`Bridge Strength`
}


converged <- na.omit(as.numeric(gsub("[^0-9]", "", list.files("out/results_noise"))))

for (j in 1:length(converged)) {
  
  result <- readRDS(paste0("out/results_noise/id_", converged[j], ".rds"))
  
  res <- roll_your_own2(result[[1]], iter = 5000, select = TRUE, 
                        FUN = bridgestrenght, communities = 1:22, progress = TRUE, 
                        cores = 16)
  saveRDS(res, paste0("out/results_noise/bs/bs_", names(result), ".rds") )
}

# Find the central symptom ----

bs_list <- lapply(paste0("out/results_noise/bs/", list.files("out/results_noise/bs")), readRDS)
id_vector <- as.numeric(gsub("[^0-9]", "", list.files("out/results_noise/bs")))
names(bs_list) <- id_vector
concussion_id <- unique(raw_data[raw_data$concussion == 1,]$subject_id)
control_id <- unique(raw_data[raw_data$concussion != 1,]$subject_id)


getCentral <- function(x, cred = 0.95) {
  
  lb <- (1-cred) / 2
  ub <- 1 - lb
  
  mu <-  apply( x$results, 1, mean)
  p <- length(mu)
  scale <- apply( x$results, 1, sd)
  ci_lb <- apply( x$results, 1, quantile, lb)
  ci_ub <- apply( x$results, 1, quantile, ub)
  
  res<- data.frame(Node = symptoms_names,
                   Post.mean = round(mu, 3),
                   Post.sd = round(scale, 3),
                   Cred.lb = round(ci_lb, 3),
                   Cred.ub = round(ci_ub, 3))
  rownames(res) <- NULL
  tail(res[order(res$Post.mean), ],1)
}

y <- purrr::map_dfr(bs_list, getCentral, .id = "subject_id")
y$concussion <- ifelse(y$subject_id %in% concussion_id, 1, 0)

table2 <- y |> 
  dplyr::with_groups(Node, 
                     dplyr::mutate,
                     freq_total = dplyr::n()) |> 
  dplyr::with_groups(c(concussion, Node), 
                     dplyr::mutate,
                     freq_group = dplyr::n()) |> 
  dplyr::with_groups(c(concussion, Node), 
                     dplyr::summarise,
                     freq_group = mean(freq_group),
                     freq_total = mean(freq_total)) |> 
  tidyr::pivot_wider(names_from = "concussion", values_from = "freq_group",
                     names_prefix = "freq_") |> 
  dplyr::mutate(pct_total = freq_total/sum(freq_total)*100, 
                pct_0 = freq_0/sum(freq_0)*100,
                pct_1 = freq_1/sum(freq_1, na.rm=TRUE)*100) |> 
  dplyr::arrange(factor(Node, levels = symptoms_names))

saveRDS(table2, "out/table2.rds")

# we might obtain table 3 estimating the within network as the subject difference
# from their own mean.


# 4 Average edge strength -------------------------------------------------

pb <- txtProgressBar(min = 0,      # Minimum value of the progress bar
                     max = length(converged), # Maximum value of the progress bar
                     style = 3,    # Progress bar style (also available style = 1 and style = 2)
                     width = 100,   # Progress bar width. Defaults to getOption("width")
                     char = "=")   # Character used to create the bar

z_list <- list()

for (j in 1:length(converged)) {
  
  result <- readRDS(paste0("out/results_noise/id_", converged[j], ".rds"))
  sel <- BGGM::select(result[[1]])
  pcor <- sel$pcor_weighted_adj
  z_list[[j]] <- BGGM::fisher_r_to_z(pcor)
  setTxtProgressBar(pb, j)
}

saveRDS(z_list, "out/z_list.rds")

z_list <- readRDS("out/z_list.rds")

sum_matrix <- Reduce("+", z_list)
mean_matrix <- sum_matrix / length(z_list)

avg_wgt <- BGGM::fisher_z_to_r(mean_matrix)

## Plot average within network ----

plotNet <- function(x) {
  
  p <- ncol(x)
  diag(x) <- 0
  net <- network::network(x)
  
  # edge weights
  network::set.edge.value(x = net, attrname = "weights",
                          value = x)
  
  # edge weights absolute
  network::set.edge.value(x = net, attrname = "abs_weights",
                          value = abs(x) * 5)
  # edge colors
  network::set.edge.attribute(x = net, attrname = "edge_color",
                              value = ifelse(net %e% "weights" < 0,
                                             "#D55E00",
                                             "#009E73"))
  e <- abs(as.numeric(x))
  
  plt <- GGally::ggnet2(net, edge.alpha = e[e != 0] / max(e),
                        edge.size = "abs_weights",
                        edge.color = "edge_color",
                        node.size = 13,
                        mode = "circle") +
    geom_point(color = "black",
               size = 14) +
    geom_point(size = 13, color = "white")  +
    guides(color = "none") +
    geom_text(label = symptoms) 
  plt
  
}


plotNet(avg_wgt) +
  labs(title = "Average within person network - Total sample")


# 4.1 Average edge strength for concussion group --------------------------

z_list_concussion <- list()
concussion <- converged[converged %in% unique(raw_data[raw_data$concussion==1, ]$subject_id)]

for (j in 1:length(concussion)) {
  
  result <- readRDS(paste0("out/results_noise/id_", concussion[j], ".rds"))
  sel <- BGGM::select(result[[1]])
  pcor <- sel$pcor_weighted_adj
  z_list_concussion[[j]] <- BGGM::fisher_r_to_z(pcor)
  setTxtProgressBar(pb, j)
}


saveRDS(z_list_concussion, "out/z_list_concussion.rds")

mean_matrix_concussion <- Reduce("+", z_list_concussion) / length(z_list_concussion)

avg_wgt_concussion <- BGGM::fisher_z_to_r(mean_matrix_concussion)

plotNet(avg_wgt_concussion) +
  labs(title = "Average within person network - Concussion")

# 4.2 Average edge strength for control group --------------------------

z_list_control <- list()
control <- converged[converged %in% unique(raw_data[raw_data$concussion==0, ]$subject_id)]

for (j in 1:length(control)) {
  
  result <- readRDS(paste0("out/results_noise/id_", control[j], ".rds"))
  sel <- BGGM::select(result[[1]])
  pcor <- sel$pcor_weighted_adj
  z_list_control[[j]] <- BGGM::fisher_r_to_z(pcor)
  setTxtProgressBar(pb, j)
}


saveRDS(z_list_control, "out/z_list_control.rds")

mean_matrix_control <- Reduce("+", z_list_control) / length(z_list_control)

avg_wgt_control <- BGGM::fisher_z_to_r(mean_matrix_control)

plotNet(avg_wgt_control) +
  labs(title = "Average within person network - Control")
