library(networktools)
library(BGGM)
library(parallel)

InAndOut <-  function(x, ...){
  res <- list(
    "InStrength" =  bridge(x, ...)$`Bridge Indegree`,
    "OutStrength" =  bridge(x, ...)$`Bridge Outdegree`
  )
   # bridge(x, ...)$`Bridge Strength`
   return(res)
}

# object <- readRDS("out/results_noise/id_138.rds")$`138`
# 
# sel <- BGGM::select(subject$`138`)
# 
# f(x = abs(sel$beta_weighted_adj), directed = TRUE, communities = 1:22)
# FUN = f; iter = NULL; select = TRUE; cred = 0.95; progress = TRUE; cores = 12

temporal_strength <- function(object,
                              FUN,
                              iter = NULL,
                              select = FALSE,
                              cred = 0.95,
                              progress = TRUE,
                              cores = 2,
                              ...) {
  
  sel <- BGGM::select(object, cred = cred)
  adj <- sel$beta_adj
  iter <- object$iter
  cn <- colnames(object$pcor_mu)
  
  # pcors are actually betas
  pcors <- object$fit$beta[,,51:(object$iter +50)]
  
  if(isTRUE(progress)){
    pb <- utils::txtProgressBar(min = 0, max = iter, style = 3)
  }
  
  getRes <- function(x) {
    
    pcors_s <- pcors[, , x] * adj
    
    est <- FUN(pcors_s, , directed = TRUE, communities = 1:22)#...)
    
    if(isTRUE(progress)){

      utils::setTxtProgressBar(pb, x)

    }
    
    est
    
  }
  
  results <- mclapply(X = 1:iter, FUN = getRes , mc.cores = cores)
  
  # Convert list of length iter (each element is a list of 2 vectors) into matrices
  InStrength <- do.call(cbind, lapply(results, `[[`, 1))  # Extract first element
  OutStrength <- do.call(cbind, lapply(results, `[[`, 2))  # Extract second element
  
  rownames(InStrength) <- cn
  rownames(OutStrength) <- cn
  
  returned_object <- list(InStrength= InStrength, OutStrength= OutStrength, iter = iter)
  
  return(returned_object)

}

print_roll_your_own <- function(x, cred = 0.95, ...) {
  cat("BGGM: Bayesian Gaussian Graphical Models \n")
  cat("--- \n")
  cat("Network Stats: Roll Your Own\n")
  cat("Posterior Samples:", x$iter, "\n")
  cat("--- \n")
  cat("Estimates: \n\n")
  
  lb <- (1 - cred) / 2
  ub <- 1 - lb
  
  process_results <- function(results) {
    dims <- dim(results)
    
    if (is.null(dims)) {
      res <- data.frame(
        Post.mean = round(mean(results), 3),
        Post.sd = round(sd(results), 3),
        Cred.lb = round(quantile(results, probs = lb), 3),
        Cred.ub = round(quantile(results, probs = ub), 3)
      )
    } else {
      mu <- apply(results, 1, mean)
      p <- length(mu)
      scale <- apply(results, 1, sd)
      ci_lb <- apply(results, 1, quantile, lb)
      ci_ub <- apply(results, 1, quantile, ub)
      
      res <- data.frame(
        Node = rownames(results),
        Post.mean = round(mu, 3),
        Post.sd = round(scale, 3),
        Cred.lb = round(ci_lb, 3),
        Cred.ub = round(ci_ub, 3)
      )
    }
    
    return(res)
  }
  
  cat("In-strength:\n")
  InStrength <- process_results(x$InStrength)
  print(InStrength, row.names = FALSE)
  
  cat("\nOut-strength:\n")
  OutStrength <- process_results(x$OutStrength)
  print(OutStrength, row.names = FALSE)
  
  cat("--- \n")
}


# Testing the function ----------------------------------------------------

# res <- temporal_strength(object, FUN = f, select = TRUE, cores = 12)
# pres <- print_roll_your_own(res)
# 
# pres
# 
# saveRDS(res, 'out/results_noise/in_and_out/ios_138.rds')

# Estimate in- and out-strength ----------------------------------------------


converged <- na.omit(as.numeric(gsub("[^0-9]", "", list.files("out/results_noise"))))

for (j in 98:length(converged)) {
  
  result <- readRDS(paste0("out/results_noise/id_", converged[j], ".rds"))
  
  cat("Estimate", j, 'out of', length(converged))
  
  res <- temporal_strength(result[[1]], iter = 5000, select = TRUE, 
                        FUN = InAndOut, communities = 1:22, progress = TRUE, 
                        cores = 12)
  
  saveRDS(res, paste0("out/results_noise/in_and_out/ios_", names(result), ".rds") )
}

# Find the central symptom ----

raw_data <- haven::read_sav("data/PCSNetworkModelAnalysis.sav") |> 
  janitor::clean_names()
data_dictionary <- readxl::read_excel("Data Dictionary.xlsx")
symptoms <- names(raw_data)[19:length(names(raw_data))]
symptoms_names <- data_dictionary[data_dictionary$`Variable Name` %in% symptoms,]$`Variable Label`

ios_list <- lapply(paste0("out/results_noise/in_and_out/", list.files("out/results_noise/in_and_out")), readRDS)
id_vector <- as.numeric(gsub("[^0-9]", "", list.files("out/results_noise/in_and_out")))
names(ios_list) <- id_vector
concussion_id <- unique(raw_data[raw_data$concussion == 1,]$subject_id)
control_id <- unique(raw_data[raw_data$concussion != 1,]$subject_id)

getCentralIn <- function(results, cred = 0.95) {
  
    lb <- (1 - cred) / 2
    ub <- 1 - lb
    
    mu <- apply(results$InStrength, 1, mean)
    p <- length(mu)
    scale <- apply(results$InStrength, 1, sd)
    ci_lb <- apply(results$InStrength, 1, quantile, lb)
    ci_ub <- apply(results$InStrength, 1, quantile, ub)
    
    res <- data.frame(
      Node = rownames(results$InStrength),
      Post.mean = round(mu, 3),
      Post.sd = round(scale, 3),
      Cred.lb = round(ci_lb, 3),
      Cred.ub = round(ci_ub, 3)
    )
    rownames(res) <- NULL
  return(tail(res[order(res$Post.mean), ],1))
}

getCentralOut <- function(results, cred = 0.95) {
  
  lb <- (1 - cred) / 2
  ub <- 1 - lb
  
  mu <- apply(results$OutStrength, 1, mean)
  p <- length(mu)
  scale <- apply(results$OutStrength, 1, sd)
  ci_lb <- apply(results$OutStrength, 1, quantile, lb)
  ci_ub <- apply(results$OutStrength, 1, quantile, ub)
  
  res <- data.frame(
    Node = rownames(results$OutStrength),
    Post.mean = round(mu, 3),
    Post.sd = round(scale, 3),
    Cred.lb = round(ci_lb, 3),
    Cred.ub = round(ci_ub, 3)
  )
  rownames(res) <- NULL
  return(tail(res[order(res$Post.mean), ],1))
}

## apply functions to all subjects
central_in <- purrr::map_dfr(ios_list, getCentralIn, .id = "subject_id")
central_out <- purrr::map_dfr(ios_list, getCentralOut, .id = "subject_id")

## create a column for identifying group
central_in$concussion <- ifelse(central_in$subject_id %in% concussion_id, 1, 0)
central_out$concussion <- ifelse(central_out$subject_id %in% concussion_id, 1, 0)

createCentralTable <- function(x) {
  x |> 
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
  
}

table_in <- createCentralTable(central_in)
table_out <- createCentralTable(central_out)

table_print <- function(x) {
  x |> 
    dplyr::transmute(
      Total = paste0(freq_total, " (", round(pct_total,2), ")"),
      Concussion = paste0(freq_1, " (", round(pct_1,2), ")"),
      Control = paste0(freq_0, " (", round(pct_0,2), ")"),
    )
  
} 

## Print to word

library(officer)
library(flextable)
library(magrittr)

# Create flextable object
ft <- flextable(data = table_print(table_out)) %>% 
  #theme_zebra %>% 
  autofit
# See flextable in RStudio viewer
ft

# Create a temp file
tmp <- tempfile(fileext = ".docx")

# Create a docx file
read_docx() %>% 
  body_add_flextable(ft) %>% 
  print(target = tmp)

# open word document
browseURL(tmp)

# G. Save as a Word table
word_document_name <-
  read_docx() %>%
  body_add_flextable(ft) %>%
  print(target = "document_with_summary.docx")

# H. Open the Word document 
browseURL(word_document_name)
