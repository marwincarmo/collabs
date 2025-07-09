
# 0 Load required packages ------------------------------------------------

library(mlVAR)
library(imputeTS)
library(dplyr)
library(ggplot2)
library(stringr)
library(qgraph)

# 1 Load in data ----------------------------------------------------------

d <- readRDS("data/SynthData14.rds")
head(d$syn)

# 2 Select the most central symptoms --------------------------------------


sxd <- d$syn |> 
  dplyr::select(Main.SubjectIDNum, obs_num, Headache, PressureInHead, FeelingSlowedDown, FeelingLikeInAFog,
                DontFeelRight, DifficultyConcentrating, FatigueOrLowEnergy)


# 3 Impute missing data ---------------------------------------------------

## 3.1 Check missing pattern

sxd |> 
  group_by(Main.SubjectIDNum) |> 
  summarise(across(everything(),
                   ~sum(is.na(.x))/n())) |> 
  summarise(across(everything(),
                   ~mean(.x, na.rm = TRUE))) |> 
  # take average across all variables
  select(!Main.SubjectIDNum) |> 
  pull() |>
  mean()

## 3.2 Imput data with inputeTS()

sxd_imp <- sxd |> 
  group_by(Main.SubjectIDNum) |> 
  mutate(across(Headache:FatigueOrLowEnergy,
                ~imputeTS::na_kalman(.x, optim.control = 5000),
                .names = "{.col}_imp")) |> 
  ungroup() |> 
  select(Main.SubjectIDNum, obs_num, contains("_imp")) |>
  rename_with(~str_remove(.x, "_imp"), contains("_imp")) |>
  mutate(id = as.numeric(as.factor(Main.SubjectIDNum)))

## 3.2 Standardize the data

sxd_imp_std <- sxd_imp |> 
  dplyr::mutate(
    dplyr::across(Headache:FatigueOrLowEnergy, scale)
  )

# 4 mlVAR -----------------------------------------------------------------


vars <- colnames(sxd_imp[,3:9])

fit_mlvar <- mlVAR(data=sxd_imp, vars = vars, idvar = "Main.SubjectIDNum",
                   contemporaneous = "orthogonal", temporal = "orthogonal", estimator = "lmer", scale = TRUE)

summary(fit_mlvar)

# 5 Obtain centrality measures --------------------------------------------

# Obtain centrality of fitted mlVAR object
centrality_mlvar <- function(fit,
                             ignore_ar = TRUE){
  
  #--- Prepare
  n_var <- length(fit$fit$var)
  n_id <- length(fit$IDs)
  n_var_temp <- n_var
  
  #--- Obtain networks
  l_beta <- lapply(1:n_id, function(i){
    mlVAR::getNet(fit,
                  subject = i,
                  type = "temporal",
                  nonsig = "show")
  })
  l_pcor <- lapply(1:n_id, function(i){
    mlVAR::getNet(fit,
                  subject = i,
                  type = "contemporaneous",
                  nonsig = "show")
  })
  
  # Obtain overall adjacency matrix
  adj_beta <- ifelse(
    mlVAR::getNet(fit,
                  type = "temporal",
                  nonsig = "hide") != 0, 1, 0)
  adj_pcor <- ifelse(
    mlVAR::getNet(fit,
                  type = "contemporaneous",
                  nonsig = "hide") != 0, 1, 0) 
  
  # Set zero fixed effects to zero
  l_beta <- lapply(l_beta, function(x){
    x * adj_beta
  })
  l_pcor <- lapply(l_pcor, function(x){
    x * adj_pcor
  })
  
  
  #--- Density
  dens_temp <- lapply(l_beta, function(x){
    sum(abs(x))/(n_var^2)
  })
  dens_cont <- lapply(l_pcor, function(x){
    x <- x
    diag(x) <- 0
    x[lower.tri(x)] <- 0L
    sum(abs(x))/(n_var * (n_var-1)/2)
  })
  
  
  #--- Centrality
  # Important: in mlVAR, the lagged vars are rows, not columns
  # but only for the getNet function!
  outstrength <- lapply(l_beta, function(x){
    if(isTRUE(ignore_ar)){
      diag(x) <- 0
      n_var_temp <- n_var - 1
    }
    rowSums(abs(x))/n_var_temp
  })
  instrength <- lapply(l_beta, function(x){
    if(isTRUE(ignore_ar)){
      diag(x) <- 0
      n_var_temp <- n_var - 1
    }
    colSums(abs(x))/n_var_temp
  })
  strength <- lapply(l_pcor, function(x){
    # because diagonal is always zero, we subtract 1
    colSums(abs(x))/(n_var - 1)
  })
  
  # return list
  return(list(outstrength = outstrength,
              instrength = instrength,
              strength = strength,
              dens_temp = dens_temp,
              dens_cont = dens_cont))
  
}

centrality_res <- centrality_mlvar(fit_mlvar)


#
# Function to summarize the most frequent central symptom
freq_central <- function(centrality_output, var_names) {
  
  # --- Temporal Network (Outstrength and Instrength) ---
  
  # Most frequent central symptom based on outstrength
  outstrength_ranks <- lapply(centrality_output$outstrength, function(x) {
    names(x) <- var_names
    order(x, decreasing = TRUE)
  })
  outstrength_most_central <- unlist(lapply(outstrength_ranks, function(x) var_names[x[1]]))
  outstrength_summary <- as.data.frame(table(outstrength_most_central))
  names(outstrength_summary) <- c("Symptom", "Frequency")
  outstrength_summary$Proportion <- outstrength_summary$Frequency / sum(outstrength_summary$Frequency)
  outstrength_summary <- outstrength_summary[order(outstrength_summary$Frequency, decreasing = TRUE), ]
  
  # Most frequent central symptom based on instrength
  instrength_ranks <- lapply(centrality_output$instrength, function(x) {
    names(x) <- var_names
    order(x, decreasing = TRUE)
  })
  instrength_most_central <- unlist(lapply(instrength_ranks, function(x) var_names[x[1]]))
  instrength_summary <- as.data.frame(table(instrength_most_central))
  names(instrength_summary) <- c("Symptom", "Frequency")
  instrength_summary$Proportion <- instrength_summary$Frequency / sum(instrength_summary$Frequency)
  instrength_summary <- instrength_summary[order(instrength_summary$Frequency, decreasing = TRUE), ]
  
  
  # --- Contemporaneous Network (Strength) ---
  
  # Most frequent central symptom based on strength
  strength_ranks <- lapply(centrality_output$strength, function(x) {
    names(x) <- var_names
    order(x, decreasing = TRUE)
  })
  strength_most_central <- unlist(lapply(strength_ranks, function(x) var_names[x[1]]))
  strength_summary <- as.data.frame(table(strength_most_central))
  names(strength_summary) <- c("Symptom", "Frequency")
  strength_summary$Proportion <- strength_summary$Frequency / sum(strength_summary$Frequency)
  strength_summary <- strength_summary[order(strength_summary$Frequency, decreasing = TRUE), ]
  
  return(list(
    temporal_outstrength_summary = outstrength_summary,
    temporal_instrength_summary = instrength_summary,
    contemporaneous_strength_summary = strength_summary
  ))
}

## 5.1 Individual centrality indices ----

most_frequent_central <- freq_central(centrality_output = centrality_res, var_names = vars)
most_frequent_central

## 5.2 Centrality of the fixed effects -----

cont <- getNet(fit_mlvar, "contemporaneous", layout = "spring", nonsig = "hide", rule = "and")
temp <- getNet(fit_mlvar, "temporal", nonsig = "hide")

centrality_auto(cont)$node.centrality["Strength"]
centrality_auto(temp)$node.centrality[c("InStrength", "OutStrength")]


# 6 Plots -----------------------------------------------------------------

L <- averageLayout(cont, temp)

p_cont <- qgraph(cont, layout = L, negDashed=FALSE, legend.cex=0.4, details=FALSE,
       legend=FALSE, theme = "colorblind")
plot(p_cont)

p_temp <- qgraph(temp, layout = L, negDashed=FALSE, legend.cex=0.4, details=FALSE,
       legend=FALSE, theme = "colorblind")
plot(p_temp)

## 6.1 Saving the plots as tiff

tiff("plot_contemporaneous.tiff", res = 300, width = 1024, height = 768, compression = "lzw")
plot(p_cont)
dev.off()

tiff("plot_temporal.tiff", res = 300, width = 1024, height = 768, compression = "lzw")
plot(p_temp)
dev.off()


# 7 dcnet ------------------------------------------------------------------
## library(dcnet)
## install.packages("data.table")
## library(data.table)

## head(sxd)

## N <- length(unique(sxd_imp_std$id))
## N

## ids <- unique(sxd_imp_std$id)
## tl <- 14
## groupvec <- rep(seq_len(N),  each = tl)
## dt <- data.table::data.table(sxd_imp_std )
## variables <- names(sxd_imp_std)[3:8]

## tsdat <- lapply(ids, function(x) dt[dt$id == x, .SD, .SDcols = variables] )
## tsdat

## #devtools::load_all(path = "./")

## fit_init <- dcnet(data = tsdat, J = N, group = groupvec, S_pred = NULL,
##                   standardize_data = FALSE,
##                   parameterization = "CCC",
##                   iterations = 30000,
##                   sampling_algorithm = "variational",
##                   meanstructure = "VAR",
##                   chains = 4,
##                   eta=0.1,
##                   threads = 8,
##                   grainsize = 5,
##                   init = 0)

## fit_init$model_fit$output()
## grep("vec_phi_fixed", colnames(fit_init$model_fit$draws()))
## fit_init$model_fit$draws()[, 3062:3070]
## variables
## phi <- matrix(colMeans(fit_init$model_fit$draws()[, grep("vec_phi_fixed", colnames(fit_init$model_fit$draws()))]), ncol = length(variables))
## round(phi, 4)
## #pdf(file = "~/Nextcloud/Documents/eric.pdf")
## qgraph::qgraph(phi, labels = variables)
## #dev.off()

## fit <- dcnet(
##     data = tsdat, parameterization = "CCC", J = N,
##     group = groupvec, standardize_data = FALSE,
##     init = fit_init$model_fit,
##     meanstructure = "VAR",
##     iterations = 50000,
##     sampling_algorithm = "variational",
##     algorithm = "fullrank", ## fullrank should be less biased
##     grad_samples = 1,
##     elbo_samples = 150,
##     eta = 0.1,
##     adapt_iter = 100,
##     grainsize = 40,
##     chains = 4,
##     threads_per_chain = 4)



##     # Simulate data:
##      Model <- mlVARsim(nPerson = 50, nNode = 3, nTime = 50, lag=1)
##      Model$Data = Model$Data*10
##      # Estimate using correlated random effects:
##      fit1 <- mlVAR(Model$Data, vars = Model$vars, idvar = Model$idvar, lags = 1, temporal = "correlated")
## summary(fit1)

## N <- 50
## N
## ids <- unique(Model$Data$ID)/10
## ids
## tl <- 50
## groupvec <- rep(seq_len(N),  each = tl)
## groupvec
## dt <- data.table::data.table(Model$Data )
## dt$ID=dt$ID/10
## dt
## variables <- names(dt)[1:3]
## variables

## tsdat <- lapply(ids, function(x) dt[dt$ID == x, .SD, .SDcols = variables] )
## tsdat

## #devtools::load_all(path = "./")

## fit_init <- dcnet(data = tsdat, J = N, group = groupvec, S_pred = NULL,
##                   standardize_data = TRUE,
##                   parameterization = "CCC",
##                   iterations = 30000,
##                   sampling_algorithm = "variational",
##                   meanstructure = "VAR",
##                   chains = 4,
##                   eta=0.1,
##                   threads = 8,
##                   grainsize = 5,
##                   init = 0)

## phi <- matrix(colMeans(fit_init$model_fit$draws()[, grep("vec_phi_fixed", colnames(fit_init$model_fit$draws()))]), ncol = length(variables))
## round(phi, 4)
