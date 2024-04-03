  # packages
  library(tidyverse)
  library(haven)
  library(psych)
  library(survey)
  library(readxl)
  library(GPArotation)
  library(factoextra)
  
  # load household data 
  #setwd("Users/bapu/Ω-personal/friends/wfp/afghanistan")
  afsms <- read_sav("afghanistan/AFSMS HH Dataset - Anonymized.sav")

  # load weights and ipc file
  weights <- read_xlsx("afghanistan/Weight file_hhs.xlsx", sheet = 3)
  ipc <- read.csv("afghanistan/ipc_2022.csv")
  
  # cleaning
  afsms <- afsms %>% 
    mutate(rCSI = as.double(rCSI),
           strat = paste0(Prov, "_", CommType)) %>% 
    left_join(weights, by = "strat") %>% 
    mutate(rep_pop = Weight_hh * avg_hh_size) %>%  # representative population
    left_join(ipc, by = "strata")
  
  # survey design
  afsms_svy <- svydesign(ids=~ClustCode, 
                         strata = ~strata,
                         weights = ~HHWeight_R, 
                         nest = T, 
                         data = afsms)
  
  pairs.panels(afsms[,c("FCS", "HDDS", "rCSI", "HHS")],
               smooth = TRUE,      # If TRUE, draws loess smooths
               scale = FALSE,      # If TRUE, scales the correlation text font
               density = TRUE,     # If TRUE, adds density plots and histograms
               ellipses = TRUE,    # If TRUE, draws ellipses
               method = "spearman", # Correlation method (also "spearman" or "kendall")
               pch = 20,           # pch symbol
               lm = FALSE,         # If TRUE, plots linear fit rather than the LOESS (smoothed) fit
               cor = TRUE,         # If TRUE, reports correlations
               hist.col = 8)       # Histograms color
  
  
  # scale data
  afsms_scale <- scale(afsms[,c("FCS", "HDDS", "rCSI", "HHS")])
  
  # factor adequacy
  KMO(afsms_scale)
  
  # Parellel analysis
  parallel_afsms <- fa.parallel(afsms_scale,
                                fm = 'ml',
                                fa = 'fa',
                                n.iter = 50,
                                SMC = TRUE,
                                quant = .95)
  parallel_afsms
  
  vss_afsms <- vss(afsms_scale,
                   rotate = 'oblimin',
                   fm = 'gls')
  
  vss_afsms
  
# constrain the residual variances
# 
  
  
  