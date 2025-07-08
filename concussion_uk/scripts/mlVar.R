
# -------------------------------------------------------------------------
# --------------- 1. Loading packages & Data ------------------------------
# -------------------------------------------------------------------------

library("lme4")
library("ggplot2")
library("dplyr")
library("tidyr")
library("qgraph")
library("tidyselect")
library("mlVAR")
library("viridis")
#library("summarytools")
library("lm.beta")

## Load in data
PCS_z <- readRDS("out/PCS_z.rds")
data_dictionary <- readxl::read_excel("Data Dictionary.xlsx")

symptoms <- names(PCS_z)[19:40]
symptoms_names <- data_dictionary[data_dictionary$`Variable Name` %in% symptoms,]$`Variable Label`
symptoms_dict <- setNames(symptoms_names, symptoms)
head(PCS_z)

pcs <- PCS_z |> 
  dplyr::select(subject_id, survey_day, survey_time,
                head:vis) |> 
  dplyr::mutate(survey_time = dplyr::case_when(
    survey_time == 10 ~ 0,
    survey_time == 3 ~ 1,
    survey_time == 8 ~ 2
  ))

## Run mlVAR (orthogonal):
# res <- mlVAR(pcs,
#              vars=symptoms,
#              idvar="subject_id",
#              dayvar="survey_day",
#              beepvar="survey_time",
#              lags = 1,
#              temporal = "orthogonal",
#              contemporaneous = "orthogonal",
#              nCores = 16)

#save(res, file="out/network_orthogonal.RData")

load("out/network_orthogonal.RData")

# Get networks:
cont <- getNet(res, "contemporaneous", layout = "circle", nonsig = "hide", rule = "and")
bet  <- getNet(res, "between", nonsig = "hide", rule = "and")
temp <- getNet(res, "temporal", nonsig = "hide")

# layout:
L <- averageLayout(cont,bet,temp)

## Contemporaneous
qgraph(cont, layout = L,
              title="Contemporaneous", theme='gray', #negDashed=FALSE,
              legend.cex=0.4, details=TRUE, legend=TRUE, nodeNames = symptoms_names)

## Temporal
qgraph(temp, layout = L,
       title="Temporal, Lag-1", theme='colorblind', negDashed=FALSE,
       legend.cex=0.4, details=TRUE, legend=TRUE, nodeNames = symptoms_names)
# temporal no diag
qgraph(temp, layout = L,
       title="Temporal, Lag-1",  theme='colorblind', negDashed=FALSE,
       legend.cex=0.4, details=TRUE, diag = FALSE, legend=TRUE, nodeNames = symptoms_names)

## Between
qgraph(bet, layout = L,
       title="Between-Subjects",theme='colorblind', negDashed=FALSE,
       legend.cex=0.4, details=TRUE, legend=TRUE, nodeNames = symptoms_names)
