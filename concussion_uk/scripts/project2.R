## Load packages --------------------------------------------------------------------

library(BGGM)
library(qgraph)
library(networktools)

## 00 Helper function
bridgestrenght <- function(x, ...){
  bridge(x, ...)$`Bridge Strength`
}

## 01 Load in data ------------------------------------------------------------------

d <- readRDS("data/BLSynth.rds")$syn

symptoms  <- names(d[,-1])

## 03 BGGM estimate ----------

net <- BGGM::estimate(d[,-1], type = "continuous", impute = TRUE
                      #, formula = ~ sex + covariate2
                      )

s_net <- summary(net)

net_sel <- BGGM::select(net)

## Edge weights --------------

s_net$dat_results

## Centrality

strength_centrality <- roll_your_own(object = net,
                                     FUN = bridgestrenght,
                                     communities = symptoms,
                                     select = TRUE,
                                     progress = TRUE)
strength_centrality

## Plot

plot(net_sel)
