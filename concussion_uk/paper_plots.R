## 01 Load packages ----------------------------------------
library(BGGM)
library(qgraph)

## 02 Load in concussion and control between estimates -----

data_dictionary <- readxl::read_excel("Data Dictionary.xlsx")
PCS_z <- readRDS("out/PCS_z.rds")
symptoms <- names(PCS_z)[19:40]
symptoms_names <- data_dictionary[data_dictionary$`Variable Name` %in% symptoms,]$`Variable Label`
symptoms_dict <- setNames(symptoms_names, symptoms)

## between concussion
fit_between_concussion <- readRDS("out/fit_between/fit_between_concussion.rds")
selected_between_concussion <- BGGM::select(fit_between_concussion)

## between control
fit_between_control <- readRDS("out/fit_between/fit_between_control.rds")
selected_between_control <- BGGM::select(fit_between_control)

## plots
# qgraph(cont, layout = L,
#        title="Contemporaneous", theme='colorblind', negDashed=FALSE,
#        groups=gr, legend.cex=0.4, details=TRUE, legend=TRUE, nodeNames = names)
# p1 <- qgraph(selected_between_concussion$pcor_adj, 
#   layout = 'circle', 
#   title="Contemporaneous",
#   labels = symptoms, 
#   theme = 'TeamFortress',
#   legend.cex=0.4,
#   label.font = 2,
#   fade = TRUE,
#   node.width = 0.85,
#   border.width = 1.5,
#   label.scale.equal = TRUE,
#   details=FALSE, 
#   legend=TRUE, 
#   nodeNames = symptoms_names
#   #,edge.width = 2.5
#       )
# 
# p0 <- qgraph(selected_between_control$pcor_adj, 
#              layout = 'circle', 
#              labels = symptoms, 
#              theme = 'gray',
#              label.font = 2,
#              fade = FALSE,
#              node.width = 1,
#              border.width = 1.9,
#              label.scale.equal = TRUE
#              #,edge.width = 2.5
# )
# 
# png(filename = "out/network_concusion.png",
#     width = 1200, height = 1200)
# plot(p1)
# dev.off()
# 
# png(filename = "out/network_control.png",
#     width = 1200, height = 1200)
# plot(p0)
# dev.off()

## Colored plots for publication
L <- averageLayout(selected_between_concussion$pcor_adj,
                   selected_between_control$pcor_adj)

png(filename = "out/figures/between_net.png",  width = 1400, height = 600,
    , res=100)

layout(matrix(c(1,1,2,2,2), nc=5, byrow = TRUE)) # 40% vs 60% widths
qgraph(selected_between_concussion$pcor_adj, 
       layout = L, 
       title="1A. Concussion",
       title.cex = 2,
       labels = symptoms, 
       theme = 'colorblind',
       # legend.cex=0.4,
       # label.font = 2,
       fade = TRUE,
       # node.width = 0.85,
       # border.width = 1.5,
       label.scale.equal = TRUE,
       vsize=12,
       details=FALSE, 
       legend=FALSE, 
       nodeNames = symptoms_names
       #,edge.width = 2.5
)
qgraph(selected_between_control$pcor_adj, 
       layout = L, 
       title="1B. Control",
       title.cex = 2,
       labels = symptoms, 
       theme = 'colorblind',
       label.scale.equal = TRUE,
       vsize=10,
       details=FALSE, 
       legend=TRUE, 
       legend.cex=0.6,
       asize=6,
       nodeNames = symptoms_names
)
dev.off()

## Plots for selected subjects ----

# First randomly select two participants from each group that completed
# at least 80% of the measurements

raw_data <- haven::read_sav("data/PCSNetworkModelAnalysis.sav") |> 
  janitor::clean_names()
converged <- as.numeric(gsub("[^0-9]", "", list.files("out/results_noise")))

set.seed(115)
ids_sample <- raw_data |> 
  dplyr::filter(subject_id %in% converged) |> 
  dplyr::with_groups(c(concussion,subject_id),
                     dplyr::count,head) |> 
  dplyr::filter(n >= max(n)*.8) |> 
  dplyr::slice_sample(n = 2, by = concussion)

# Full file paths
file_paths <- file.path("out/results_noise", paste0("id_", ids_sample$subject_id, ".rds"))
ex_files <- lapply(file_paths, readRDS)
names(ex_files) <- ids_sample$subject_id

selected <- lapply(ex_files, function(x) {
  BGGM::select(x[[1]])$pcor_weighted_adj
})

png(filename = "out/figures/ex_nets.png",  width = 6, height = 5, units = "in", res = 300)

layout(matrix(c(1,2,3,4), 2, 2))

# Set outer margins to allow room for the titles
par(oma = c(0, 0, 3, 0))  # bottom, left, top, right

for (i in 1:length(selected)) {
  qgraph(selected[[i]], 
         title=names(selected[i]),
         layout = L,
         title.cex = .95,
         labels = symptoms, 
         theme = 'colorblind',
         label.scale.equal = TRUE,
         vsize=10,
         details=FALSE, 
         legend=FALSE, 
         legend.cex=0.6,
         asize=6,
         nodeNames = symptoms_names
  )
}

# Add column titles (adjust line and adj as needed)
mtext("Control", side = 3, outer = TRUE, line = 1, at = 0.25, cex = 1)
mtext("Concussion", side = 3, outer = TRUE, line = 1, at = 0.75, cex = 1)

dev.off()
