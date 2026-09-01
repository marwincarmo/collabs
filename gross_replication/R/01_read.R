## -----------------------------------------------------------------------------
## 01  Read the raw exports and look at what is in them.
##
## Run this first, on its own, whenever a new export arrives. The point is to
## SEE the questions and their answer options before recoding anything: if the
## survey changed, it shows up here.
## -----------------------------------------------------------------------------

source("R/00_setup.R")

person_raw <- read_expiwell(raw_person_file)
esm_raw    <- read_expiwell(raw_esm_file)

cat("\n---- Day-1 survey ----------------------------------------------------\n")
cat(nrow(person_raw$data), "rows,",
    n_distinct(person_raw$data$`Participant ID`), "participants\n\n")
print(person_raw$codebook, n = Inf)

cat("\n---- ESM prompts -----------------------------------------------------\n")
cat(nrow(esm_raw$data), "rows,",
    n_distinct(esm_raw$data$`Participant ID`), "participants\n\n")
print(esm_raw$codebook, n = Inf)

## The answer options, in the order the survey defines them. Recoding in 02
## relies on this order, so it is worth reading once.
cat("\n---- ESM answer options ----------------------------------------------\n")
for (q in names(esm_raw$choices)) {
  if (is.null(esm_raw$choices[[q]])) next
  cat("\n", q, ": ", str_trunc(esm_raw$codebook$question[esm_raw$codebook$column == q], 60), "\n", sep = "")
  cat(paste0("   ", seq_along(esm_raw$choices[[q]]), " = ",
             str_trunc(esm_raw$choices[[q]], 65)), sep = "\n")
}

cat("\n---- Metadata columns ------------------------------------------------\n")
print(names(esm_raw$data)[!grepl("^Q\\d+$", names(esm_raw$data))])

## Day of Survey: 1 is the day-1 survey, 2-8 are the seven ESM days.
cat("\nESM prompts per Day of Survey:\n")
print(table(esm_raw$data$`Day of Survey`, useNA = "ifany"))

## The platform is supposed to number the prompts within each day. In the pilot
## it does not -- every row says 1 -- so 02 derives the order from Start Date.
cat("\nOccasion within Day:\n")
print(table(esm_raw$data$`Occasion within Day`, useNA = "ifany"))

## Personal data. These columns are dropped in 02 and never reach any output,
## but they are in the raw files, so those must not be shared or committed.
cat("\nColumns holding personal data (dropped in 02):\n")
print(intersect(c("ip address", "Location - Lat", "Location - Long"),
                names(esm_raw$data)))

saveRDS(list(person = person_raw, esm = esm_raw),
        file.path(dir_derived, "01_raw.rds"))
cat("\nSaved:", file.path(dir_derived, "01_raw.rds"), "\n")
