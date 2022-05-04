library(here)
source(here("R/utils_sl_varimp.R"))
source(here("R/util.R"))
cpus <- 20
plan(multisession, workers = cpus)

set.seed(5929922)
Casestodate <- fit_sl_varimp(outcome = "TotalCasesUpToDate", label = "Total COVID-19 Cases To-Date")

