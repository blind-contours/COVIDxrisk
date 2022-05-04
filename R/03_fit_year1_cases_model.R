library(here)
source(here("R/utils_sl_varimp.R"))
source(here("R/util.R"))
cpus <- 20
plan(multisession, workers = cpus)

set.seed(5929922)
Cases1year <- fit_sl_varimp(outcome = "Casesat1year", label = "COVID-19 Cases at 1 Year", num_boot = 10)

