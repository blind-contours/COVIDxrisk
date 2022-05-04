library(here)
source(here("R/utils_sl_varimp.R"))
source(here("R/util.R"))
cpus <- 20
plan(multisession, workers = cpus)

set.seed(5929922)
Deaths1year <- fit_sl_varimp(outcome = "Deathsat1year", label = "COVID-19 Deaths at 1 Year")

