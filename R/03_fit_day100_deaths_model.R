library(here)
source(here("R/utils_sl_varimp.R"))
source(here("R/util.R"))
cpus <- 20
plan(multisession, workers = cpus)

set.seed(5929922)
Deathsday100 <- fit_sl_varimp(outcome = "CountyRelativeDay100Deaths", label = "COVID-19 Deaths at Day 100", num_boot = 10)

