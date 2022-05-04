library(here)
source(here("R/utils_sl_varimp.R"))
source(here("R/util.R"))
cpus <- 20
plan(multisession, workers = cpus)

set.seed(5929922)
Casesday100 <- fit_sl_varimp(outcome = "CountyRelativeDay100Cases", label = "COVID-19 Cases at Day 100")

