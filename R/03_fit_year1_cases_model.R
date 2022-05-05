library(here)
source(here("R/utils_sl_varimp.R"))
source(here("R/util.R"))
cpus <- 20
plan(multisession, workers = cpus)

set.seed(5929922)
# Cases1year <- fit_sl_varimp(outcome = "Casesat1year", label = "COVID-19 Cases at 1 Year", num_boot = 10)

# set the fit_sl_varimp args
outcome <- "Casesat1year"
label <- "COVID-19 Cases at 1 Year"
num_boot <- 10

# extract super learner fitting procedure
## load data
covid_data_processed <- read.csv(PROCESSED_DATA_PATH("cleaned_covid_data_final_Mar_31_22.csv"), check.names = FALSE)
covid_data_processed <- covid_data_processed[, -1]

Data_Dictionary <- read_excel(PROCESSED_DATA_PATH("Data_Dictionary.xlsx"))
Data_Dictionary <- Data_Dictionary[Data_Dictionary$`Variable Name` %in% colnames(covid_data_processed), ]

covid_data_processed$CountyRelativeDay100Cases <- covid_data_processed$CountyRelativeDay100Cases / covid_data_processed$Population
covid_data_processed$TotalCasesUpToDate <- covid_data_processed$TotalCasesUpToDate / covid_data_processed$Population
covid_data_processed$CountyRelativeDay100Deaths <- covid_data_processed$CountyRelativeDay100Deaths / covid_data_processed$Population
covid_data_processed$TotalDeathsUpToDate <- covid_data_processed$TotalDeathsUpToDate / covid_data_processed$Population
covid_data_processed$Deathsat1year <- covid_data_processed$Deathsat1year / covid_data_processed$Population
covid_data_processed$Casesat1year <- covid_data_processed$Casesat1year / covid_data_processed$Population

outcomes <- c(
  "CountyRelativeDay100Cases",
  "TotalCasesUpToDate",
  "CountyRelativeDay100Deaths",
  "TotalDeathsUpToDate",
  "Deathsat1year",
  "Casesat1year"
)


covars <- colnames(covid_data_processed)[-which(names(covid_data_processed) %in% c(
  outcomes,
  "fips",
  "county_names"
))]

task <- make_sl3_Task(
  data = covid_data_processed,
  covariates = covars,
  outcome = outcome,
  folds = origami::make_folds(covid_data_processed, fold_fun = folds_vfold, V = 10)
)

discrete_sl <- source(here("R/utils_create_sl.R"))
discrete_sl <- discrete_sl$value()
## fit the sl3 object
sl_fit <- discrete_sl$train(task)

# run varimp function
var_importance <- run_varimp(
  fit = sl_fit,
  loss = loss_squared_error,
  covars = covars,
  outcome = outcome,
  data = covid_data_processed,
  Data_Dictionary = Data_Dictionary,
  label = label,
  num_boot = num_boot
)

# save the results
SL_results <- list("fit" = sl_fit, "var_imp" = var_importance)
saveRDS(SL_results, here(paste("Models/", outcome, ".RDS", sep = "")))
