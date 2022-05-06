library(here)
source(here("R/utils_sl_varimp.R"))
source(here("R/util.R"))
cpus <- 5
plan(multisession, workers = cpus)

set.seed(5929922)
# Cases1year <- fit_sl_varimp(outcome = "Casesat1year", label = "COVID-19 Cases at 1 Year", num_boot = 10)

# set the fit_sl_varimp args
outcome <- "Casesat1year"

all_outcomes <- c(
  "CountyRelativeDay100Cases",
  "TotalCasesUpToDate",
  "CountyRelativeDay100Deaths",
  "TotalDeathsUpToDate",
  "Deathsat1year",
  "Casesat1year"
)
label <- "COVID-19 Cases at 1 Year"
num_boot <- 10
var_combn <- 3

start <- proc.time()

# load in the data
load_data_results <-  load_data(path_data = "cleaned_covid_data_final_Mar_31_22.csv", path_data_dict = "Data_Dictionary.xlsx")
data <- load_data_results$data
data_dictionary <- load_data_results$data_dictionary

# fit the SL model
sl_fit <- create_sl(data = data, outcome = outcome, all_outcomes = all_outcomes)
sl <- sl_fit$sl_fit
covars <- sl_fit$covars

plan(multicore, workers = cpus)

var_imp_results <- run_varimp(fit = sl,
           loss = loss_squared_error,
           covars= covars,
           outcome = outcome,
           data = data,
           Data_Dictionary = data_dictionary,
           label = label,
           num_boot = num_boot,
           m = var_combn)

proc.time() - start

# save model
saveRDS(var_imp_results$fit, here(paste("Models/", outcome, ".RDS", sep = "")))

# save ind var risk data
saveRDS(var_imp_results$var_imp$Var_Risk_Results, here(paste("data/", outcome, "_ind_var_imp_risk.RDS", sep = "")))

# save subgroup risk data
saveRDS(var_imp_results$var_imp$Subgroup_Risk_Results, here(paste("data/", outcome, "_subgroup_imp_risk.RDS", sep = "")))

# save ind quantile pred data
saveRDS(var_imp_results$var_imp$Var_Quantile_Results, here(paste("data/", outcome, "_ind_var_imp_quantile.RDS", sep = "")))

# save subgroup quantile pred data
saveRDS(var_imp_results$var_imp$Subgroup_Quantile_Results, here(paste("data/", outcome, "_subgroup_imp_risk.RDS", sep = "")))

# save subgroup quantile pred data
saveRDS(var_imp_results$var_imp$Intxn_Risk_Results, here(paste("data/", outcome, "_intxn_imp_risk.RDS", sep = "")))

# print model risk scaled back to units
print(var_imp_results$var_imp$model_risk)



