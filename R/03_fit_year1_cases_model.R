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
num_boot <- 3
var_combn <- 3

start_time <- proc.time()

# load in the data
load_data_results <- load_data(path_data = "cleaned_covid_data_final_Mar_31_22.csv", path_data_dict = "Data_Dictionary.xlsx")
data <- load_data_results$data
data_dictionary <- load_data_results$data_dictionary

# fit the SL model
sl_fit <- create_sl(data = data, outcome = outcome, all_outcomes = all_outcomes)
sl <- sl_fit$sl_fit
covars <- sl_fit$covars

fit_model_time <- proc.time()

print("Finished Model Fitting")

loaded_list <- load_model(
  fit = sl,
  loss = loss_squared_error,
  covars = covars,
  outcome = outcome,
  data = data,
  Data_Dictionary = data_dictionary
)

X <- loaded_list$X
Y <- loaded_list$Y
outcome <- loaded_list$outcome
subcategories <- loaded_list$Subcategories
variable_list <- loaded_list$Variable_list
total_outcome <- loaded_list$total
load_model_time <- proc.time()

X <- sample(X, 4)

plan(multicore, workers = cpus)

var_imp_risk_results <- var_imp_risk(
  X = X, data = data, outcome = outcome,
  covars = covars, fit = sl,
  loss = loss_squared_error, Y = Y, num_boot = num_boot,
  Data_Dictionary = data_dictionary
)

variable_imp_risk_time <- proc.time()

print("Finished Risk Variable Importance")

subcat_imp_risk_results <- subcat_imp_risk(
  subcategories = subcategories,
  data = data, outcome = outcome,
  covars = covars,
  fit = sl,
  loss = loss_squared_error,
  Y = Y,
  num_boot = num_boot,
  variable_list = variable_list
)

subcat_imp_risk_time <- proc.time()

print("Finished Risk Sub-Category Importance")

var_imp_quantile_results <- var_imp_quantile(
  X = X,
  data = data,
  outcome = outcome,
  covars = covars,
  fit = sl,
  loss = loss_squared_error,
  Y = Y,
  num_boot = num_boot,
  total,
  Data_Dictionary = data_dictionary,
  total = total_outcome,
  p_val_fun = p_val_fun
)

subcat_imp_quantile_results <- subcat_imp_quantile(subcategories,
  data = data,
  outcome = outcome,
  covars = covars,
  fit = sl,
  Y = Y,
  num_boot = num_boot,
  variable_list = variable_list,
  total = total_outcome,
  Data_Dictionary = data_dictionary,
  p_val_fun = p_val_fun
)

subcat_imp_quantile_results <- subcat_imp_quantile(subcategories,
                                                   data = data,
                                                   outcome = outcome,
                                                   covars = covars,
                                                   fit = sl,
                                                   Y = Y,
                                                   num_boot = num_boot,
                                                   variable_list = variable_list,
                                                   total = total_outcome,
                                                   Data_Dictionary = data_dictionary,
                                                   p_val_fun = p_val_fun)

mips_results <- mips_imp_risk(risk_importance,
                          data,
                          outcome,
                          covars,
                          fit,
                          loss,
                          Y,
                          num_boot,
                          m,
                          Data_Dictionary,
                          p_val_fun)

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
