library(here)
source(here("R/utils_sl_varimp.R"))
source(here("R/util.R"))
cpus <- 25
plan(multisession, workers = cpus)

set.seed(5929922)

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
var_combn <- 2

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

saveRDS(sl, here(paste("Models/", outcome, ".RDS", sep = "")))

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

subcategories <- loaded_list$Subcategories
variable_list <- loaded_list$Variable_list
total_outcome <- loaded_list$total

risk_rescaled <- loaded_list$risk_rescaled
risk <- loaded_list$risk

load_model_time <- proc.time()

plan(multicore, workers = cpus)

var_imp_risk_results <- var_imp_risk(
  X = X, data = data, outcome = outcome,
  covars = covars, fit = sl,
  loss = loss_squared_error, Y = Y, num_boot = num_boot,
  Data_Dictionary = data_dictionary
)

variable_imp_risk_time <- proc.time()

var_imp_risk_results$Label <- data_dictionary$`Nice Label`[match(var_imp_risk_results$Variable, data_dictionary$`Variable Name`)]

saveRDS(var_imp_risk_results, here(paste("data/", outcome, "_ind_var_imp_risk.RDS", sep = "")))

print("Finished Risk Variable Importance")

subcat_imp_risk_results <- subcat_imp_risk(
  subcategories = subcategories,
  data = data, outcome = outcome,
  covars = covars,
  fit = sl,
  loss = loss_squared_error,
  Y = Y,
  num_boot = num_boot,
  variable_list = variable_list)

subcat_imp_risk_time <- proc.time()

saveRDS(subcat_imp_risk_results, here(paste("data/", outcome, "_subgroup_imp_risk.RDS", sep = "")))

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

saveRDS(var_imp_quantile_results, here(paste("data/", outcome, "_ind_var_imp_quantile.RDS", sep = "")))

print("Finished Quantile-Based Variable Importance")

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

print("Finished Sub-Category Quantile Importance")

saveRDS(subcat_imp_quantile_results, here(paste("data/", outcome, "_subgroup_imp_quant.RDS", sep = "")))


mips_results <- mips_imp_risk(risk_importance = var_imp_risk_results,
                          data = data,
                          outcome = outcome,
                          covars = covars,
                          fit = sl,
                          loss = loss_squared_error,
                          Y= Y,
                          num_boot = num_boot,
                          m = var_combn,
                          Data_Dictionary = data_dictionary,
                          p_val_fun = p_val_fun,
                          risk = risk)

print("Finished MIPS")

saveRDS(mips_results, here(paste("data/", outcome, "_intxn_imp_risk.RDS", sep = "")))

print(risk)
