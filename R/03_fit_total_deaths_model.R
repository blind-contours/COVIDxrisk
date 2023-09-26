library(here)
source(here("R/utils_sl_varimp.R"))
source(here("R/util.R"))

# Configuration
cpus <- 20
set.seed(59293322)
outcome <- "TotalDeathsUpToDate"
num_boot <- 10
var_combn <- 2
n_folds <- 10
quantile_threshold <- 0.75


# Initialize parallel processing
plan(multisession, workers = cpus, gc = TRUE)

# Load data and prepare for CV
load_data_results <- load_data(path_data = "cleaned_covid_data_final.csv", path_data_dict = "Data_Dictionary.xlsx")
data <- load_data_results$data
data$fold <- sample(rep(1:n_folds, length.out = nrow(data)))  # Assign fold
data_dictionary <- load_data_results$data_dictionary

# Lists to save results
var_imp_results_list <- vector("list", n_folds)
subcat_imp_results_list <- vector("list", n_folds)
mips_imp_results_list <- vector("list", n_folds)

all_outcomes <- c(
  "CountyRelativeDay100Cases",
  "TotalCasesUpToDate",
  "CountyRelativeDay100Deaths",
  "TotalDeathsUpToDate",
  "Deathsat1year",
  "Casesat1year"
)

# Cross-Validation Loop
for (fold in 1:n_folds) {
  cat("Starting fold", fold, "\n")

  # Data Splitting
  train_data <- data[data$fold != fold,]
  validation_data <- data[data$fold == fold,]

  # Model Fitting
  sl_fit <- create_sl(data = train_data, outcome = outcome, all_outcomes = all_outcomes,
                      quantile_threshold = quantile_threshold)

  top_predictors <- sl_fit$top_vars
  covars <- sl_fit$covars
  best_learner <- sl_fit$best_learner

  # Load Model
  loaded_list <- load_model(
    fit = sl_fit$sl_fit,
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

  var_imp_fold_results <- var_imp_quantile(X = top_predictors,
                                           data = validation_data,
                                           outcome = outcome,
                                           covars = covars,
                                           fit = best_learner,
                                           loss = loss_squared_error,
                                           Y = Y,
                                           num_boot = num_boot,
                                           Data_Dictionary = data_dictionary,
                                           p_val_fun = p_val_fun)

  var_imp_fold_results$fold <- fold


  # Variable Importance
  var_imp_results_list[[fold]] <- var_imp_fold_results

  cat("Finished Quantile-Based Variable Importance for fold", fold, "\n")

  # Sub-Category Importance
  subcat_fold_results <- subcat_imp_quantile(subcategories,
                                             data = validation_data,
                                             outcome = outcome,
                                             covars = covars,
                                             fit = best_learner,
                                             Y = Y,
                                             num_boot = num_boot,
                                             variable_list = variable_list,
                                             Data_Dictionary = data_dictionary,
                                             p_val_fun = p_val_fun)

  subcat_fold_results$fold <- fold
  subcat_imp_results_list[[fold]] <- subcat_fold_results

  cat("Finished Sub-Category Quantile Importance for fold", fold, "\n")

  # Interaction Quantile Importance
  mips_fold_results <- mips_imp_quantile(quantile_importance = var_imp_fold_results,
                                         data = validation_data,
                                         outcome = outcome,
                                         covars = covars,
                                         fit = best_learner,
                                         loss = loss_squared_error,
                                         Y = Y,
                                         num_boot = num_boot,
                                         m = var_combn,
                                         Data_Dictionary = data_dictionary,
                                         p_val_fun = p_val_fun,
                                         total = total_outcome)

  mips_fold_results$fold <- fold

  mips_imp_results_list[[fold]] <- mips_fold_results

  cat("Finished Interaction Quantile Importance for fold", fold, "\n")

  gc()  # Memory cleanup
}

var_imp_fold_df <- do.call(rbind, var_imp_results_list)
subcat_fold_df <- do.call(rbind, subcat_imp_results_list)
mips_imp_results_df <- do.call(rbind, mips_imp_results_list)

# Save results
saveRDS(var_imp_fold_df, here("data/all_var_imp_results.RDS"))
saveRDS(subcat_fold_df, here("data/all_subcat_imp_results.RDS"))
saveRDS(mips_imp_results_df, here("data/all_mips_imp_results.RDS"))



