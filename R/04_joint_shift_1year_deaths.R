library(here)
source(here("R/util.R"))
cpus <- 5
mcoptions <- list(preschedule=FALSE, set.seed=FALSE)

doParallel::registerDoParallel(cpus)
future::plan(future::multisession, gc = TRUE)

## SL result

ML_pipeline_results <- readRDS(here("Models/Deathsat1year.RDS"))
outcome <- "Deathsat1year_PopScale"
################ Define Global Variables ##################

SCALE <- FALSE

data_original <- read.csv(PROCESSED_DATA_PATH("cleaned_covid_data_final_Mar_31_22.csv"), check.names = FALSE)
data_original <- data_original[, -1] # remove the column with empty name

## read in data dictionary for identifying subgroups of top variables to isolate the different control conditions
Data_Dictionary <- read_excel(PROCESSED_DATA_PATH("Data_Dictionary.xlsx"))
Data_Dictionary_Used <- Data_Dictionary %>%
  filter(Keep == "Yes") %>%
  select(`Variable Name`, `Label`, `Nice Label`)

## remove from the list covariates that had too many NAs and were then dropped before analysis, FIPS, and the outcome data:
vars_rmv_na <- read.csv(PROCESSED_DATA_PATH("vars_removed_na_thresh_combined.csv"))
vars_rmv_na <- as.vector(vars_rmv_na$x)

# 6 types of outcomes
initial_outcomes <- c(
  "CountyRelativeDay100Cases",
  "TotalCasesUpToDate",
  "CountyRelativeDay100Deaths",
  "TotalDeathsUpToDate",
  "Deathsat1year",
  "Casesat1year"
)

removing <- c(
  vars_rmv_na,
  initial_outcomes
)

Data_Dictionary_Used <- Data_Dictionary_Used[-match(removing, Data_Dictionary_Used$`Variable Name`), ]
variable_list <- Data_Dictionary_Used$`Variable Name`
subcategory_list <- Data_Dictionary_Used$Label


## create the per-capita outcomes
data_original <- data_original %>% mutate(across(contains("RelativeDay") | contains("UpToDate") | contains("at1year"),
  .fns = list(PopScale = function(x) x / data_original$Population),
  .names = "{col}_{fn}"
))

# test perc reduced
percents <- seq(0, 0.9, by = 0.1)

pop_scales <- names(data_original)[grepl("PopScale", names(data_original))]

all_outcomes <- c(
  initial_outcomes,
  pop_scales
)


target_outcomes <- c(pop_scales)

## set up covar features
covars <- colnames(data_original)[-which(names(data_original) %in% c(
  all_outcomes,
  "X1",
  "fips",
  "county_names",
  "countyFIPS"
))]


get_top_variables <- function(ML_results) {
  varimp_data <- ML_results$var_imp
  highest_rr_idx <- which.max(varimp_data$joint_results$diff) # TODO: use diff from joint_results
  top_var <- varimp_data$joint_results$`Variable Combo`[highest_rr_idx]
  return(top_var)
}

top_vars <- get_top_variables(ML_pipeline_results)

###############################################
##### TOP VARS & CAT FOR EACH OUTCOME #########
###############################################

## get the top variables, their subcategories, and accompanying variables in same category for marginal predictions
top_vars <- strsplit(top_vars, " & ")

# convert them from nice label to variable label
top_vars <- variable_list[match(top_vars[[1]], Data_Dictionary_Used$`Nice Label`)]
top_var_subgroups <- subcategory_list[match(top_vars, variable_list)]
top_var_subcat_vars <- purrr::map(.x = top_var_subgroups, ~ variable_list[subcategory_list %in% .x])


## set up bootstrap CI function
bootstrapCI <- function(target_variable,
                        data_original,
                        ML_pipeline_results,
                        covars,
                        outcome,
                        perc,
                        marginal_directions) {

  future::plan(future::sequential, gc = TRUE)

  sl <- ML_pipeline_results$fit

  nr <- nrow(data_original)
  data_tmp <- data_original
  resampled_data <- data_tmp[sample(1:nr, size = nr, replace = TRUE), ]

  task <- make_sl3_Task(
    data = resampled_data,
    outcome = outcome,
    covariates = covars,
    folds = origami::make_folds(resampled_data, fold_fun = folds_vfold, V = 2)
  )

  sl_fit_full_resampled <- suppressMessages(sl$train(task))

  ## get the original data and reduce the target variable by perc
  data_resampled_reduced <- resampled_data

  if (length(target_variable) > 1) {
    for (t_var_index in 1:length(target_variable)) {
      t_var <- target_variable[[t_var_index]]
      t_var_dir <- marginal_directions[[t_var_index]]

      if (t_var_dir == 1) {
        t_thresh <- min(data_resampled_reduced[, t_var])
        data_resampled_reduced[[t_var]] <- ifelse(t_thresh < data_resampled_reduced[[t_var]], data_resampled_reduced[[t_var]] - (data_resampled_reduced[[t_var]] * perc), t_thresh)
      } else {
        t_thresh <- max(data_resampled_reduced[, t_var])
        perc <- -perc
        data_resampled_reduced[[t_var]] <- ifelse(t_thresh > data_resampled_reduced[[t_var]], data_resampled_reduced[[t_var]] - (data_resampled_reduced[[t_var]] * perc), t_thresh)
      }
    }
  } else {
    t_var <- target_variable
    t_thresh <- min(data_resampled_reduced[, t_var])
    data_resampled_reduced[[t_var]] <- ifelse(t_thresh < data_resampled_reduced[[t_var]], data_resampled_reduced[[t_var]] - (data_resampled_reduced[[t_var]] * perc), t_thresh)
  }

  reduced_tasks <- make_sl3_Task(
    data = data_resampled_reduced,
    outcome = outcome,
    covariates = covars,
    folds = origami::make_folds(data_resampled_reduced, fold_fun = folds_vfold, V = 2)
  )

  ## predict through superlearner for reduced data on resampled models
  sl_preds_reduced_full <- sl_fit_full_resampled$predict(reduced_tasks)

  results <- data.frame(sl_preds_reduced_full)
  colnames(results) <- c("SL_full_model")

  return(results)
}

bootstrap_marginal_predictions <- function(target_variable,
                                           ML_pipeline_results,
                                           outcome,
                                           data_original = data_original,
                                           covars = covars,
                                           percents = percents,
                                           boot_num) {
  pop <- data_original$Population

  target_variable <- as.list(target_variable)
  target_variable[[length(target_variable) + 1]] <- target_variable

  boot_df_SL_full <- replicate(
    n = length(target_variable),
    expr = {
      as.data.frame(matrix(nrow = length(percents), ncol = 5))
    },
    simplify = F
  )

  boot_array_list <- list()
  marginal_direction_list <- list()
  marginal_directions <- list()

  for (var_index in 1:length(target_variable)) {
    var <- unlist(target_variable[[var_index]])

    if (length(var) > 1) {
      var_label <- paste(var, collapse = " & ")
    } else {
      var_label <- var
    }

    for (i in 1:length(percents)) {
      perc <- percents[i]

      remaining <- boot_num

      # while (remaining > 0) {

        # bootstrap for boot_num number of times
        boot_updates <- foreach(
          this_iter = seq(remaining),
          .errorhandling = "pass",
          .options.multicore=mcoptions
        ) %dopar% {
          bootstrapCI(
            target_variable = var,
            data_original = data_original,
            ML_pipeline_results = ML_pipeline_results,
            covars = covars,
            outcome = outcome,
            perc = perc,
            marginal_directions
          )
        }

        boot_updates <- boot_updates[lapply(boot_updates, is.null) == FALSE]

      #   remaining <- boot_num - length(boot_updates)
      #
      #   print(remaining)
      # }

      boot_data_array_SL_full <- colSums(do.call(cbind, boot_updates) * pop)
      boot_data_vector_SL_full <- do.call(cbind, boot_updates) * pop

      boot_array_list[[var_index]] <- boot_data_vector_SL_full

      probs <- c(0.025, 0.50, 0.975)
      SL_quantiles <- quantile(boot_data_array_SL_full, probs = probs, na.rm = TRUE)

      colnames(boot_df_SL_full[[var_index]]) <- c("Perc", "Boot Low", "Boot Pred", "Boot High", "Variable")

      boot_df_SL_full[[var_index]][i, 1] <- perc
      boot_df_SL_full[[var_index]][i, 2] <- SL_quantiles[1]
      boot_df_SL_full[[var_index]][i, 3] <- SL_quantiles[2]
      boot_df_SL_full[[var_index]][i, 4] <- SL_quantiles[3]
      boot_df_SL_full[[var_index]][i, 5] <- var_label

      max_diff <- boot_df_SL_full[[var_index]]$`Boot Pred`[1] - boot_df_SL_full[[var_index]]$`Boot Pred`[10]
      direction <- ifelse(max_diff > 0, 1, -1)

      marginal_directions[[var_index]] <- direction
    }
  }

  boot_results <- do.call(rbind, boot_df_SL_full)

  actual_outcome <- sum(data_original[[strsplit(outcome, "_")[[1]][1]]])

  plot <- ggplot(boot_results, aes(x = `Perc`, y = `Boot Pred`, colour = Variable)) +
    geom_line(aes(linetype = `Variable`, colour = `Variable`)) +
    geom_point(aes(y = `Boot Pred`, colour = `Variable`)) +
    geom_errorbar(
      aes(ymin = `Boot Low`, ymax = `Boot High`, group = `Variable`),
      width = 0.05
    ) +
    xlab("Percent Reduced") +
    ylab("Model Predictions") +
    geom_hline(yintercept = actual_outcome)

  figure_file_name <- paste(outcome, "_marginal_predictions", ".png", sep = "")
  CI_data_file_name <- here("data/processed", paste(paste(outcome, "CI_data", sep = "_"), ".csv", sep = ""))
  array_data_file_name <- here("data/processed", paste(paste(outcome, "full_data", sep = "_"), ".RDS", sep = ""))

  ggsave(
    filename = figure_file_name,
    plot = plot,
    device = NULL,
    path = here("Figures"),
    scale = 1,
    width = 8,
    height = 5,
    units = c("in"),
    dpi = 300,
    limitsize = TRUE
  )

  write.csv(boot_results, file = CI_data_file_name)
  saveRDS(boot_array_list, file = array_data_file_name)

  return(NULL)
}

joint_impact_day100_cases <- bootstrap_marginal_predictions(
  target_variable = top_vars,
  ML_pipeline_results = ML_pipeline_results,
  outcome = outcome,
  data_original = data_original,
  covars = covars,
  percents = percents,
  boot_num = 5
)
