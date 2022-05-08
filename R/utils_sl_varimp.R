create_sl <- function(data = covid_data_processed, outcome = outcome, all_outcomes = all_outcomes) {
  covars <- colnames(data)[-which(names(data) %in% c(
    all_outcomes,
    "fips",
    "county_names"
  ))]

  task <- make_sl3_Task(
    data = data,
    covariates = covars,
    outcome = outcome,
    folds = origami::make_folds(data, fold_fun = folds_vfold, V = 10)
  )

  discrete_sl <- source(here("R/utils_create_sl.R"))
  discrete_sl <- discrete_sl$value()
  ## fit the sl3 object
  sl_fit <- discrete_sl$train(task)

  return(list("sl_fit" = sl_fit, "covars" = covars))
}

load_data <- function(path_data = "cleaned_covid_data_final_Mar_31_22.csv", path_data_dict = "Data_Dictionary.xlsx") {
  ## load data
  covid_data_processed <- read.csv(PROCESSED_DATA_PATH(path_data), check.names = FALSE)
  covid_data_processed <- covid_data_processed[, -1]

  Data_Dictionary <- read_excel(PROCESSED_DATA_PATH(path_data_dict))
  Data_Dictionary <- Data_Dictionary[Data_Dictionary$`Variable Name` %in% colnames(covid_data_processed), ]

  covid_data_processed$CountyRelativeDay100Cases <- covid_data_processed$CountyRelativeDay100Cases / covid_data_processed$Population
  covid_data_processed$TotalCasesUpToDate <- covid_data_processed$TotalCasesUpToDate / covid_data_processed$Population
  covid_data_processed$CountyRelativeDay100Deaths <- covid_data_processed$CountyRelativeDay100Deaths / covid_data_processed$Population
  covid_data_processed$TotalDeathsUpToDate <- covid_data_processed$TotalDeathsUpToDate / covid_data_processed$Population
  covid_data_processed$Deathsat1year <- covid_data_processed$Deathsat1year / covid_data_processed$Population
  covid_data_processed$Casesat1year <- covid_data_processed$Casesat1year / covid_data_processed$Population

  return(list("data" = covid_data_processed, "data_dictionary" = Data_Dictionary))
}

p_val_fun <- function(x) {
  SE <- (x[3] - x[1]) / (2 * 1.96)
  z <- x[2] / SE
  P <- exp(-0.717 * z - 0.416 * z^2)
  return(P)
}

set_quantiles <- function(data, X, target, target_q, nontarget_q, subcategory_flag = FALSE) {
  if (subcategory_flag == FALSE) {
    if (is.null(nontarget_q)) {
      data[[target]] <- quantile(data[[target]], target_q)
    } else {
      for (i in X) {
        if (i == target) {
          data[[i]] <- quantile(data[[i]], target_q)
        } else {
          data[[i]] <- quantile(data[[i]], nontarget_q)
        }
      }
    }
  } else {
    for (i in colnames(data)) {
      if (i %in% target) {
        data[[i]] <- quantile(data[[i]], target_q)
      } else {
        data[[i]] <- quantile(data[[i]], nontarget_q)
      }
    }
  }

  return(data)
}

load_model <- function(fit,
                       loss,
                       covars,
                       outcome,
                       data = covid_data_processed,
                       Data_Dictionary = Data_Dictionary) {

  task <- fit$training_task
  dat <- task$data
  X <- task$nodes$covariates
  Y <- task$Y
  preds <- fit$predict_fold(task, fold_number = "validation")
  risk <- mean(loss(preds, Y))
  risk_rescaled <- mean(sqrt((Y - preds)^2) * dat$Population)
  total <- sum(dat[[outcome]] * dat$Population)


  Data_Dictionary_Used <- Data_Dictionary %>%
    filter(Keep == "Yes") %>%
    select(`Variable Name`, `Label`)
  Data_Dictionary_Used <- Data_Dictionary_Used[Data_Dictionary_Used$`Variable Name` %in% covars, ]
  subcategories <- Data_Dictionary_Used$Label
  variable_list <- Data_Dictionary_Used$`Variable Name`

  return_list <- list("Task" = task,
                      "Data" = data,
                      "X" = X,
                      "Y" = Y,
                      "risk" = risk,
                      "total" = total,
                      "Data_Dictionary_Used" = Data_Dictionary_Used,
                      "Subcategories" = subcategories,
                      "Variable_list" = variable_list)
}

var_imp_risk <- function(X, data, outcome, covars, fit, loss, Y, num_boot, Data_Dictionary) {

  ##########################################################################################
  ######################## VARIABLE IMPORTANCE BASED ON RISK ###############################
  ##########################################################################################

  risk_importance <- furrr::future_map_dfr(X, function(i) {
    boot_results_list <- list()
    for (boot in seq(num_boot)) {
      nr <- nrow(data)

      resampled_data <- data[sample(1:nr, size = nr, replace = TRUE), ]
      resampled_data_perm <- resampled_data
      resampled_data_perm[[i]] <- sample(resampled_data[[i]], size = nr)

      task_no_perm <- make_sl3_Task(
        data = resampled_data,
        outcome = outcome,
        covariates = covars
      )

      task_perm <- make_sl3_Task(
        data = resampled_data_perm,
        outcome = outcome,
        covariates = covars
      )

      resampled_sl_preds <- fit$predict_fold(task_no_perm, fold_number = "validation")
      resampled_perm_sl_preds <- fit$predict_fold(task_perm, fold_number = "validation")

      varimp_metric <- mean(loss(resampled_perm_sl_preds, Y)) / mean(loss(resampled_sl_preds, Y))

      boot_results_list[[boot]] <- varimp_metric
    }

    pval <- (1 + sum(unlist(boot_results_list) <= 1)) / (num_boot + 1)
    quantiles <- quantile(unlist(boot_results_list), probs <- c(0.025, 0.50, 0.975))

    results_list <- list("Variable" = i,
                         "Lower_CI" = quantiles[[1]],
                         "Est" = quantiles[[2]],
                         "Upper_CI" = quantiles[[3]],
                         "P_Value" = pval)

    return(results_list)
  }, .options = furrr::furrr_options(seed = TRUE))

  risk_importance$Variable <- Data_Dictionary$`Nice Label`[match(risk_importance$Variable, Data_Dictionary$`Variable Name`)]

  return(risk_importance)
}

subcat_imp_risk <- function(subcategories, data,
                            outcome, covars,
                            fit, loss,
                            Y, num_boot,
                            variable_list) {

  ##########################################################################################
  ######################## SUBGROUP IMPORTANCE BASED ON RISK ###############################
  ##########################################################################################

  X <- unique(subcategories)
  X <- X[X != "outcome"]

  subgroup_risk_importance <- furrr::future_map_dfr(X, function(i) {
    boot_results_list <- list()
    for (boot in seq(num_boot)) {
      nr <- nrow(data)
      subcat_vars <- variable_list[which(subcategories %in% i)]

      resampled_data_perm <- as.data.frame(data[sample(1:nr, size = nr, replace = TRUE), ])
      resampled_data_no_perm <- as.data.frame(data[sample(1:nr, size = nr, replace = TRUE), ])
      resampled_data <- as.data.frame(data[sample(1:nr, size = nr, replace = TRUE), ])

      resampled_data_perm[, subcat_vars] <- resampled_data_no_perm[, subcat_vars]

      task_no_perm <- make_sl3_Task(
        data = resampled_data,
        outcome = outcome,
        covariates = covars
      )

      task_perm <- make_sl3_Task(
        data = resampled_data_perm,
        outcome = outcome,
        covariates = covars
      )

      resampled_sl_preds <- fit$predict_fold(task_no_perm, fold_number = "validation")
      resampled_perm_sl_preds <- fit$predict_fold(task_perm, fold_number = "validation")

      varimp_metric <- mean(loss(resampled_perm_sl_preds, Y)) / mean(loss(resampled_sl_preds, Y))

      boot_results_list[[boot]] <- varimp_metric
    }

    pval <- (1 + sum(unlist(boot_results_list) <= 1)) / (num_boot + 1)
    quantiles <- quantile(unlist(boot_results_list), probs <- c(0.025, 0.50, 0.975))

    results_list <- list("Variable" = i, "Lower_CI" = quantiles[[1]], "Est" = quantiles[[2]], "Upper_CI" = quantiles[[3]], "P_Value" = pval)

    return(results_list)
  }, .options = furrr::furrr_options(seed = TRUE))

  return(subgroup_risk_importance)
}

var_imp_quantile <- function(X, data,
                            outcome, covars,
                            fit, loss,
                            Y, num_boot,
                            variable_list,
                            total, Data_Dictionary, p_val_fun) {


  ##############################################################################
  ######################## QUANTILE INTERACTIONS ###############################
  ##############################################################################

  quantile_importance <- furrr::future_map_dfr(X, function(i) {
    quantile_boot_results_list <- list()
    for (boot in seq(num_boot)) {
      nr <- nrow(data)

      resampled_data <- as.data.frame(data[sample(1:nr, size = nr, replace = TRUE), ])

      # all non-target at 25%
      target_25_nontarget_25 <- set_quantiles(data = resampled_data, X, target = i, target_q = 0.25, nontarget_q = 0.25)
      target_75_nontarget_25 <- set_quantiles(data = resampled_data, X, target = i, target_q = 0.75, nontarget_q = 0.25)

      # all non-target at 75%
      target_25_nontarget_75 <- set_quantiles(data = resampled_data, X, target = i, target_q = 0.25, nontarget_q = 0.75)
      target_75_nontarget_75 <- set_quantiles(data = resampled_data, X, target = i, target_q = 0.75, nontarget_q = 0.75)

      # all non-target at 75%
      target_25_nontarget_50 <- set_quantiles(data = resampled_data, X, target = i, target_q = 0.25, nontarget_q = 0.50)
      target_75_nontarget_50 <- set_quantiles(data = resampled_data, X, target = i, target_q = 0.75, nontarget_q = 0.50)

      task_target_25_nontarget_25 <- make_sl3_Task(
        data = target_25_nontarget_25,
        covariates = covars,
        outcome = outcome
      )

      task_target_75_nontarget_25 <- make_sl3_Task(
        data = target_75_nontarget_25,
        covariates = covars,
        outcome = outcome
      )

      # all non-target at 75%

      task_target_25_nontarget_75 <- make_sl3_Task(
        data = target_25_nontarget_75,
        covariates = covars,
        outcome = outcome
      )

      task_target_75_nontarget_75 <- make_sl3_Task(
        data = target_75_nontarget_75,
        covariates = covars,
        outcome = outcome
      )

      # all non-target at 50%

      task_target_25_nontarget_50 <- make_sl3_Task(
        data = target_25_nontarget_50,
        covariates = covars,
        outcome = outcome
      )

      task_target_75_nontarget_50 <- make_sl3_Task(
        data = target_75_nontarget_50,
        covariates = covars,
        outcome = outcome
      )

      y_25_25_predictions <- fit$predict_fold(task = task_target_25_nontarget_25, fold_number = "validation")
      y_75_25_predictions <- fit$predict_fold(task = task_target_75_nontarget_25, fold_number = "validation")

      y_25_75_predictions <- fit$predict_fold(task = task_target_25_nontarget_75, fold_number = "validation")
      y_75_75_predictions <- fit$predict_fold(task = task_target_75_nontarget_75, fold_number = "validation")

      y_25_50_predictions <- fit$predict_fold(task = task_target_25_nontarget_50, fold_number = "validation")
      y_75_50_predictions <- fit$predict_fold(task = task_target_75_nontarget_50, fold_number = "validation")

      # y_25_obs_predictions <- fit$predict_fold(task = task_target_25_nontarget_obs, fold_number = "validation")
      # y_75_obs_predictions <- fit$predict_fold(task = task_target_75_nontarget_obs, fold_number = "validation")

      delta_rest_high <- mean(y_75_75_predictions - y_25_75_predictions)
      delta_rest_medium <- mean(y_75_50_predictions - y_25_50_predictions)
      delta_rest_low <- mean(y_75_25_predictions - y_25_25_predictions)

      # delta_rest_obs <- mean(y_75_obs_predictions - y_25_obs_predictions)
      varimp_metric <- delta_rest_high - delta_rest_low

      results_list <- list("Delta_High" = delta_rest_high, "Delta_Low" = delta_rest_low, "Delta_Medium" = delta_rest_medium, "Interaction" = varimp_metric)
      quantile_boot_results_list[[boot]] <- results_list
    }

    quantile_boot_results <- bind_rows(quantile_boot_results_list)

    quantile_boot_results_CIs <- t(sapply(quantile_boot_results, quantile, probs <- c(0.025, 0.50, 0.975)))
    pvals <- as.data.frame(apply(quantile_boot_results_CIs, 1, p_val_fun))
    result <- bind_cols(i, quantile_boot_results_CIs, pvals)
    result$Condition <- rownames(result)

    rownames(result) <- NULL
    return(result)
  }, .options = furrr::furrr_options(seed = TRUE))

  colnames(quantile_importance) <- c("Variable", "Lower_CI", "Est", "Upper_CI", "P_Value", "Condition")
  quantile_importance$Variable <- Data_Dictionary$`Nice Label`[match(quantile_importance$Variable, Data_Dictionary$`Variable Name`)]

  quantile_importance[, 2:4] <- as.data.frame(sapply(quantile_importance[, 2:4], as.numeric) * total)
 return(quantile_importance)
}

subcat_imp_quantile <- function(subcategories,
                                data,
                                outcome,
                                covars,
                                fit,
                                loss,
                                Y,
                                num_boot,
                                variable_list,
                                total,
                                Data_Dictionary,
                                p_val_fun) {

  #####################################################################################
  ######################## SUBGROUP QUANTILE IMPORTANCE ###############################
  #####################################################################################

  X <- unique(subcategories)
  X <- X[X != "outcome"]

  subgroup_quantile_importance <- furrr::future_map_dfr(X, function(i) {
    quantile_subcat_boot_results_list <- list()
    for (boot in seq(num_boot)) {
      nr <- nrow(data)

      resampled_data <- as.data.frame(data[sample(1:nr, size = nr, replace = TRUE), ])

      subcat_vars <- variable_list[which(subcategories %in% i)]

      subcat_25_nontarget_obs <- set_quantiles(data = resampled_data, X, target = subcat_vars, target_q = 0.25, nontarget_q = 0.5, subcategory_flag = TRUE)
      subcat_75_nontarget_obs <- set_quantiles(data = resampled_data, X, target = subcat_vars, target_q = 0.75, nontarget_q = 0.5, subcategory_flag = TRUE)

      task_sub_cat_25_nontarget_obs <- make_sl3_Task(
        data = subcat_25_nontarget_obs,
        covariates = covars,
        outcome = outcome
      )

      task_sub_cat_75_nontarget_obs <- make_sl3_Task(
        data = subcat_75_nontarget_obs,
        covariates = covars,
        outcome = outcome
      )

      subcat_75_obs_predictions <- fit$predict_fold(task = task_sub_cat_75_nontarget_obs, fold_number = "validation")
      subcat_25_obs_predictions <- fit$predict_fold(task = task_sub_cat_25_nontarget_obs, fold_number = "validation")

      varimp_metric <- mean(subcat_75_obs_predictions - subcat_25_obs_predictions)

      quantile_subcat_boot_results_list[[boot]] <- varimp_metric
    }

    quantiles <- quantile(unlist(quantile_subcat_boot_results_list), probs <- c(0.025, 0.50, 0.975))
    p_val <- p_val_fun(quantiles)

    results_list <- list("Variable" = i, "Lower_CI" = quantiles[[1]], "Est" = quantiles[[2]], "Upper_CI" = quantiles[[3]], "P_Value" = p_val)

    return(results_list)
  }, .options = furrr::furrr_options(seed = TRUE))

  subgroup_quantile_importance[, c(2:4)] <- as.data.frame(sapply(subgroup_quantile_importance[, c(2:4)], as.numeric) * total)

  return(subgroup_quantile_importance)
}


mips_imp_risk <- function(risk_importance,
                                data,
                                outcome,
                                covars,
                                fit,
                                loss,
                                Y,
                                num_boot,
                                m,
                                Data_Dictionary,
                                p_val_fun) {

  ##############################################################################
  ######################## JOINT PERM INTERACTIONS #############################
  ##############################################################################

  variable_combinations <- combn(subset(risk_importance, risk_importance$`Lower_CI` > 1.0)$Variable, m = m)
  ### Create list with all intxn_size interactions for the intxn_list variable set of interest:
  variable_combinations <- as.data.frame(variable_combinations)
  ### Run the additive vs. joint error calculation for each set of possible interactions of selected size:

  X <- as.vector(apply(variable_combinations, 2, paste, collapse = " & "))

  permuted_importance <- furrr::future_map_dfr(X, function(i) {
    target_vars <- unlist(str_split(i, " & "))
    mips_boot_results_list <- list()
    for (boot in seq(num_boot)) {
      nr <- nrow(dat)
      resampled_data_perm <- as.data.frame(dat[sample(1:nr, size = nr, replace = TRUE), ])
      resampled_data_no_perm <- as.data.frame(dat[sample(1:nr, size = nr, replace = TRUE), ])

      resampled_data_perm[, target_vars] <- resampled_data_no_perm[, target_vars]

      ## compute the additive risk for this set of variables
      additives <- risk_importance %>%
        filter(Variable %in% target_vars) %>%
        select(Est)

      additive_risk <- sum(unlist(additives) - 1) + 1

      task_perm <- make_sl3_Task(
        data = resampled_data_perm,
        outcome = outcome,
        covariates = covars
      )

      # resampled_sl_preds <- fit$predict_fold(task_no_perm, fold_number = "validation")
      resampled_perm_sl_preds <- fit$predict_fold(task_perm, fold_number = "validation")

      risk_scrambled <- mean(loss(resampled_perm_sl_preds, Y))
      varimp_metric <- risk_scrambled / risk

      target_vars_nice <- Data_Dictionary$`Nice Label`[match(target_vars, Data_Dictionary$`Variable Name`)]

      results_list <- list(
        "Variable_Combo" = paste(target_vars_nice, collapse = " & "),
        "Joint_Risk" = varimp_metric,
        "Additive_Risk" = additive_risk
      )

      mips_boot_results_list[[boot]] <- results_list
    }

    mips_boot_results <- bind_rows(mips_boot_results_list)
    mips_boot_results$MIP <- mips_boot_results$Joint_Risk - mips_boot_results$Additive_Risk

    MIPS_CI <- as.data.frame(t(quantile(mips_boot_results$MIP, probs <- c(0.025, 0.50, 0.975))))
    MIPS_pval <- apply(MIPS_CI, 1, p_val_fun)

    mips_result <- list(
      "Variable_Combo" = paste(target_vars_nice, collapse = " & "),
      "Lower_CI" = MIPS_CI[[1]],
      "Est" = MIPS_CI[[2]],
      "Upper_CI" = MIPS_CI[[3]],
      "P_Val" = MIPS_pval
    )

    return(mips_result)
  }, .options = furrr::furrr_options(seed = TRUE))

  return(permuted_importance)
}

