create_sl <- function(data = covid_data_processed,
                      outcome = outcome,
                      all_outcomes = all_outcomes) {
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

load_data <- function(path_data = "cleaned_covid_data_final_Mar_31_22.csv",
                      path_data_dict = "Data_Dictionary.xlsx") {
  ## load data
  covid_data_processed <- read.csv(PROCESSED_DATA_PATH(path_data),
    check.names = FALSE
  )

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

  data <- data[1, ]

  return(data)
}

set_cond_quantiles <- function(data, target) {

  quantiles <- quantile(data[[target]])

  q1_group <- data[data[[target]] <= quantiles[2], ]
  q4_group <- data[data[[target]] >= quantiles[4], ]

  pairwise_dist <- dist2(
    x = scale(q1_group[, which(colnames(q1_group) != target)]),
    y = scale(q4_group[, which(colnames(q4_group) != target)])
  )

  min_dist <- min(pairwise_dist)
  min_element_index <- which(pairwise_dist == min_dist, arr.ind = TRUE)

  selected_q1 <- q1_group[min_element_index[, 1], ]
  selected_q4 <- q4_group[min_element_index[, 2], ]

  selected_q1[[target]] <- quantiles[2]
  selected_q4[[target]] <- quantiles[4]

  selected_q1 <- selected_q1[1, ]
  selected_q4 <- selected_q4[1, ]

  return(list("Q4_data" = selected_q4, "Q1_data" = selected_q1, "Dist" = min_dist))
}

set_cond_pair_quantiles <- function(data, targets, cor_dir_ind) {

  # rownames(data) <- data$fips
  data_scaled <- as.data.frame(scale(data))

  quantiles_t1 <- quantile(data_scaled[[targets[1]]])
  quantiles_t2 <- quantile(data_scaled[[targets[2]]])

  if (cor_dir_ind == 1) {
    q11_group_ind <- data_scaled[[targets[1]]] <= quantiles_t1[2] & data_scaled[[targets[2]]] <= quantiles_t2[2]
    q44_group_ind <- data_scaled[[targets[1]]] >= quantiles_t1[4] & data_scaled[[targets[2]]] >= quantiles_t2[4]

    no_support_flag <- any(length(unique(q11_group_ind)) == 1 | length(unique(q44_group_ind)) == 1)

    if (no_support_flag == TRUE) {
      return(NULL)
    } else {
      q11_group <- data_scaled[q11_group_ind, ]
      q44_group <- data_scaled[q44_group_ind, ]
    }

    pairwise_dist <- dist2(
      x = q11_group[, which(colnames(q11_group) %notin% targets)],
      y = q44_group[, which(colnames(q44_group) %notin% targets)]
    )

    min_dist <- min(pairwise_dist)
    min_element_index <- which(pairwise_dist == min_dist, arr.ind = TRUE)

    if (dim(min_element_index)[2] == 1) {
      selected_q11_obs <- q11_group[min_element_index[, 1], ]
      selected_q44_obs <- q44_group[min_element_index[, 2], ]
    } else {
      selected_q11_obs <- q11_group[min_element_index[1, 1], ]
      selected_q44_obs <- q44_group[min_element_index[1, 2], ]
    }

    selected_q11_obs[[targets[1]]] <- quantiles_t1[2]
    selected_q11_obs[[targets[2]]] <- quantiles_t2[2]

    selected_q44_obs[[targets[1]]] <- quantiles_t1[4]
    selected_q44_obs[[targets[2]]] <- quantiles_t2[4]

    return(list("Data_1" = selected_q44_obs, "Data_2" = selected_q11_obs))
  } else {
    q14_group_ind <- data_scaled[[targets[1]]] <= quantiles_t1[2] & data_scaled[[targets[2]]] >= quantiles_t2[4]
    q41_group_ind <- data_scaled[[targets[1]]] >= quantiles_t1[4] & data_scaled[[targets[2]]] <= quantiles_t2[2]

    no_support_flag <- any(length(unique(q14_group_ind)) == 1 | length(unique(q41_group_ind)) == 1)

    if (no_support_flag == TRUE) {
      return(NULL)
    } else {
      q14_group <- data_scaled[q14_group_ind, ]
      q41_group <- data_scaled[q41_group_ind, ]
    }

    pairwise_dist <- dist2(
      x = q14_group[, which(colnames(q14_group) %notin% targets)],
      y = q41_group[, which(colnames(q41_group) %notin% targets)]
    )

    min_dist <- min(pairwise_dist)
    min_element_index <- which(pairwise_dist == min_dist, arr.ind = TRUE)

    if (dim(min_element_index)[2] == 1) {
      selected_q14_obs <- q14_group[min_element_index[, 1], ]
      selected_q41_obs <- q41_group[min_element_index[, 2], ]
    } else {
      selected_q14_obs <- q14_group[min_element_index[1, 1], ]
      selected_q41_obs <- q41_group[min_element_index[1, 2], ]
    }

    selected_q14_obs[[targets[1]]] <- quantiles_t1[2]
    selected_q14_obs[[targets[2]]] <- quantiles_t2[4]

    selected_q41_obs[[targets[1]]] <- quantiles_t1[4]
    selected_q41_obs[[targets[2]]] <- quantiles_t2[2]

    return(list("Data_1" = selected_q14_obs, "Data_2" = selected_q41_obs))
  }
}

set_mips_quantiles <- function(data, targets, target_qs, nontarget_q) {
  threshs <- list()
  for (var_i in seq(targets)) {
    target <- targets[var_i]
    thresh_i <- quantile(data[[target]], probs = target_qs[var_i])
    threshs[var_i] <- thresh_i
  }
  threshs <- unlist(threshs)

  for (thresh in seq(threshs)) {
    target <- targets[thresh]
    if (target_qs[thresh] == 0.25) {
      data <- data[data[[target]] <= threshs[thresh], ]
    } else {
      data <- data[data[[target]] >= threshs[thresh], ]
    }
  }

  medians <- as.data.frame(matrix(sapply(data, quantile, probs = nontarget_q), nrow = 1))
  colnames(medians) <- colnames(data)

  for (i in seq(targets)) {
    target <- targets[i]
    medians[[target]] <- threshs[i]
  }
  # medians <- sapply(medians, rep.int, times=10)
  # medians <- as.data.frame(medians)
  return(medians)
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

  return_list <- list(
    "Task" = task,
    "Data" = data,
    "X" = X,
    "Y" = Y,
    "risk" = risk,
    "risk_rescaled" = risk_rescaled,
    "total" = total,
    "Data_Dictionary_Used" = Data_Dictionary_Used,
    "Subcategories" = subcategories,
    "Variable_list" = variable_list
  )
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

    results_list <- list(
      "Variable" = i,
      "Lower_CI" = quantiles[[1]],
      "Est" = quantiles[[2]],
      "Upper_CI" = quantiles[[3]],
      "P_Value" = pval
    )

    return(results_list)
  }, .options = furrr::furrr_options(seed = TRUE))

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

var_imp_quantile <- function(X,
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


  ##############################################################################
  ######################## QUANTILE INTERACTIONS ###############################
  ##############################################################################

  quantile_importance <- furrr::future_map_dfr(X, function(i) {
    quantile_boot_results_list <- list()
    for (boot in seq(num_boot)) {
      nr <- nrow(data)

      resampled_data <- as.data.frame(data[sample(1:nr, size = nr, replace = TRUE), ])

      # all non-target at 25%
      counterfactual_data <- set_cond_quantiles(data = resampled_data, target = i)

      Q1_data <- counterfactual_data$Q1_data
      Q4_data <- counterfactual_data$Q4_data

      Q1_task <- make_sl3_Task(
        data = Q1_data,
        covariates = covars,
        outcome = outcome
      )

      Q4_task <- make_sl3_Task(
        data = Q4_data,
        covariates = covars,
        outcome = outcome
      )

      Q1_predictions <- fit$predict_fold(task = Q1_task, fold_number = "full") * total
      Q4_predictions <- fit$predict_fold(task = Q4_task, fold_number = "full") * total

      delta_Q4_Q1 <- Q4_predictions - Q1_predictions

      results_list <- list("Delta_Q4_Q1" = delta_Q4_Q1, "Ref_Dist" = counterfactual_data$Dist)
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
  quantile_importance$Label <- Data_Dictionary$`Nice Label`[match(quantile_importance$Variable, Data_Dictionary$`Variable Name`)]

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

      task_sub_cat_75_nontarget_obs <- make_sl3_Task(
        data = subcat_75_nontarget_obs,
        covariates = covars,
        outcome = outcome
      )

      task_sub_cat_25_nontarget_obs <- make_sl3_Task(
        data = subcat_25_nontarget_obs,
        covariates = covars,
        outcome = outcome
      )

      subcat_75_obs_predictions <- fit$predict_fold(task = task_sub_cat_75_nontarget_obs, fold_number = "full")
      subcat_25_obs_predictions <- fit$predict_fold(task = task_sub_cat_25_nontarget_obs, fold_number = "full")

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
                          p_val_fun,
                          risk) {

  ##############################################################################
  ######################## JOINT PERM INTERACTIONS #############################
  ##############################################################################

  cut_off <- quantile(risk_importance$Est, 0.75)
  variable_combinations <- combn(subset(risk_importance, risk_importance$Est > cut_off)$Variable, m = m)
  ### Create list with all intxn_size interactions for the intxn_list variable set of interest:
  X <- as.data.frame(variable_combinations)
  ### Run the additive vs. joint error calculation for each set of possible interactions of selected size:

  permuted_importance <- furrr::future_map_dfr(1:dim(X)[2], function(i) {
    target_vars <- X[, i]
    mips_boot_results_list <- list()
    for (boot in seq(num_boot)) {
      nr <- nrow(data)
      resampled_data_perm <- as.data.frame(data[sample(1:nr, size = nr, replace = TRUE), ])
      resampled_data_no_perm <- as.data.frame(data[sample(1:nr, size = nr, replace = TRUE), ])

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

mips_imp_quantile <- function(quantile_importance,
                              data,
                              outcome,
                              covars,
                              fit,
                              loss,
                              Y,
                              num_boot,
                              m,
                              Data_Dictionary,
                              p_val_fun,
                              total) {

  ##############################################################################
  ######################## QUANTILE BASED INTERACTIONS #########################
  ##############################################################################

  quantile_importance <- quantile_importance %>% filter(Condition == "Delta_Q4_Q1")

  cut_off <- quantile(quantile_importance$Est, 0.97)
  variable_combinations <- combn(subset(quantile_importance, quantile_importance$Est > cut_off)$Variable, m = 2)
  ### Create list with all intxn_size interactions for the intxn_list variable set of interest:
  X <- as.data.frame(variable_combinations)
  ### Run the additive vs. joint error calculation for each set of possible interactions of selected size:

  print(dim(X))

  mips_quantile_importance_results <- furrr::future_map_dfr(1:dim(X)[2], function(i) {
    target_vars <- X[, i]
    quantile_mips_boot_results_list <- list()

    additives <- quantile_importance %>%
      filter(Variable %in% target_vars) %>%
      select(Est)

    additive_sum <- sum(additives)

    for (boot in seq(num_boot)) {
      nr <- nrow(data)

      resampled_data <- as.data.frame(data[sample(1:nr, size = nr, replace = TRUE), ])

      corr_value <- cor(resampled_data[target_vars])[1, 2]
      corr_pos_ind <- ifelse(corr_value > 0, 1, 0)

      counterfactual_data <- set_cond_pair_quantiles(data = resampled_data, targets = target_vars, cor_dir_ind = corr_pos_ind)

      if (is.null(counterfactual_data)) {
        target_vars_nice <- Data_Dictionary$`Nice Label`[match(target_vars, Data_Dictionary$`Variable Name`)]
        result <- data.frame(mat = matrix(ncol = 5, nrow = 1))
        result$Variables <- paste(target_vars_nice, collapse = " & ")
        result$Corr <- NA

        return(result)
      } else {
        task_data1 <- make_sl3_Task(
          data = counterfactual_data$Data_1,
          outcome = outcome,
          covariates = covars
        )

        task_data2 <- make_sl3_Task(
          data = counterfactual_data$Data_2,
          outcome = outcome,
          covariates = covars
        )

        data1_sl_preds <- fit$predict_fold(task_data1, fold_number = "full") * total
        data2_sl_preds <- fit$predict_fold(task_data2, fold_number = "full") * total

        if (corr_pos_ind == 1) {
          varimp_metric <- mean(data1_sl_preds - data2_sl_preds)
        } else {
          varimp_metric <- mean(abs(data1_sl_preds - data2_sl_preds))
        }

        joint_additive_diff <- varimp_metric - additive_sum

        results_list <- list("Joint_Diff" = varimp_metric, "Additive_Diff" = additive_sum, "Interaction" = joint_additive_diff)

        quantile_mips_boot_results_list[[boot]] <- results_list
      }

      quantile_mips_boot_results <- bind_rows(quantile_mips_boot_results_list)
      quantile_mips_results_CIs <- t(sapply(quantile_mips_boot_results, quantile, probs <- c(0.025, 0.50, 0.975)))
      pvals <- as.data.frame(apply(quantile_mips_results_CIs, 1, p_val_fun))

      target_vars_nice <- Data_Dictionary$`Nice Label`[match(target_vars, Data_Dictionary$`Variable Name`)]

      result <- bind_cols(quantile_mips_results_CIs, pvals)
      result$Condition <- rownames(result)
      result$Variables <- paste(target_vars_nice, collapse = " & ")

      rownames(result) <- NULL

      result$Corr <- corr_value
    }


    return(result)
  }, .options = furrr::furrr_options(seed = TRUE))

  colnames(mips_quantile_importance_results) <- c("Lower_CI", "Est", "Upper_CI", "P_Value", "Condition", "Variable_Comb", "Correlation")

  return(mips_quantile_importance_results)
}
