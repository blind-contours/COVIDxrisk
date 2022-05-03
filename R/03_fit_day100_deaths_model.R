library(here)
source(here("R/util.R"))
plan(multisession)
cpus <- 20

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
    for (j in target) {
      data[[j]] <- quantile(data[[j]], target_q)
    }
  }

  return(data)
}

run_varimp <- function(fit,
                       loss,
                       covars,
                       outcome,
                       data = covid_data_processed,
                       Data_Dictionary = Data_Dictionary,
                       label = label,
                       num_boot = num_boot,
                       m = 3) {
  task <- fit$training_task
  dat <- task$data
  X <- task$nodes$covariates
  Y <- task$Y
  preds <- fit$predict_fold(task, fold_number = "validation")
  risk <- mean(loss(preds, Y))

  nworkers <- cpus
  doParallel::registerDoParallel(nworkers)

  Data_Dictionary_Used <- Data_Dictionary %>%
    filter(Keep == "Yes") %>%
    select(`Variable Name`, `Label`)
  Data_Dictionary_Used <- Data_Dictionary_Used[Data_Dictionary_Used$`Variable Name` %in% covars, ]
  subcategories <- Data_Dictionary_Used$Label
  variable_list <- Data_Dictionary_Used$`Variable Name`

  #################################################################################
  ######################## VARIABLE IMPORTANCE RISK ###############################
  #################################################################################

  remaining <- X
  risk_importance_list <- list()
  iter <- 1

  while (length(remaining) > 0) {
    risk_importance <- foreach(i = remaining, .combine = "rbind", .errorhandling = "pass") %dopar% {
      boot_results_list <- list()
      for (boot in seq(num_boot)) {
        nr <- nrow(dat)

        resampled_data <- dat[sample(1:nr, size = nr, replace = TRUE), ]
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

      result <- cbind(i, t(as.data.frame(quantile(unlist(boot_results_list), probs = c(0.25, 0.5, 0.75)))), pval)
      rownames(result) <- NULL
      return(result)
    }
    risk_importance_list[[iter]] <- risk_importance
    remaining <- remaining[remaining %notin% risk_importance[, 1]]
    print(remaining)
    iter <- iter + 1
    print(iter)
  }

  risk_importance <- do.call(rbind, risk_importance_list)
  risk_results <- data.frame(risk_importance)
  # X = names(risk_importance), risk_ratio = unlist(risk_importance))
  colnames(risk_results) <- c("X", "Lower CI", "Est", "Upper CI", "P-Value")
  risk_results[, 2:5] <- sapply(risk_results[, 2:5], as.numeric)
  risk_results_ordered <- risk_results[order(-risk_results$Est), ]

  print("Finished LOO-Risk Importance")

  ##############################################################################
  ######################## SUBGROUP INTERACTIONS ###############################
  ##############################################################################

  remaining <- unique(subcategories)
  remaining <- remaining[remaining != "outcome"]
  subgroup_importance_list <- list()
  iter <- 1

  while (length(remaining) > 0) {
    subgroup_risk_importance <- foreach(i = remaining, .combine = "rbind", .errorhandling = "pass") %dopar% {
      subcat_boot_results_list <- list()
      for (boot in seq(num_boot)) {
        nr <- nrow(dat)
        subcat_vars <- variable_list[which(subcategories %in% i)]

        resampled_data_perm <- as.data.frame(dat[sample(1:nr, size = nr, replace = TRUE), ])
        resampled_data_no_perm <- as.data.frame(dat[sample(1:nr, size = nr, replace = TRUE), ])
        resampled_data <- as.data.frame(dat[sample(1:nr, size = nr, replace = TRUE), ])

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

        subcat_boot_results_list[[boot]] <- varimp_metric
      }

      pval <- (1 + sum(unlist(subcat_boot_results_list) <= 1)) / (num_boot + 1)

      result <- cbind(i, t(as.data.frame(quantile(unlist(subcat_boot_results_list), probs = c(0.25, 0.5, 0.75)))), pval)
      rownames(result) <- NULL
      return(result)
    }

    subgroup_importance_list[[iter]] <- subgroup_risk_importance
    remaining <- remaining[remaining %notin% subgroup_risk_importance[, 1]]
    print(remaining)
    iter <- iter + 1
    print(iter)
  }

  subgroup_importance <- as.data.frame(do.call(rbind, subgroup_importance_list))

  colnames(subgroup_importance) <- c("X", "Lower CI", "Est", "Upper CI", "P-Value")
  subgroup_importance[, 2:5] <- sapply(subgroup_importance[, 2:5], as.numeric)
  subgroup_importance_ordered <- subgroup_importance[order(-subgroup_importance$Est), ]

  ##############################################################################
  ######################## QUANTILE INTERACTIONS ###############################
  ##############################################################################

  remaining <- X
  quantile_importance_list <- list()
  iter <- 1

  while (length(remaining) > 0) {
    quantile_importance <- foreach(i = remaining, .combine = "rbind", .errorhandling = "pass") %dopar% {
      quantile_boot_results_list <- list()
      for (boot in seq(num_boot)) {
        dat <- fit$training_task$data
        nr <- nrow(dat)

        resampled_data <- as.data.frame(dat[sample(1:nr, size = nr, replace = TRUE), ])

        # all non-target at 25%
        target_25_nontarget_25 <- set_quantiles(data = resampled_data, X, target = i, target_q = 0.25, nontarget_q = 0.25)
        target_75_nontarget_25 <- set_quantiles(data = resampled_data, X, target = i, target_q = 0.75, nontarget_q = 0.25)

        # all non-target at 75%
        target_25_nontarget_75 <- set_quantiles(data = resampled_data, X, target = i, target_q = 0.25, nontarget_q = 0.75)
        target_75_nontarget_75 <- set_quantiles(data = resampled_data, X, target = i, target_q = 0.75, nontarget_q = 0.75)

        # all non-target at 75%
        target_25_nontarget_50 <- set_quantiles(data = resampled_data, X, target = i, target_q = 0.25, nontarget_q = 0.50)
        target_75_nontarget_50 <- set_quantiles(data = resampled_data, X, target = i, target_q = 0.75, nontarget_q = 0.50)

        # all non-target at observed
        # target_25_nontarget_obs <-set_quantiles(data = dat, X, target = i, target_q = 0.25, nontarget_q = NULL)
        # target_75_nontarget_obs <-set_quantiles(data = dat, X, target = i, target_q = 0.75, nontarget_q = NULL)


        # all non-target at 25%

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

        result <- cbind(i, delta_rest_high, delta_rest_medium, delta_rest_low, varimp_metric)
        quantile_boot_results_list[[boot]] <- result
      }

      quantile_boot_results <- as.data.frame(do.call(rbind, quantile_boot_results_list))
      quantile_boot_results[, 2:5] <- sapply(quantile_boot_results[, 2:5], as.numeric)

      quantile_boot_results_CIs <- t(sapply(quantile_boot_results[, 2:5], quantile, probs = c(0.25, 0.5, 0.75)))
      result <- cbind(i, quantile_boot_results_CIs, as.data.frame(apply(quantile_boot_results_CIs, 1, p_val_fun)))
      result$condition <- rownames(result)
      rownames(result) <- NULL
      return(result)
    }

    quantile_importance_list[[iter]] <- quantile_importance
    remaining <- remaining[remaining %notin% quantile_importance[, 1]]
    print(remaining)
    iter <- iter + 1
    print(iter)
  }

  quantile_importance <- do.call(rbind, quantile_importance_list)
  quantile_results <- data.table(quantile_importance)
  colnames(quantile_results) <- c("X", "Lower CI", "Est", "Upper CI", "P-Value", "Condition")
  quantile_results_ordered <- quantile_results[order(-quantile_results$Est)]

  print("Finished Quantile Interaction Search")

  #####################################################################################
  ######################## SUBGROUP QUANTILE IMPORTANCE ###############################
  #####################################################################################

  remaining <- unique(subcategories)
  remaining <- remaining[remaining != "outcome"]
  subgroup_quantile_importance_list <- list()
  iter <- 1

  while (length(remaining) > 0) {
    subgroup_quantile_importance <- foreach(i = remaining, .combine = "rbind", .errorhandling = "pass") %dopar% {
      quantile_subcat_boot_results_list <- list()
      for (boot in seq(num_boot)) {
        dat <- fit$training_task$data
        nr <- nrow(dat)

        resampled_data <- as.data.frame(dat[sample(1:nr, size = nr, replace = TRUE), ])

        subcat_vars <- variable_list[which(subcategories %in% i)]

        subcat_25_nontarget_obs <- set_quantiles(data = resampled_data, X, target = subcat_vars, target_q = 0.25, nontarget_q = NULL, subcategory_flag = TRUE)
        subcat_75_nontarget_obs <- set_quantiles(data = resampled_data, X, target = subcat_vars, target_q = 0.75, nontarget_q = NULL, subcategory_flag = TRUE)

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
        result <- cbind(i, varimp_metric)
        quantile_subcat_boot_results_list[[boot]] <- result
      }

      quantile_subcat_boot_results <- as.data.frame(do.call(rbind, quantile_subcat_boot_results_list))
      quantile_subcat_boot_results$varimp_metric <- as.numeric(quantile_subcat_boot_results$varimp_metric)

      quantile_subcat_boot_CIs <- t(as.data.frame(quantile(quantile_subcat_boot_results$varimp_metric, probs = c(0.25, 0.5, 0.75))))
      pvals <- apply(quantile_subcat_boot_CIs, 1, p_val_fun)

      result <- cbind(i, quantile_subcat_boot_CIs, pvals)
      colnames(result) <- c("X", "Lower CI", "Est", "Upper CI", "P-vals")
      rownames(result) <- NULL
      return(result)
    }

    subgroup_quantile_importance_list[[iter]] <- subgroup_quantile_importance
    remaining <- remaining[remaining %notin% subgroup_quantile_importance[, 1]]
    print(remaining)
    iter <- iter + 1
    print(iter)
  }

  subgroup_quantile_importance <- do.call(rbind, subgroup_quantile_importance_list)
  subgroup_quantile_importance <- data.table(subgroup_quantile_importance)
  # colnames(subgroup_quantile_importance) <- c("X", "Sub-Category Delta")
  subgroup_quantile_results_ordered <- subgroup_quantile_importance[order(-subgroup_quantile_importance$Est)]

  total <- sum(dat[[outcome]] * dat$Population)

  subgroup_quantile_results_ordered[, c(2:4)] <- as.data.frame(sapply(subgroup_quantile_results_ordered[, c(2:4)], as.numeric) * total)

  risk_results_ordered$X <- Data_Dictionary$`Nice Label`[match(risk_results_ordered$X, Data_Dictionary$`Variable Name`)]
  quantile_results_ordered$X <- Data_Dictionary$`Nice Label`[match(quantile_results_ordered$X, Data_Dictionary$`Variable Name`)]

  quantile_results_ordered[, 2:4] <- as.data.frame(sapply(quantile_results_ordered[, 2:4], as.numeric) * total)

  risk_rescaled <- mean(sqrt((Y - preds)^2) * dat$Population)

  variable_combinations <- combn(subset(risk_results, risk_results$`Lower CI` > 1)$X, m = m)
  ### Create list with all intxn_size interactions for the intxn_list variable set of interest:
  variable_combinations <- as.data.frame(variable_combinations)
  ### Run the additive vs. joint error calculation for each set of possible interactions of selected size:

  ##############################################################################
  ######################## JOINT PERM INTERACTIONS #############################
  ##############################################################################

  joint_perm_importance_list <- list()
  iter <- 1
  remaining <- as.vector(apply(variable_combinations, 2, paste, collapse = " & "))

  while (length(remaining) > 0) {
    permuted_importance <- foreach(i = 1:length(remaining), .combine = "rbind") %dopar% {
      target_vars <- remaining[i]
      target_vars <- unlist(str_split(target_vars, " & "))
      mips_boot_results_list <- list()
      for (boot in seq(num_boot)) {
        nr <- nrow(dat)
        resampled_data_perm <- as.data.frame(dat[sample(1:nr, size = nr, replace = TRUE), ])
        resampled_data_no_perm <- as.data.frame(dat[sample(1:nr, size = nr, replace = TRUE), ])

        resampled_data_perm[, target_vars] <- resampled_data_no_perm[, target_vars]

        ## compute the additive risk for this set of variables
        additives <- risk_results %>%
          filter(X %in% target_vars) %>%
          select(Est)
        additive_risk <- sum(unlist(additives) - 1) + 1

        # task_no_perm <- make_sl3_Task(data = resampled_data_no_perm,
        #                               outcome = outcome,
        #                               covariates = covars)

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

        result <- cbind(paste(target_vars_nice, collapse = " & "), varimp_metric, additive_risk, remaining[[i]])
        mips_boot_results_list[[boot]] <- result
      }

      mips_boot_results <- as.data.frame(do.call(rbind, mips_boot_results_list))
      mips_boot_results[, 2:3] <- sapply(mips_boot_results[, 2:3], as.numeric)
      mips_boot_results$MIP <- mips_boot_results$varimp_metric - mips_boot_results$additive_risk

      MIPS_CI <- as.data.frame(t(quantile(mips_boot_results$MIP, probs = c(0.25, 0.5, 0.75))))
      MIPS_pval <- apply(MIPS_CI, 1, p_val_fun)

      mips_result <- cbind(paste(target_vars_nice, collapse = " & "), MIPS_CI, MIPS_pval, remaining[[i]])

      colnames(mips_result) <- c("Nice Label", "Lower CI", "Est", "Upper CI", "P-Value", "Raw Vars")
      return(mips_result)
    }

    joint_perm_importance_list[[iter]] <- permuted_importance
    remaining <- remaining[remaining %notin% permuted_importance[, 6]]
    print(remaining)
    iter <- iter + 1
    print(iter)
  }

  permuted_importance <- do.call(rbind, joint_perm_importance_list)

  print("Finished Joint Permutation")

  permuted_importance <- as.data.frame(permuted_importance)
  permuted_importance <- permuted_importance[1:5]

  permuted_importance <- permuted_importance[order(permuted_importance$Est, decreasing = TRUE), ]

  permuted_importance <- subset(permuted_importance, permuted_importance$`Lower CI` > 0)

  ## Individual Variable Importance From Model Risk Ratio:

  risk_plot <- risk_results_ordered %>%
    arrange(Est) %>%
    # First sort by val. This sort the dataframe but NOT the factor levels
    filter(Est > 1.0) %>%
    mutate(name = factor(X, levels = X)) %>%
    # This trick update the factor levels
    ggplot(aes(x = name, y = Est)) +
    geom_point(size = 4, color = "orange") +
    geom_errorbar(aes(ymin = `Lower CI`, ymax = `Upper CI`), width = .1) +
    coord_flip() +
    theme_bw(base_size = 12) +
    ylab("Model Risk Ratio") +
    xlab("County Feature")

  ## Individual Variable Importance From Model Quantile Predictions:

  quantile_results_nonzero <- quantile_results_ordered[rowSums(quantile_results_ordered[, 2:4]) != 0, ]

  quantile_preds_plot <- quantile_results_nonzero %>%
    filter(Condition != "varimp_metric") %>%
    # mutate(name=factor(X, levels=X)) %>%   # This trick update the factor levels
    ggplot(aes(x = X, y = Est, colour = factor(Condition))) +
    geom_point(size = 2) +
    geom_errorbar(aes(ymin = `Lower CI`, ymax = `Upper CI`), width = .1) +
    coord_flip() +
    theme_bw(base_size = 12) +
    ylab("Quantile-Based Predictions") +
    xlab("County Feature") +
    scale_y_continuous(breaks = scales::pretty_breaks(n = 10))

  ## Interaction Variable Importance From Model Quantile Predictions:

  quantile_intxn_plot <- quantile_results_nonzero %>%
    filter(Condition == "varimp_metric") %>%
    mutate(name = factor(X, levels = X)) %>%
    # This trick update the factor levels
    ggplot(aes(x = name, y = Est)) +
    geom_point(size = 4, color = "blue") +
    geom_errorbar(aes(ymin = `Lower CI`, ymax = `Upper CI`), width = .1) +
    coord_flip() +
    theme_bw(base_size = 12) +
    ylab("Quantile-Based Interactions") +
    xlab("County Feature")

  joint_permutation_plot <- permuted_importance %>%
    ggplot(aes(`Nice Label`, Est)) +
    geom_point(size = 4, color = "chartreuse4") +
    geom_errorbar(aes(ymin = `Lower CI`, ymax = `Upper CI`), width = .1) +
    coord_flip() +
    theme_bw(base_size = 12) +
    ylab("Joint vs. Additive Individual Risk Difference") +
    xlab("Risk Combination")

  sub_group_risk_plot <- subgroup_importance_ordered %>%
    arrange(Est) %>% # First sort by val. This sort the dataframe but NOT the factor levels
    mutate(name = factor(X, levels = X)) %>% # This trick update the factor levels
    ggplot(aes(x = name, y = Est)) +
    geom_errorbar(aes(ymin = `Lower CI`, ymax = `Upper CI`), size = 0.8, width = 0.2) +
    geom_point(size = 4, color = "orange") +
    coord_flip() +
    theme_bw(base_size = 12) +
    ylab("Model Risk Ratio") +
    xlab("County-Level Risk Subgroups")

  sub_group_quantile_plot <- subgroup_quantile_results_ordered %>%
    arrange(Est) %>%
    mutate(name = factor(X, levels = X)) %>%
    ggplot(aes(x = name, y = `Est`)) +
    geom_errorbar(aes(ymin = `Lower CI`, ymax = `Upper CI`), size = 0.8, width = 0.2) +
    geom_point(size = 3, color = "red") +
    coord_flip() +
    theme_bw(base_size = 12) +
    ylab("All Sub-Group 75th - 25th Quantile Predictions") +
    xlab("County Feature")

  # p <- plot_grid(risk_plot, quantiles_plot, interaction_plot, labels = label, vjust = -0.1, nrow = 1)
  ggsave(here(paste("Figures/", "varimp_", label, ".png", sep = "")), risk_plot, width = 8, height = 6)
  ggsave(here(paste("Figures/", "quantile_imp_", label, ".png", sep = "")), quantile_preds_plot, width = 8, height = 6)
  ggsave(here(paste("Figures/", "interaction_imp_", label, ".png", sep = "")), quantile_intxn_plot, width = 8, height = 6)
  ggsave(here(paste("Figures/", "jointimp_", label, ".png", sep = "")), joint_permutation_plot, width = 14, height = 6)
  ggsave(here(paste("Figures/", "subgroup_imp_", label, ".png", sep = "")), sub_group_risk_plot, width = 8, height = 6)
  ggsave(here(paste("Figures/", "subgroup_quantile_imp_", label, ".png", sep = "")), sub_group_quantile_plot, width = 8, height = 6)


  return(list(
    "Var_Risk_Results" = risk_results_ordered,
    "Subgroup_Risk_Results" = subgroup_importance_ordered,
    "Var_Quantile_Results" = quantile_results_ordered,
    "Subgroup_Quantile_Results" = subgroup_quantile_results_ordered,
    "Intxn_Risk_Results" = permuted_importance,
    "model_risk" = risk_rescaled
  ))
}

fit_sl_varimp <- function(outcome, label, num_boot) {

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

  ## get variable importance from the sl3 object
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

  SL_results <- list("fit" = sl_fit, "var_imp" = var_importance)

  saveRDS(SL_results, here(paste("Models/", outcome, ".RDS", sep = "")))

  return(NULL)
}

Deathsday100 <- fit_sl_varimp(outcome = "CountyRelativeDay100Deaths", label = "COVID-19 Deaths at Day 100", num_boot = 10)
