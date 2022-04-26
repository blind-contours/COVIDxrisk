library(here)
source(here("R/util.R"))
plan(multisession)
cpus <- 20

set_quantiles <- function(data, X, target, target_q, nontarget_q, subcategory_flag = FALSE){
  if (subcategory_flag == FALSE) {
    if (is.null(nontarget_q)) {
      data[[target]] <- quantile(data[[target]], target_q)
    }else{
      for (i in X) {
        if (i == target) {
          data[[i]] <- quantile(data[[i]], target_q)
        }
        else{
          data[[i]] <- quantile(data[[i]], nontarget_q)
        }
      }
    }

  }else{
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
                       thresh = 1.01,
                       m = 3)
{

  task <- fit$training_task
  dat <- task$data
  X <- task$nodes$covariates
  Y <- task$Y
  preds <- fit$predict_fold(task, fold_number = "validation")
  risk <- mean(loss(preds, Y))

  nworkers <- cpus
  doParallel::registerDoParallel(nworkers)

  #################################################################################
  ######################## VARIABLE IMPORTANCE RISK ###############################
  #################################################################################

  remaining <- X
  risk_importance_list <- list()
  iter <- 1

  while (length(remaining) > 0) {
    risk_importance <- foreach(i = remaining, .combine = 'rbind', .errorhandling = "pass") %dopar% {

      scrambled_col <- data.table(sample(unlist(dat[, i, with = FALSE]), nrow(dat)))
      names(scrambled_col) <- i
      scrambled_col_names <- task$add_columns(scrambled_col)
      scrambled_col_task <- task$next_in_chain(column_names = scrambled_col_names)
      scrambled_sl_preds <- fit$predict_fold(scrambled_col_task, fold_number = "validation")
      no_i_risk <- mean(loss(scrambled_sl_preds, Y))
      varimp_metric <- no_i_risk/risk

      result <- cbind(i, varimp_metric)

      return(result)
    }
    risk_importance_list[[iter]] <- risk_importance
    remaining <- remaining[remaining %notin% risk_importance[,1]]
    print(remaining)
    iter <- iter + 1
    print(iter)
  }

  risk_importance <- do.call(rbind, risk_importance_list)
  risk_results <- data.frame(risk_importance)
  # X = names(risk_importance), risk_ratio = unlist(risk_importance))
  colnames(risk_results) <- c("X", "risk_ratio")
  risk_results$risk_ratio <- as.numeric(risk_results$risk_ratio)
  risk_results_ordered <- risk_results[order(-risk_results$risk_ratio),]

  print("Finished LOO-Risk Importance")

  ##############################################################################
  ######################## SUBGROUP INTERACTIONS ###############################
  ##############################################################################

  Data_Dictionary_Used <- Data_Dictionary %>% filter(Keep == "Yes") %>% select(`Variable Name`, `Label`)
  Data_Dictionary_Used <- Data_Dictionary_Used[Data_Dictionary_Used$`Variable Name` %in% covars ,]
  subcategories <- Data_Dictionary_Used$Label
  variable_list <- Data_Dictionary_Used$`Variable Name`

  remaining <- unique(subcategories)
  remaining <- remaining[remaining!="outcome"]
  subgroup_importance_list <- list()
  iter <- 1

  while (length(remaining) > 0) {
    subgroup_risk_importance <- foreach(i = remaining, .combine = 'rbind', .errorhandling = "pass") %dopar% {

      subcat_vars <- variable_list[which(subcategories %in% i)]
      scrambled_rows <- dat[sample(nrow(dat)), ]
      scrambled_rows_selection <- scrambled_rows %>% dplyr::select(!!subcat_vars)
      scrambled_col_names <- task$add_columns(scrambled_rows_selection)
      scrambled_col_task <- task$next_in_chain(column_names = scrambled_col_names)
      scrambled_sl_preds <- fit$predict_fold(scrambled_col_task,
                                             fold_number = "validation")

      risk_scrambled <- mean(loss(scrambled_sl_preds, Y))
      varimp_metric <- risk_scrambled/risk
      result <- cbind(i, varimp_metric)
      return(result)
    }

    subgroup_importance_list[[iter]] <- subgroup_risk_importance
    remaining <- remaining[remaining %notin% subgroup_risk_importance[,1]]
    print(remaining)
    iter <- iter + 1
    print(iter)
  }

  subgroup_importance <- as.data.frame(do.call(rbind, subgroup_importance_list))

  colnames(subgroup_importance) <- c("X", "subgroup_risk_ratio")
  subgroup_importance$subgroup_risk_ratio <- as.numeric(subgroup_importance$subgroup_risk_ratio)
  subgroup_importance_ordered <- subgroup_importance[order(-subgroup_importance$subgroup_risk_ratio),]


  ##############################################################################
  ######################## QUANTILE INTERACTIONS ###############################
  ##############################################################################

  remaining <- X
  quantile_importance_list <- list()
  iter <- 1

  while (length(remaining) > 0) {
    quantile_importance <- foreach(i = remaining, .combine = 'rbind',  .errorhandling = "pass") %dopar% {

      dat <- fit$training_task$data

      # all non-target at 25%
      target_25_nontarget_25 <-set_quantiles(data = dat, X, target = i, target_q = 0.25, nontarget_q = 0.25)
      target_75_nontarget_25 <-set_quantiles(data = dat, X, target = i, target_q = 0.75, nontarget_q = 0.25)

      # all non-target at 75%
      target_25_nontarget_75 <-set_quantiles(data = dat, X, target = i, target_q = 0.25, nontarget_q = 0.75)
      target_75_nontarget_75 <-set_quantiles(data = dat, X, target = i, target_q = 0.75, nontarget_q = 0.75)

      # all non-target at 75%
      target_25_nontarget_50 <-set_quantiles(data = dat, X, target = i, target_q = 0.25, nontarget_q = 0.50)
      target_75_nontarget_50 <-set_quantiles(data = dat, X, target = i, target_q = 0.75, nontarget_q = 0.50)

      # all non-target at observed
      # target_25_nontarget_obs <-set_quantiles(data = dat, X, target = i, target_q = 0.25, nontarget_q = NULL)
      # target_75_nontarget_obs <-set_quantiles(data = dat, X, target = i, target_q = 0.75, nontarget_q = NULL)


      # all non-target at 25%

      task_target_25_nontarget_25 <- make_sl3_Task(
        data = target_25_nontarget_25,
        covariates = covars,
        outcome = outcome)

      task_target_75_nontarget_25 <- make_sl3_Task(
        data = target_75_nontarget_25,
        covariates = covars,
        outcome = outcome)

      # all non-target at 75%

      task_target_25_nontarget_75 <- make_sl3_Task(
        data = target_25_nontarget_75,
        covariates = covars,
        outcome = outcome)

      task_target_75_nontarget_75 <- make_sl3_Task(
        data = target_75_nontarget_75,
        covariates = covars,
        outcome = outcome)

      # all non-target at 50%

      task_target_25_nontarget_50 <- make_sl3_Task(
        data = target_25_nontarget_50,
        covariates = covars,
        outcome = outcome)

      task_target_75_nontarget_50 <- make_sl3_Task(
        data = target_75_nontarget_50,
        covariates = covars,
        outcome = outcome)

      # all non-target at observed

      # task_target_25_nontarget_obs <- make_sl3_Task(
      #   data = target_25_nontarget_obs,
      #   covariates = covars,
      #   outcome = outcome)
      #
      # task_target_75_nontarget_obs <- make_sl3_Task(
      #   data = target_75_nontarget_obs,
      #   covariates = covars,
      #   outcome = outcome)


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

      return(result)
    }

    quantile_importance_list[[iter]] <- quantile_importance
    remaining <- remaining[remaining %notin% quantile_importance[,1]]
    print(remaining)
    iter <- iter + 1
    print(iter)
  }

  quantile_importance <- do.call(rbind, quantile_importance_list)
  quantile_results <- data.table(quantile_importance)
  colnames(quantile_results) <- c("X", "Delta Rest High", "Delta Rest Medium", "Delta Rest Low", "Interaction Metric")
  quantile_results_ordered <- quantile_results[order(-quantile_results$`Interaction Metric`)]


  print("Finished Quantile Interaction Search")

  #####################################################################################
  ######################## SUBGROUP QUANTILE IMPORTANCE ###############################
  #####################################################################################

  remaining <- unique(subcategories)
  remaining <- remaining[remaining!="outcome"]
  subgroup_quantile_importance_list <- list()
  iter <- 1

  while (length(remaining) > 0) {
    subgroup_quantile_importance <- foreach(i = remaining, .combine = 'rbind', .errorhandling = "pass") %dopar% {

      subcat_vars <- variable_list[which(subcategories %in% i)]

      subcat_25_nontarget_obs <-set_quantiles(data = dat, X, target = subcat_vars, target_q = 0.25, nontarget_q = NULL, subcategory_flag = TRUE)
      subcat_75_nontarget_obs <-set_quantiles(data = dat, X, target = subcat_vars, target_q = 0.75, nontarget_q = NULL, subcategory_flag = TRUE)

      task_sub_cat_25_nontarget_obs <- make_sl3_Task(
        data = subcat_25_nontarget_obs,
        covariates = covars,
        outcome = outcome)

      task_sub_cat_75_nontarget_obs <- make_sl3_Task(
        data = subcat_75_nontarget_obs,
        covariates = covars,
        outcome = outcome)

      subcat_75_obs_predictions <- fit$predict_fold(task = task_sub_cat_75_nontarget_obs, fold_number = "validation")
      subcat_25_obs_predictions <- fit$predict_fold(task = task_sub_cat_25_nontarget_obs, fold_number = "validation")


      varimp_metric <- mean(subcat_75_obs_predictions - subcat_25_obs_predictions)
      result <- cbind(i, varimp_metric)
      return(result)
    }

    subgroup_quantile_importance_list[[iter]] <- subgroup_quantile_importance
    remaining <- remaining[remaining %notin% subgroup_quantile_importance[,1]]
    print(remaining)
    iter <- iter + 1
    print(iter)
  }

  subgroup_quantile_importance <- do.call(rbind, subgroup_quantile_importance_list)
  subgroup_quantile_importance <- data.table(subgroup_quantile_importance)
  colnames(subgroup_quantile_importance) <- c("X", "Sub-Category Delta")
  subgroup_quantile_results_ordered <- subgroup_quantile_importance[order(-subgroup_quantile_importance$`Sub-Category Delta`)]

  total <- sum(dat[[outcome]] * dat$Population)

  subgroup_quantile_results_ordered$`Sub-Category Delta` <- as.numeric(subgroup_quantile_results_ordered$`Sub-Category Delta`) * total

  merged_results <- merge(risk_results_ordered, quantile_results_ordered, by= "X", all.x = TRUE, all.y = TRUE)
  merged_results <- subset(merged_results, merged_results$X != "CentroidLon" & merged_results$X != "CentroidLat" & merged_results$X != "Latitude"  & merged_results$X != "Longitude")
  risk_results <- subset(risk_results, risk_results$X != "CentroidLon" & risk_results$X != "CentroidLat" & risk_results$X != "Latitude"  & risk_results$X != "Longitude")

  merged_results$X<- Data_Dictionary$`Nice Label`[match(merged_results$X, Data_Dictionary$`Variable Name`)]

  merged_results[,2:6] <- sapply(merged_results[,2:6], as.numeric)
  merged_results[3:6] <- merged_results[3:6] * total

  risk <- risk * total
  variable_combinations <- combn(subset(risk_results, risk_ratio > quantile(merged_results$risk_ratio, .95))$X, m = m)
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

    permuted_importance <- foreach(i = 1:length(remaining), .combine = 'rbind') %dopar% {
      target_vars <- remaining[i]
      target_vars <- unlist(str_split(target_vars, " & "))

      ## compute the additive risk for this set of variables
      additives <- risk_results %>% filter(X %in% target_vars) %>% select(risk_ratio)
      additive_risk <- sum(unlist(additives) - 1) + 1

      ## calculate the permuted risk for this set of variables
      scrambled_rows <- dat[sample(nrow(dat)), ]
      scrambled_rows_selection <- scrambled_rows %>% dplyr::select(!!target_vars)
      scrambled_col_names <- task$add_columns(scrambled_rows_selection)
      scrambled_col_task <- task$next_in_chain(column_names = scrambled_col_names)
      scrambled_sl_preds <- fit$predict_fold(scrambled_col_task,
                                             fold_number = "validation")

      risk_scrambled <- mean(loss(scrambled_sl_preds, Y))
      varimp_metric <- risk_scrambled/risk

      target_vars_nice <- Data_Dictionary$`Nice Label`[match(target_vars, Data_Dictionary$`Variable Name`)]

      result <- cbind(paste(target_vars_nice, collapse = " & "), varimp_metric, additive_risk, remaining[[i]])
      result
    }

    joint_perm_importance_list[[iter]] <- permuted_importance
    remaining <- remaining[remaining %notin% permuted_importance[,4]]
    print(remaining)
    iter <- iter + 1
    print(iter)
  }

  permuted_importance <- do.call(rbind, joint_perm_importance_list)

  print("Finished Joint Permutation")

  permuted_importance <- as.data.frame(permuted_importance)
  permuted_importance <- permuted_importance[1:3]
  permuted_importance$diff <- round(as.numeric(permuted_importance$varimp_metric) - as.numeric(permuted_importance$additive_risk), 3)
  colnames(permuted_importance)[1] <- "Variable Combo"

  permuted_importance <- permuted_importance[order(permuted_importance$diff,decreasing = TRUE),]

  test <- permuted_importance[1:10, ]#subset(permuted_importance, diff >= quantile(permuted_importance$diff , .95))
  test <- test[!is.na(test$diff),]
  test <- melt(test, id.vars=c("Variable Combo", "diff"))
  test$value <- round(as.numeric(test$value),3)

  test$variable <- factor(test$variable, levels=c("varimp_metric", "additive_risk"), labels=c("Joint Risk", "Additive Risk"))
  colnames(test)[3] <- "Type"


  risk_plot <- merged_results %>%
    arrange(risk_ratio) %>%    # First sort by val. This sort the dataframe but NOT the factor levels
    filter(risk_ratio > thresh)  %>%
    mutate(name=factor(X, levels=X)) %>%   # This trick update the factor levels
    ggplot( aes(x=name, y=risk_ratio)) +
    geom_segment( aes(xend=name, yend=1)) +
    geom_point( size=4, color="orange") +
    coord_flip() +
    theme_bw(base_size = 12) +
    ylab("Model Risk Ratio") +
    xlab("County Feature")

  # quantile plots

  merged_results_plot <- merged_results %>%
    arrange(risk_ratio) %>%    # First sort by val. This sort the dataframe but NOT the factor levels
    filter(risk_ratio > thresh)  %>%
    mutate(name=factor(X, levels=X))

  quantile_results_long <- melt(merged_results_plot[c(3:7)], id.vars="name")

  quantiles_plot <- quantile_results_long %>%
    filter(variable !=  "Interaction Metric") %>%
    ggplot(aes(name, value, col=variable)) +
    geom_point(size=4) +
    coord_flip() +
    theme_bw(base_size = 12)

  interaction_plot <- quantile_results_long %>%
    filter(variable ==  "Interaction Metric") %>%
    ggplot(aes(name, value)) +
    geom_segment( aes(xend=name, yend=0)) +
    geom_point(size=4, color="blue") +
    coord_flip() +
    theme_bw(base_size = 12)

  joint_permutation_plot <- test %>%
    group_by(`Variable Combo`) %>%
    ggplot(aes(x= value, y= reorder(`Variable Combo`,value))) +
    geom_line(aes(group = `Variable Combo`),color="grey") +
    geom_point(aes(color=Type), size=4) +
    labs(y="Combination") +
    ylab("County Features") +
    xlab("Model Risk Ratio") +
    theme_bw(base_size = 12)

  sub_group_risk_plot <- subgroup_importance_ordered %>%
    arrange(subgroup_risk_ratio) %>%    # First sort by val. This sort the dataframe but NOT the factor levels
    mutate(name=factor(X, levels=X)) %>%   # This trick update the factor levels
    ggplot( aes(x=name, y=subgroup_risk_ratio)) +
    geom_segment( aes(xend=name, yend=1)) +
    geom_point( size=4, color="orange") +
    coord_flip() +
    theme_bw(base_size = 12) +
    ylab("Model Risk Ratio") +
    xlab("County Feature")

  sub_group_quantile_plot <- subgroup_quantile_results_ordered %>%
    arrange(`Sub-Category Delta`) %>%    # First sort by val. This sort the dataframe but NOT the factor levels
    mutate(name=factor(X, levels=X)) %>%   # This trick update the factor levels
    ggplot( aes(x=name, y=`Sub-Category Delta`)) +
    geom_segment( aes(xend=name, yend=0)) +
    geom_point( size=4, color="orange") +
    coord_flip() +
    theme_bw(base_size = 12) +
    ylab("Model Risk Ratio") +
    xlab("County Feature")


  # p <- plot_grid(risk_plot, quantiles_plot, interaction_plot, labels = label, vjust = -0.1, nrow = 1)
  ggsave(here(paste("Figures/", "varimp_", label, ".png", sep = "")), risk_plot,  width = 8, height = 6)
  ggsave(here(paste("Figures/", "quantile_imp_", label, ".png", sep = "")), quantiles_plot,  width = 8, height = 6)
  ggsave(here(paste("Figures/", "interaction_imp_", label, ".png", sep = "")), interaction_plot,  width = 8, height = 6)
  ggsave(here(paste("Figures/", "jointimp_", label, ".png", sep = "")), joint_permutation_plot, width = 14, height = 6)
  ggsave(here(paste("Figures/", "subgroup_imp_", label, ".png", sep = "")), sub_group_risk_plot, width = 8, height = 6)
  ggsave(here(paste("Figures/", "subgroup_quantile_imp_", label, ".png", sep = "")), sub_group_quantile_plot, width = 8, height = 6)


  return(list("indiv_results" = merged_results, "joint_results"= test, "model_risk" = risk))
}


fit_sl_varimp <- function(outcome,label) {

  ## load data
  covid_data_processed <- read.csv(PROCESSED_DATA_PATH("cleaned_covid_data_final_Mar_31_22.csv"), check.names = FALSE)
  covid_data_processed <- covid_data_processed[,-1]

  Data_Dictionary <- read_excel(PROCESSED_DATA_PATH("Data_Dictionary.xlsx"))
  Data_Dictionary <- Data_Dictionary[Data_Dictionary$`Variable Name` %in% colnames(covid_data_processed),]

  covid_data_processed$CountyRelativeDay100Cases <- covid_data_processed$CountyRelativeDay100Cases / covid_data_processed$Population
  covid_data_processed$TotalCasesUpToDate <- covid_data_processed$TotalCasesUpToDate / covid_data_processed$Population
  covid_data_processed$CountyRelativeDay100Deaths <- covid_data_processed$CountyRelativeDay100Deaths / covid_data_processed$Population
  covid_data_processed$TotalDeathsUpToDate <- covid_data_processed$TotalDeathsUpToDate / covid_data_processed$Population
  covid_data_processed$Deathsat1year <- covid_data_processed$Deathsat1year / covid_data_processed$Population
  covid_data_processed$Casesat1year <- covid_data_processed$Casesat1year / covid_data_processed$Population

  outcomes <- c("CountyRelativeDay100Cases",
                "TotalCasesUpToDate",
                "CountyRelativeDay100Deaths" ,
                "TotalDeathsUpToDate",
                "Deathsat1year",
                "Casesat1year")


  covars <- colnames(covid_data_processed)[-which(names(covid_data_processed) %in% c(
    outcomes,
    "fips",
    "county_names"
  ))]

  task <- make_sl3_Task(
    data = covid_data_processed,
    covariates = covars,
    outcome = outcome,
    folds = origami::make_folds(covid_data_processed, fold_fun = folds_vfold, V = 10))

  discrete_sl <- source(here("R/utils_create_sl.R"))
  discrete_sl <- discrete_sl$value()
  ## fit the sl3 object
  sl_fit <- discrete_sl$train(task)

  ## get variable importance from the sl3 object
  var_importance <- run_varimp(fit = sl_fit,
                               loss = loss_squared_error,
                               covars = covars,
                               outcome = outcome,
                               data = covid_data_processed,
                               Data_Dictionary = Data_Dictionary,
                               label = label,
                               thresh =  1.01)

  SL_results <- list('fit' = sl_fit, 'var_imp' = var_importance)

  saveRDS(SL_results, here(paste("Models/", outcome, ".RDS", sep = "")))

  return(NULL)
}

Cases1year <- fit_sl_varimp(outcome = "Casesat1year", label = "COVID-19 Cases at 1 Year")

