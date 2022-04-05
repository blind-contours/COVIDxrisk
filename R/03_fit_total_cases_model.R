library(here)
source(here("R/util.R"))
plan(multisession)
cpus <- 20

set_quantiles <- function(data, X, target, target_q, nontarget_q){

  for (i in X) {
    if (i == target) {
      data[[i]] <- quantile(data[[i]], target_q)
    }
    else{
      data[[i]] <- quantile(data[[i]], nontarget_q)
    }
  }
  return(data)

}

run_varimp <- function(fit,
                       loss,
                       covars,
                       outcome,
                       data = covid_data_processed,
                       data_dictionary = Data_Dictionary,
                       label = label,
                       thresh = 1.02,
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

  risk_importance <- foreach(i = X, .combine = 'c', .errorhandling = "pass") %dopar% {
    scrambled_col <- data.table(sample(unlist(dat[, i, with = FALSE]), nrow(dat)))
    names(scrambled_col) <- i
    scrambled_col_names <- task$add_columns(scrambled_col)
    scrambled_col_task <- task$next_in_chain(column_names = scrambled_col_names)
    scrambled_sl_preds <- fit$predict_fold(scrambled_col_task, fold_number = "validation")
    #   i_removed_learner <- fit$reparameterize(list(covariates = setdiff(X, i)))
    #   i_removed_fit <- i_removed_learner$train(task)
    #   i_removed_pred <- i_removed_fit$predict_fold(task, fold_number = "validation")
    no_i_risk <- mean(loss(scrambled_sl_preds, Y))
    varimp_metric <- no_i_risk/risk
    #
    return(varimp_metric)
  }

  quantile_importance <- foreach(i = X, .combine = 'c') %dopar% {

    dat <- fit$training_task$data

    target_25_nontarget_25 <-set_quantiles(data = dat, X, target = i, target_q = 0.25, nontarget_q = 0.25)
    target_75_nontarget_25 <-set_quantiles(data = dat, X, target = i, target_q = 0.75, nontarget_q = 0.25)
    target_25_nontarget_75 <-set_quantiles(data = dat, X, target = i, target_q = 0.25, nontarget_q = 0.75)
    target_75_nontarget_75 <-set_quantiles(data = dat, X, target = i, target_q = 0.75, nontarget_q = 0.75)

    task_target_25_nontarget_25 <- make_sl3_Task(
      data = target_25_nontarget_25,
      covariates = covars,
      outcome = outcome)

    task_target_75_nontarget_25 <- make_sl3_Task(
      data = target_75_nontarget_25,
      covariates = covars,
      outcome = outcome)

    task_target_25_nontarget_75 <- make_sl3_Task(
      data = target_25_nontarget_75,
      covariates = covars,
      outcome = outcome)

    task_target_75_nontarget_75 <- make_sl3_Task(
      data = target_75_nontarget_75,
      covariates = covars,
      outcome = outcome)


    y_25_25_predictions <- fit$predict_fold(task = task_target_25_nontarget_25, fold_number = "validation")
    y_75_25_predictions <- fit$predict_fold(task = task_target_75_nontarget_25, fold_number = "validation")

    y_25_75_predictions <- fit$predict_fold(task = task_target_25_nontarget_75, fold_number = "validation")
    y_75_75_predictions <- fit$predict_fold(task = task_target_75_nontarget_75, fold_number = "validation")

    delta_rest_high <- mean(y_75_75_predictions - y_25_75_predictions)
    delta_rest_low <- mean(y_75_25_predictions - y_25_25_predictions)

    varimp_metric <- delta_rest_high - delta_rest_low
    return(varimp_metric)

  }

  names(risk_importance) <- X
  names(quantile_importance) <- X

  risk_results <- data.table(X = names(risk_importance), risk_ratio = unlist(risk_importance))
  risk_results_ordered <- risk_results[order(-risk_results$risk_ratio)]

  quantile_results <- data.table(X = names(quantile_importance), risk_difference = unlist(quantile_importance))
  quantile_results_ordered <- quantile_results[order(-quantile_results$risk_difference)]

  merged_results <- merge(risk_results_ordered, quantile_results_ordered, by= "X")
  merged_results <- subset(merged_results, merged_results$X != "CentroidLon" & merged_results$X != "CentroidLat" & merged_results$X != "Latitude"  & merged_results$X != "Longitude")
  risk_results <- subset(risk_results, risk_results$X != "CentroidLon" & risk_results$X != "CentroidLat" & risk_results$X != "Latitude"  & risk_results$X != "Longitude")

  merged_results$X<- data_dictionary$`Nice Label`[match(merged_results$X, data_dictionary$`Variable Name`)]

  variable_combinations <- combn(subset(risk_results, risk_ratio > quantile(merged_results$risk_ratio, .97))$X, m = m)
  ### Create list with all intxn_size interactions for the intxn_list variable set of interest:
  variable_combinations <- as.data.frame(variable_combinations)
  ### Run the additive vs. joint error calculation for each set of possible interactions of selected size:

  permuted_importance <- foreach(i = 1:dim(variable_combinations)[2], .combine = 'rbind') %dopar% {
    target_vars <- variable_combinations[,i]
    ## compute the additive risk for this set of variables
    additives <- risk_importance[target_vars]
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

    target_vars <- data_dictionary$`Nice Label`[match(target_vars, data_dictionary$`Variable Name`)]

    result <- cbind(paste(target_vars, collapse = " & "), varimp_metric, additive_risk)
    result
  }

  permuted_importance <- as.data.frame(permuted_importance)
  permuted_importance$diff <- round(as.numeric(permuted_importance$varimp_metric) - as.numeric(permuted_importance$additive_risk), 3)
  colnames(permuted_importance)[1] <- "Variable Combo"

  test <- subset(permuted_importance, diff >= quantile(permuted_importance$diff , .95))
  test <- melt(test, id.vars=c("Variable Combo", "diff"))
  test$value <- round(as.numeric(test$value),3)

  test$variable <- factor(test$variable, levels=c("varimp_metric", "additive_risk"), labels=c("Joint Risk", "Additive Risk"))
  colnames(test)[3] <- "Type"

  risk_plot <- merged_results %>%
    arrange(risk_ratio) %>%    # First sort by val. This sort the dataframe but NOT the factor levels
    filter(risk_ratio > 1.01)  %>%
    mutate(name=factor(X, levels=X)) %>%   # This trick update the factor levels
    ggplot( aes(x=name, y=risk_ratio)) +
    geom_segment( aes(xend=name, yend=1)) +
    geom_point( size=4, color="orange") +
    coord_flip() +
    theme_bw() +
    ylab("Model Risk Ratio") +
    xlab("County Feature")

  total <- sum(dat[[outcome]] * dat$Population)
  merged_results$risk_difference <- merged_results$risk_difference * total

  quantile_plot <- merged_results %>%
    arrange(risk_ratio) %>%    # First sort by val. This sort the dataframe but NOT the factor levels
    filter(risk_ratio > 1.01)  %>%
    mutate(name=factor(X, levels=X)) %>%   # This trick update the factor levels
    ggplot( aes(x=name, y=risk_difference)) +
    geom_segment( aes(xend=name, yend=0)) +
    geom_point( size=4, color="blue") +
    coord_flip() +
    theme_bw() +
    ylab("(75 - 25 | 75 ) − (75 - 25 | 25 )") +
    xlab("")


  joint_permutation_plot <- test %>%
    group_by(`Variable Combo`) %>%
    ggplot(aes(x= value, y= reorder(`Variable Combo`,value))) +
    geom_line(aes(group = `Variable Combo`),color="grey") +
    geom_point(aes(color=Type), size=4) +
    labs(y="Combination") +
    ylab("County Features") +
    xlab("Model Risk Ratio") +
    theme_bw()


  p <- plot_grid(risk_plot, quantile_plot, labels = label, vjust = -0.1)
  ggsave(here(paste("Figures/", "varimp_", label, ".png", sep = "")), p,  width = 15, height = 12)
  ggsave(here(paste("Figures/", "jointimp_", label, ".png", sep = "")), joint_permutation_plot, width = 15, height = 6)

  return(list("indiv_results" = merged_results, "joint_results"= test))
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

  discrete_sl <- source("R/utils_create_sl.R")
  discrete_sl <- discrete_sl$value()
  ## fit the sl3 object
  sl_fit <- discrete_sl$train(task)

  ## get variable importance from the sl3 object
  var_importance <- run_varimp(fit = sl_fit,
                               loss = loss_squared_error,
                               covars = covars,
                               outcome = outcome,
                               data = covid_data_processed,
                               data_dictionary = Data_Dictionary,
                               label = label,
                               thresh = 1.02)

  SL_results <- list('fit' = sl_fit, 'var_imp' = var_importance)

  saveRDS(SL_results, here(paste("Models/processed/", outcome, ".RDS", sep = "")))

  return(NULL)
}


Casestodate <- fit_sl_varimp(outcome = "TotalCasesUpToDate", label = "Total COVID-19 Cases To-Date")

