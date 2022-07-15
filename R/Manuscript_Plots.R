library(here)
source(here("R/util.R"))
library(usmap)
library(ggplot2)

COVIDxRisk_ATE_plot <- function(path, time_frame){
  print(path)
  raw <- readRDS(path)
  ATE <- raw %>% filter(Condition == "ATE")
  ATE_pos <- ATE %>% filter(Est > 0) %>% arrange(desc(Est)) %>% top_n(20, Est)

  if (grepl("case", tolower(path))) {
    color <- "orange"
    outcome <- "cases"
  }else{
    color <- "red"
    outcome <- "deaths"
  }

  label <- paste("ATE",outcome, time_frame, sep = "_")

  plot <- ATE_pos %>%
    arrange(Est) %>%
    # First sort by val. This sort the dataframe but NOT the factor levels
    mutate(name = factor(Label, levels = Label)) %>%
    # This trick update the factor levels
    ggplot(aes(x = name, y = Est)) +
    geom_point(size = 2, color = color) +
    geom_errorbar(aes(ymin = `Lower_CI`, ymax = `Upper_CI`), width = .3, color = color ) +
    coord_flip() +
    theme_bw(base_size = 12)

  ggsave(here(paste("Figures/", "varimp_", label, ".png", sep = "")), plot, width = 8, height = 6)
}

paths <- c("~/COVIDxrisk/data/TotalCasesUpToDate_ind_var_imp_quantile.RDS",
           "~/COVIDxrisk/data/TotalDeathsUpToDate_ind_var_imp_quantile.RDS",
           "~/COVIDxrisk/data/Deathsat1year_ind_var_imp_quantile.RDS",
           "~/COVIDxrisk/data/Casesat1year_ind_var_imp_quantile.RDS",
           "~/COVIDxrisk/data/CountyRelativeDay100Cases_ind_var_imp_quantile.RDS",
           "~/COVIDxrisk/data/CountyRelativeDay100Deaths_ind_var_imp_quantile.RDS")

time_frames <- c("total",
                 "total",
                 "1year",
                 "1year",
                 "day100",
                 "day100")

mapply(COVIDxRisk_ATE_plot, paths, time_frames)

################################################################################
############################## VTE PLOT ########################################
################################################################################

COVIDxRisk_VTE_plot <- function(path, time_frame){
  print(path)
  raw <- readRDS(path)
  VTE <- raw %>% filter(Condition == "VTE")
  VTE_pos <- VTE %>% filter(Est > 0) %>% arrange(desc(Est)) %>% top_n(20, Est)

  if (grepl("case", tolower(path))) {
    color <- "blue"
    outcome <- "cases"
  }else{
    color <- "darkgreen"
    outcome <- "deaths"
  }

  label <- paste("VTE",outcome, time_frame, sep = "_")

  plot <- VTE_pos %>%
    arrange(Est) %>%
    # First sort by val. This sort the dataframe but NOT the factor levels
    mutate(name = factor(Label, levels = Label)) %>%
    # This trick update the factor levels
    ggplot(aes(x = name, y = Est)) +
    geom_point(size = 2, color = color) +
    geom_errorbar(aes(ymin = `Lower_CI`, ymax = `Upper_CI`), width = .3, color = color ) +
    coord_flip() +
    theme_bw(base_size = 12)

  ggsave(here(paste("Figures/", "varimp_", label, ".png", sep = "")), plot, width = 8, height = 6)
}

paths <- c("~/COVIDxrisk/data/TotalCasesUpToDate_ind_var_imp_quantile.RDS",
           "~/COVIDxrisk/data/TotalDeathsUpToDate_ind_var_imp_quantile.RDS",
           "~/COVIDxrisk/data/Deathsat1year_ind_var_imp_quantile.RDS",
           "~/COVIDxrisk/data/Casesat1year_ind_var_imp_quantile.RDS",
           "~/COVIDxrisk/data/CountyRelativeDay100Cases_ind_var_imp_quantile.RDS",
           "~/COVIDxrisk/data/CountyRelativeDay100Deaths_ind_var_imp_quantile.RDS")

time_frames <- c("total",
                 "total",
                 "1year",
                 "1year",
                 "day100",
                 "day100")

mapply(COVIDxRisk_VTE_plot, paths, time_frames)


################################################################################
################## MAX QUARTILE ATE VARIANCE ###################################
################################################################################

COVIDxRisk_max_var_plot <- function(path, time_frame){
  print(path)
  raw <- readRDS(path)
  Blip_intxn <- raw %>% filter(Condition == "Blip_Intxn")
  Max_ATE_var <- Blip_intxn %>% filter(Est > 0) %>% arrange(desc(Est)) %>% top_n(20, Est)
  Data_Dictionary <- read_excel(PROCESSED_DATA_PATH( "Data_Dictionary.xlsx"))

  Max_ATE_var$EM_label <- Data_Dictionary$`Nice Label`[match(Max_ATE_var$Intxn_Var, Data_Dictionary$`Variable Name`)]
  Max_ATE_var$Var_Intxn <- paste(Max_ATE_var$Label, Max_ATE_var$EM_label, sep = "-in-")

  if (grepl("case", tolower(path))) {
    color <- "turquoise1"
    outcome <- "cases"
  }else{
    color <- "violet"
    outcome <- "deaths"
  }

  label <- paste("Blip_Intxn",outcome, time_frame, sep = "_")

  plot <- Max_ATE_var %>%
    arrange(Est) %>%
    # First sort by val. This sort the dataframe but NOT the factor levels
    mutate(name = factor(Var_Intxn, levels = Var_Intxn)) %>%
    # This trick update the factor levels
    ggplot(aes(x = name, y = Est)) +
    geom_point(size = 2, color = color) +
    geom_errorbar(aes(ymin = `Lower_CI`, ymax = `Upper_CI`), width = .3, color = color ) +
    coord_flip() +
    theme_bw(base_size = 12)

  ggsave(here(paste("Figures/", "varimp_", label, ".png", sep = "")), plot, width = 10, height = 6)
  }

paths <- c("~/COVIDxrisk/data/TotalCasesUpToDate_ind_var_imp_quantile.RDS",
           "~/COVIDxrisk/data/TotalDeathsUpToDate_ind_var_imp_quantile.RDS",
           "~/COVIDxrisk/data/Deathsat1year_ind_var_imp_quantile.RDS",
           "~/COVIDxrisk/data/Casesat1year_ind_var_imp_quantile.RDS",
           "~/COVIDxrisk/data/CountyRelativeDay100Cases_ind_var_imp_quantile.RDS",
           "~/COVIDxrisk/data/CountyRelativeDay100Deaths_ind_var_imp_quantile.RDS")

time_frames <- c("total",
                 "total",
                 "1year",
                 "1year",
                 "day100",
                 "day100")

mapply(COVIDxRisk_max_var_plot, paths, time_frames)


################################################################################
################## MAX QUARTILE ATE VARIANCE ###################################
################################################################################

COVIDxRisk_ATE_qua_spec_plot <- function(path, time_frame, target_var){
  print(path)
  raw <- readRDS(path)
  ATE_in_Quarts <- raw %>% filter(Target_var == target_var) %>% filter(Condition == "Blip_preds")
  Data_Dictionary <- read_excel(PROCESSED_DATA_PATH( "Data_Dictionary.xlsx"))

  ATE_in_Quarts$target_var <- Data_Dictionary$`Nice Label`[match(ATE_in_Quarts$Target_var, Data_Dictionary$`Variable Name`)]
  ATE_in_Quarts$em_var <- Data_Dictionary$`Nice Label`[match(ATE_in_Quarts$EM_var, Data_Dictionary$`Variable Name`)]

  if (grepl("case", tolower(path))) {
    color <- "gold1"
    outcome <- "cases"
  }else{
    color <- "navy"
    outcome <- "deaths"
  }

  label <- paste("ATE_by_Quartile", outcome, time_frame, sep = "_")

  EM_quartile_plot <- ATE_in_Quarts %>%
    ggplot(aes(x = EM_quantile, y = `50%`)) +
    geom_point(size = 2, color = color) +
    geom_errorbar(aes(ymin = `2.5%`, ymax = `97.5%`), width = .3, color = color ) +
    theme_bw(base_size = 12) +
    labs(x = paste("Quartiles of", ATE_in_Quarts$em_var, sep = " "), y = paste(ATE_in_Quarts$target_var, " ATE"))

  ggsave(here(paste("Figures/", "varimp_", label, ".png", sep = "")), EM_quartile_plot, width = 10, height = 6)
}


target_vars <- c("pm10", "Life.expectancy.raw.value", "SIMAZINE", "SIMAZINE", "whitewhite", "FirstCaseDay")


paths <- c("~/COVIDxrisk/data/TotalCasesUpToDate_intxn_imp_quantile.RDS",
           "~/COVIDxrisk/data/TotalDeathsUpToDate_intxn_imp_quantile.RDS",
           "~/COVIDxrisk/data/Casesat1year_intxn_imp_quantile.RDS",
           "~/COVIDxrisk/data/Deathsat1year_intxn_imp_quantile.RDS",
           "~/COVIDxrisk/data/CountyRelativeDay100Cases_intxn_imp_quantile.RDS",
           "~/COVIDxrisk/data/CountyRelativeDay100Deaths_intxn_imp_quantile.RDS")

time_frames <- c("total",
                 "total",
                 "1year",
                 "1year",
                 "day100",
                 "day100")

mapply(COVIDxRisk_ATE_qua_spec_plot, paths, time_frames, target_vars)



US_county_ATE_plot <- function(outcome,
                               path_model,
                               target_var,
                               states){

  path_data <- "cleaned_covid_data_final.csv"

  data <- read.csv(PROCESSED_DATA_PATH(path_data),
                                     check.names = FALSE)
  data <- data[, -1]

  sl <- readRDS(path_model)

  all_outcomes <- c(
    "CountyRelativeDay100Cases",
    "TotalCasesUpToDate",
    "CountyRelativeDay100Deaths",
    "TotalDeathsUpToDate",
    "Deathsat1year",
    "Casesat1year"
  )

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

  quantiles <- quantile(data[[target_var]])

  data_Q1 <- data_Q4 <- data

  data_Q1[[target_var]] <- quantiles[2]
  data_Q4[[target_var]] <- quantiles[4]

  Q1_task <- make_sl3_Task(
    data = data_Q1,
    covariates = covars,
    outcome = outcome
  )

  Q4_task <- make_sl3_Task(
    data = data_Q4,
    covariates = covars,
    outcome = outcome
  )

  total <- sum(data[[outcome]])

  Q1_predictions <- sl$predict_fold(task = Q1_task, fold_number = "full") * total
  Q4_predictions <- sl$predict_fold(task = Q4_task, fold_number = "full") * total

  blip <- Q4_predictions - Q1_predictions

  df <- data.frame(
    fips = data$fips,
    values = blip
  )

  us_counties_select_plot <- plot_usmap(data = df, include = c(states)) + scale_fill_continuous(type = "viridis")
  us_counties_all_plot <- plot_usmap(data = df) + scale_fill_continuous(type = "viridis")

  ggsave(here(paste("Figures/", "ATE_region_", outcome,"_", target_var, ".png", sep = "")), us_counties_select_plot, width = 10, height = 6)
  ggsave(here(paste("Figures/", "ATE_all_",  outcome, "_", target_var, ".png", sep = "")), us_counties_all_plot, width = 10, height = 6)
  }

# Total Cases -------------------------
US_county_ATE_plot(outcome = "TotalCasesUpToDate",
                   target_var = "DIMETHENAMID-P",
                   path_model <- "~/COVIDxrisk/Models/TotalCasesUpToDate.RDS",
                   states = c("DC", "DE", "FL", "GA", "MD", "NC", "SC", "VA",
                              "WV", "AL", "KY", "MS", "TN", "AR", "LA", "OK",
                              "TX", "PA"))

US_county_ATE_plot(outcome = "TotalCasesUpToDate",
                   target_var = "pm10",
                   path_model <- "~/COVIDxrisk/Models/Models/TotalCasesUpToDate.RDS",
                   states = c("DC", "DE", "FL", "GA", "MD", "NC", "SC", "VA",
                              "WV", "AL", "KY", "MS", "TN", "AR", "LA", "OK",
                              "TX", "PA"))

US_county_ATE_plot(outcome = "TotalCasesUpToDate",
                   target_var = "2,4-D",
                   path_model <- "~/COVIDxrisk/Models/Models/TotalCasesUpToDate.RDS",
                   states = c("DC", "DE", "FL", "GA", "MD", "NC", "SC", "VA",
                              "WV", "AL", "KY", "MS", "TN", "AR", "LA", "OK",
                              "TX", "PA"))


# Total Deaths -------------------------
US_county_ATE_plot(outcome = "TotalDeathsUpToDate",
                   target_var = "Life.expectancy.raw.value",
                   path_model <- "~/COVIDxrisk/Models/Models/TotalDeathsUpToDate.RDS",
                   states = c("DC", "DE", "FL", "GA", "MD", "NC", "SC", "VA",
                              "WV", "AL", "KY", "MS", "TN", "AR", "LA", "OK",
                              "TX", "PA"))

US_county_ATE_plot(outcome = "TotalDeathsUpToDate",
                   target_var = "all_heartdisease_deathrate",
                   path_model <- "~/COVIDxrisk/Models/Models/TotalDeathsUpToDate.RDS",
                   states = c("DC", "DE", "FL", "GA", "MD", "NC", "SC", "VA",
                              "WV", "AL", "KY", "MS", "TN", "AR", "LA", "OK",
                              "TX", "PA"))

US_county_ATE_plot(outcome = "TotalDeathsUpToDate",
                   target_var = "Median.household.income.raw.value",
                   path_model <- "~/COVIDxrisk/Models/Models/TotalDeathsUpToDate.RDS",
                   states = c("DC", "DE", "FL", "GA", "MD", "NC", "SC", "VA",
                              "WV", "AL", "KY", "MS", "TN", "AR", "LA", "OK",
                              "TX", "PA"))

# Day 100 Cases -------------------------

US_county_ATE_plot(outcome = "CountyRelativeDay100Cases",
                   target_var = "rep_ratio",
                   path_model <- "~/COVIDxrisk/Models/Models/TotalCasesUpToDate.RDS",
                   states = c("DC", "DE", "FL", "GA", "MD", "NC", "SC", "VA",
                              "WV", "AL", "KY", "MS", "TN", "AR", "LA", "OK",
                              "TX", "PA"))


US_county_ATE_plot(outcome = "CountyRelativeDay100Cases",
                   target_var = "whitewhite",
                   path_model <- "~/COVIDxrisk/Models/Models/CountyRelativeDay100Cases.RDS",
                   states = c("DC", "DE", "FL", "GA", "MD", "NC", "SC", "VA",
                              "WV", "AL", "KY", "MS", "TN", "AR", "LA", "OK",
                              "TX", "PA"))


# Day 100 Deaths -------------------------

US_county_ATE_plot(outcome = "CountyRelativeDay100Deaths",
                   target_var = "BWI",
                   path_model <- "~/COVIDxrisk/Models/Models/TotalDeathsUpToDate.RDS",
                   states = c("DC", "DE", "FL", "GA", "MD", "NC", "SC", "VA",
                              "WV", "AL", "KY", "MS", "TN", "AR", "LA", "OK",
                              "TX", "PA"))


df <- data.frame(
  fips = data$fips,
  values = data$avg_arsenic_levels
)

us_counties_select_plot <- plot_usmap(data = df, include = c( "FL", "GA", "NC", "SC", "VA")) + scale_fill_continuous(type = "viridis")

