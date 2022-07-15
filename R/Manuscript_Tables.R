library(here)
source(here("R/utils_sl_varimp.R"))
source(here("R/util.R"))

create_manuscript_tables <- function(path, outcome){

  data <- readRDS(path)
  ate_data <- data %>% filter(Condition == "ATE")
  vte_data <- data %>% filter(Condition == "VTE")
  max_var_ate_data <- data %>% filter(Condition == "Blip_Intxn")

  write.csv(ate_data, here("results", paste("ATE",outcome, ".csv", sep = "")))
  write.csv(vte_data, here("results", paste("VTE",outcome, ".csv", sep = "")))
  write.csv(max_var_ate_data, here("results", paste("INTXN",outcome, ".csv", sep = "")))
}

paths <- c("~/COVIDxrisk/data/TotalCasesUpToDate_ind_var_imp_quantile.RDS",
           "~/COVIDxrisk/data/TotalDeathsUpToDate_ind_var_imp_quantile.RDS",
           "~/COVIDxrisk/data/Deathsat1year_ind_var_imp_quantile.RDS",
           "~/COVIDxrisk/data/Casesat1year_ind_var_imp_quantile.RDS",
           "~/COVIDxrisk/data/CountyRelativeDay100Cases_ind_var_imp_quantile.RDS",
           "~/COVIDxrisk/data/CountyRelativeDay100Deaths_ind_var_imp_quantile.RDS")

outcomes <- c("totalcases",
                 "totaldeaths",
                 "deaths1year",
                 "cases1year",
                 "casesday100",
                 "deathsday100")

mapply(create_manuscript_tables, paths, outcomes)

