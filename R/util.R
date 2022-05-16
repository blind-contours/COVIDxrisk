check_packages = function(p){
  if(!require(p, character.only = TRUE)){
    install.packages(p)
  }
  library(p, character.only = TRUE)
}

packages <- c("flexclust", "furrr", "caret", "rvest", "dplyr",  "origami", "knitr", "ggplot2", "R6", "gam", "gbm", "remotes",
              "tidyverse", "here", "tidyr", "readxl", "pROC", "imputeTS", "data.table", "xgboost", "zeallot",
              "pheatmap", "ggfortify", "gtools", "readxl", "readr", "origami", "ggplot2", "R6", "data.table", "cowplot", "future")

lapply(packages, check_packages)

# install superLearner
if(!require("sl3", character.only = TRUE)) {
  remotes::install_github("tlverse/sl3")
  library("sl3")

}

`%notin%` <- Negate(`%in%`)

RAW_DATA_PATH = function(file_name){
  return(paste(here('data/raw/'), file_name, sep=''))
}

PROCESSED_DATA_PATH = function(file_name){
  return(paste(here('data/processed/'), file_name, sep=''))
}
