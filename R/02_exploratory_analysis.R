library(tidyverse)
library(imputeTS)
library(pheatmap)
library(ggfortify)
library(gtools)
library(readxl)
library(here)

## load data after preprocessing
#covid_data_processed <- read_excel(here("Analysis/update_data/data/processed/cleaned_covid_data_final.xlsx"),sheet = 1)

covid_data_processed <- read_csv(PROCESSED_DATA_PATH("cleaned_covid_data_final_Mar_31_22.csv"))

# get data dictionary
Data_Dictionary <- read_excel(PROCESSED_DATA_PATH("Data_Dictionary.xlsx"))

vars_2_keep <- Data_Dictionary %>%
  filter(Keep == "Yes")

## create quantiles of the outcome data to use as labels in the clustering analyses
outcomes <- c("CountyRelativeDay100Cases",
              "TotalCasesUpToDate",
              "CountyRelativeDay100Deaths",
              "TotalDeathsUpToDate",
              "FirstCaseDay",
              "Deathsat1year",
              "Casesat1year")

features_data <- covid_data_processed %>%
  select(-c(outcomes, '...1', 'fips'))

outcomes_data <- covid_data_processed %>%
  select(outcomes)

## map quantcut from gtools across each column of outcomes_data
outcome_quantiles <- purrr::map(outcomes_data, quantcut, q=4)

## relabel the quantiles
for (i in seq_along(1:length(outcome_quantiles))) {
  levels(outcome_quantiles[[i]]) <-  c("Q1", "Q2", "Q3", "Q4")
}

## Heatmap

## make dataframe of each outcome, row names FIPS to do multiple annotations in the heatmap

## day 25 cases
cases_25_df = data.frame("Quantile_Day_25_Cases" = outcome_quantiles$CountyRelativeDay100Cases)
## total cases
cases_total_df = data.frame("Quantile_Total_Cases" = outcome_quantiles$TotalCasesUpToDate)
## total cases
deaths_100_df = data.frame("Quantile_Day_100_Deaths" = outcome_quantiles$CountyRelativeDay100Deaths)
## TotalDeathsUpToDate
deaths_total_df = data.frame("Quantile_Total_Deaths" = outcome_quantiles$TotalDeathsUpToDate)
## TotalDeathsUpToDate
first_day_df = data.frame("Quantile_Day_First_Case" = outcome_quantiles$FirstCaseDay)

deathsat1year = data.frame("Quantile_One_Year_Deaths" = outcome_quantiles$Deathsat1year)

Casesat1year = data.frame("Quantile_One_Year_Cases" = outcome_quantiles$Casesat1year)


oucome_annotation <- cbind(cases_25_df,
                           cases_total_df,
                           deaths_100_df,
                           deaths_total_df,
                           first_day_df,
                           deathsat1year,
                           Casesat1year)


covid_numerical <- as.matrix(features_data)
rownames(covid_numerical) <- covid_data_processed$fips

#scale the data
covid_num_scale <- scale(covid_numerical)

write.csv(x = colnames(covid_num_scale), file = PROCESSED_DATA_PATH("column_dendrogram_labeling.csv"))

col_labels <- read.csv(PROCESSED_DATA_PATH("column_dendrogram_labeled.csv"))
col_labels <- col_labels$X.1

col_labels <- as.data.frame(col_labels)
row.names(col_labels) <- colnames(covid_numerical)
col_labels$col_labels <- as.factor(col_labels$col_labels)

my.breaks <- c(seq(-5, 0, by=0.1), seq(0.1, 5, by=0.1))
my.colors <- c(colorRampPalette(colors = c("blue", "white"))(length(my.breaks)/2), colorRampPalette(colors = c("white", "orange", "red", "purple"))(length(my.breaks)/2))

## set rownames for row dendrogram and column dendrogram
rownames(oucome_annotation) <-  rownames(covid_num_scale)

#variable_labels <- as.data.frame(t(variable_labels))
#rownames(variable_labels) <-  colnames(covid_num_scale)
annoCol<-list(Quantile_Day_25_Cases=c(Q1="blue", Q2="red", Q3="orange", Q4="grey"),
              Quantile_Total_Cases=c(Q1="blue", Q2="red", Q3="orange", Q4="grey"),
              Quantile_Day_100_Deaths=c(Q1="blue", Q2="red", Q3="orange", Q4="grey"),
              Quantile_Total_Deaths=c(Q1="blue", Q2="red", Q3="orange", Q4="grey"),
              Quantile_Day_First_Case=c(Q1="blue", Q2="red", Q3="orange", Q4="grey"),
              Quantile_One_Year_Deaths=c(Q1="blue", Q2="red", Q3="orange", Q4="grey"),
              Quantile_One_Year_Cases=c(Q1="blue", Q2="red", Q3="orange", Q4="grey"),
              col_labels = c(`Air Pollution` = "brown",
                             `Climate Change` = "red",
                             `Demography/Geography` = "orange",
                             `Education/Employment` = "cyan",
                              Health = "darkorchid1",
                             Incarceration = "dimgray",
                             Pesticides = "gold",
                             Racism = "yellow",
                             Segregation = "green1",
                             Transit = "lightseagreen",
                             `Vulnerable Populations` = "midnightblue",
                             `Water Contaminants` = "aliceblue"))

## create the final heatmap
covid_factors_heatmap <- pheatmap(covid_num_scale,main = "COVID-19 Heatmap",
         annotation_names_row = FALSE,
         annotation_names_col = FALSE,
         show_colnames = FALSE,
         show_rownames = FALSE,
         scale = "row",
         color = my.colors,
         breaks = my.breaks,
         annotation_row = oucome_annotation,
         annotation_col = col_labels,
         cutree_rows = 4,
         cutree_cols = 5,
         annotation_colors = annoCol)

## get the clusters of variables
clusters_col <- cutree(covid_factors_heatmap$tree_col, k = 5)
clusters_row <- cutree(covid_factors_heatmap$tree_row, k = 4)

clusters_col <- as.data.frame(clusters_col)
clusters_row <- as.data.frame(clusters_row)

colnames_dendro_reordered <- colnames(covid_num_scale)[covid_factors_heatmap$tree_col$order]
rownames_dendro_reordered <- rownames(covid_num_scale)[covid_factors_heatmap$tree_row$order]

high_values_cluster_group1 <- colnames_dendro_reordered[28:32]
high_values_cluster_group2 <- colnames_dendro_reordered[46:53]
low_values_cluster_group1 <- colnames_dendro_reordered[54:64]


raw_data_w_states <- read.csv(here("Analysis/update_data/data/processed/CountiesMergedData_July_15.csv"))

states_cluster_1 <- table(raw_data_w_states$State[match(rownames(subset(clusters_row, clusters_row == 1)), raw_data_w_states$FIPS)])

## from the reordered columns from the dendrogram we now are indexing the variable names for where we see features related to outcome






## JUNK


table(clusters)

hot_spot_1 <- colnames_dendro_reordered[1:24] ## 2
hot_spot_2 <- colnames_dendro_reordered[25:51] ## 4
hot_spot_3 <- colnames_dendro_reordered[53:83] ## 3
hot_spot_4 <- colnames_dendro_reordered[84:93] ## 1
hot_spot_5 <- colnames_dendro_reordered[93:110] ## 5

hot_spot_high_all <- colnames_dendro_reordered[72:101] ## 5

## ugly change later! don't write code like this, you know better
clusters[rownames(clusters) %in% hot_spot_1,] ## so it appears that the cluster on the far right is cluster 2 from cut-tree
clusters[rownames(clusters) %in% hot_spot_2,] ## so it appears that the cluster on the second left is cluster 4 from cut-tree
clusters[rownames(clusters) %in% hot_spot_3,] ## so it appears that the cluster on the second left is cluster 4 from cut-tree
clusters[rownames(clusters) %in% hot_spot_4,] ## so it appears that the cluster on the second left is cluster 4 from cut-tree
clusters[rownames(clusters) %in% hot_spot_5,] ## so it appears that the cluster on the second left is cluster 4 from cut-tree


clusters[rownames(clusters) %in% hot_spot_high_all,] ## so it appears that the cluster on the second left is cluster 4 from cut-tree


subset(clusters, clusters == 4)


