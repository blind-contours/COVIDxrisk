library(here)
source(here("R/utils_sl_varimp.R"))
source(here("R/util.R"))
cpus <- 5
plan(multisession, workers = cpus)

covid_data_processed <- read_csv(PROCESSED_DATA_PATH("cleaned_covid_data_final_Mar_31_22.csv"))

# get data dictionary
Data_Dictionary <- read_excel(PROCESSED_DATA_PATH("Data_Dictionary.xlsx"))

vars_2_keep <- Data_Dictionary %>%
  filter(Keep == "Yes")

subcategories <- Data_Dictionary[Data_Dictionary$`Variable Name`
                                 %in% colnames(covid_data_processed), ]$Label
subcategories <- subcategories[subcategories != "outcome"]
subcategories <- subcategories[subcategories != "identifier"]

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

cases_100_df = data.frame("Quantile_Day_100_Cases"
                          = outcome_quantiles$CountyRelativeDay100Cases)
cases_total_df = data.frame("Quantile_Total_Cases" =
                              outcome_quantiles$TotalCasesUpToDate)
deaths_100_df = data.frame("Quantile_Day_100_Deaths" =
                             outcome_quantiles$CountyRelativeDay100Deaths)
deaths_total_df = data.frame("Quantile_Total_Deaths" =
                               outcome_quantiles$TotalDeathsUpToDate)
first_day_df = data.frame("Quantile_Day_First_Case" =
                            outcome_quantiles$FirstCaseDay)
deathsat1year = data.frame("Quantile_One_Year_Deaths" =
                             outcome_quantiles$Deathsat1year)
Casesat1year = data.frame("Quantile_One_Year_Cases" =
                            outcome_quantiles$Casesat1year)

oucome_annotation <- cbind(cases_100_df,
                           cases_total_df,
                           deaths_100_df,
                           deaths_total_df,
                           first_day_df,
                           deathsat1year,
                           Casesat1year)

covid_numerical <- as.matrix(features_data)
rownames(covid_numerical) <- covid_data_processed$fips

covid_num_scale <- scale(covid_numerical)

col_labels <- as.data.frame(subcategories)
row.names(col_labels) <- colnames(covid_numerical)
colnames(col_labels) <- "col_labels"
col_labels$col_labels <- as.factor(col_labels$col_labels)

my.breaks <- c(seq(-5, 0, by=0.1), seq(0.1, 5, by=0.1))
my.colors <- c(colorRampPalette(colors =
                                  c("blue", "white"))(length(my.breaks)/2),
               colorRampPalette(colors =
                                  c("white", "orange", "red", "purple"))
               (length(my.breaks)/2))

rownames(oucome_annotation) <-  rownames(covid_num_scale)

annoCol<-list(Quantile_Day_100_Cases=c(Q1="blue", Q2="red",
                                       Q3="orange", Q4="grey"),
              Quantile_Total_Cases=c(Q1="blue", Q2="red",
                                     Q3="orange", Q4="grey"),
              Quantile_Day_100_Deaths=c(Q1="blue", Q2="red",
                                        Q3="orange", Q4="grey"),
              Quantile_Total_Deaths=c(Q1="blue", Q2="red",
                                      Q3="orange", Q4="grey"),
              Quantile_Day_First_Case=c(Q1="blue", Q2="red",
                                        Q3="orange", Q4="grey"),
              Quantile_One_Year_Deaths=c(Q1="blue", Q2="red",
                                         Q3="orange", Q4="grey"),
              Quantile_One_Year_Cases=c(Q1="blue", Q2="red",
                                        Q3="orange", Q4="grey"),
              col_labels = c(`Air Pollution` = "brown",
                             `Climate Change` = "red",
                             `Demography` = "orange",
                             `Geography` = "bisque",
                             `Housing` = "blue",
                             `Education/Employment` = "cyan",
                             Health = "darkorchid1",
                             `Incarceration Metrics` = "dimgray",
                             Occupation = "pink",
                             Pesticides = "gold",
                             `Physical Environment` = "yellow",
                             Segregation = "green1",
                             Transit = "lightseagreen",
                             `Vulnerable Populations` = "midnightblue",
                             `Water Toxicants` = "aliceblue"))

covid_factors_heatmap <- pheatmap(covid_num_scale,
                                  main = "COVID-19 Heatmap",
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
                                  cutree_cols = 20,
                                  annotation_colors = annoCol,
                                  cluster_rows = TRUE,
                                  cluster_cols = TRUE)

## get the clusters of variables
clusters_col <- cutree(covid_factors_heatmap$tree_col, k = 20)
clusters_row <- cutree(covid_factors_heatmap$tree_row, k = 4)

clusters_col <- as.data.frame(clusters_col)
clusters_row <- as.data.frame(clusters_row)

cluster_id_ordered_to_hm <- clusters_col$clusters_col[
  covid_factors_heatmap$tree_col$order]

county_cluster_ordered_to_hm <- clusters_row$clusters_row[
  covid_factors_heatmap$tree_row$order]

county_i_ordered_to_hm <- rownames(clusters_row)[
  covid_factors_heatmap$tree_row$order]

variable_ordered_to_hm <-rownames(clusters_col)[
  covid_factors_heatmap$tree_col$order]

rownames_dendro_reordered <- rownames(covid_num_scale)[
  covid_factors_heatmap$tree_row$order]


vars_clust_hm_ordered <- cbind.data.frame(variable_ordered_to_hm, cluster_id_ordered_to_hm)
counties_clust_hm_ordered <- cbind.data.frame(county_i_ordered_to_hm, county_cluster_ordered_to_hm)


subset(vars_clust_hm_ordered, cluster_id_ordered_to_hm == 1)
high_outcome_counties <- subset(counties_clust_hm_ordered, county_cluster_ordered_to_hm == 2)


fips_state_xwalk <- read_excel(here("data/processed/fips_state_crosswalk.xlsx"))

states_cluster_2 <- fips_state_xwalk$State[match(high_outcome_counties$county_i_ordered_to_hm, fips_state_xwalk$FIPS)]

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




