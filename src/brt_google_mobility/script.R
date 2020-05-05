# Loading Required Libraries
library(tidyverse); library(gbm); library(dismo); library(conflicted); library(gtools); library(lubridate)
conflict_prefer("select", "dplyr")
conflict_prefer("filter", "dplyr")
conflict_prefer("area", "patchwork")

# Loading Google Mobility Data
mob <- read.csv("data/Google_Mobility_Data.csv", stringsAsFactors = FALSE) %>%
  filter(sub_region_1 == "") %>%
  mutate(overall = 1/4 * retail_and_recreation_percent_change_from_baseline +
           1/4 * grocery_and_pharmacy_percent_change_from_baseline + 
           1/4 * transit_stations_percent_change_from_baseline + 
           1/4 * workplaces_percent_change_from_baseline) %>%
  mutate(date = as.Date(date, format="%d/%m/%Y")) %>%
  select(country_region, date, overall)

# Loading World Bank Metadata
wb_metadata <- read.csv("data/World_Bank_Country_Metadata.csv") %>%
  rename(ISO = ï..country_code) %>%
  select(ISO, income_group, region) %>%
  filter(region != "")

# Loading in ACAPs Data
ACAPs_measure <- read.csv("data/ACAPS_Intervention_Data.csv") %>%
  rename(country = COUNTRY, measure = MEASURE, type = LOG_TYPE, date = DATE_IMPLEMENTED) %>%
  select(country, ISO, measure, type, date) %>%
  filter(country != "", measure != "", type != "") %>%
  mutate(measure = as.numeric(factor(measure)), type = as.numeric(factor(type))) %>%
  mutate(combined = paste0(type, "_", measure)) %>%
  mutate(date = as.Date(date, format="%d/%m/%Y")) %>%
  select(country, ISO, date, combined) %>%
  pivot_wider(names_from = combined, values_from = combined) %>%
  filter(!is.na(date))

# Renaming Columns of ACAPs Data and Inputting In Number of Each Measure Put In On Each Date
measures <- colnames(ACAPs_measure)[-(1:3)]
for (i in 1:length(measures)) {
  index <- measures[i]
  ACAPs_measure[, index] <- unlist(lapply(seq_along(ACAPs_measure$country), function(x) {
    length <- length(ACAPs_measure[, index][[1]][[x]])
  }))
  print(i)
}

# Tracking Cumulative Number of Each Type of Measure Implemented and Creating New Rows for Dates 
# Not Present in Each Country
new_ACAPs_measure <- ACAPs_measure %>%
  group_by(country) %>%
  arrange(country, date) %>%
  mutate_at(vars(-country, -ISO, -date), funs(cumsum(.))) %>%    
  complete(date = seq.Date(min(ACAPs_measure$date), max(ACAPs_measure$date), by = "days")) %>%
  mutate_at(vars(-country, -ISO, -date), funs(replace(., row_number() == 1, 0))) 

# Filling Those Newly Created Rows With the Value In the Row Above (As We're Tracking Cumulative)
column_names <- colnames(new_ACAPs_measure)[-c(1:3)]
new_ACAPs_cat <- fill(new_ACAPs_measure, column_names, .direction = c("down"))
new_ACAPs_cat <- fill(new_ACAPs_cat, "ISO", .direction = "updown")

# Joining This Data to the World Bank and Google Mobility Data
overall <- new_ACAPs_cat %>%
  left_join(mob, by = c("country" = "country_region", "date" = "date")) %>% 
  filter(!is.na(overall)) %>%
  left_join(wb_metadata, by = "ISO")

overall_test <- overall %>%
  ungroup(country) %>%
  select(overall, everything(), -country, -ISO, -date)

# Running the BRT 
tree_complexity <- 8
bag_fraction <- 0.5
max_trees <- 3000
learning_rate <- 0.05
x <- as.data.frame(overall_test)
brt <- gbm.step(data = x, 
                gbm.x = 2:66,
                gbm.y = 1,
                family = "gaussian", 
                tree.complexity = tree_complexity, 
                learning.rate = learning_rate, 
                bag.fraction = bag_fraction, 
                max.trees = max_trees, n.folds = 5)

# predicted <- predict.gbm(brt, x[, c(2:66)], n.trees = brt$gbm.call$best.trees, type = "response")
# plot(overall$overall, predicted, ylim = c(-100, 20), xlim = c(-100, 20), pch = 20, cex = 2, ylab = "")

# Get Countries Present in ACAPs But Missing From Google, Infer Google Mobility from ACAPs
google_country_codes <- unique(overall$ISO)
ACAPs_country_codes <- as.character(unique(ACAPs_measure$ISO))
non_google_country_codes <- ACAPs_country_codes[!(ACAPs_country_codes %in% google_country_codes)]
ACAPs_subset <- new_ACAPs_cat[new_ACAPs_cat$ISO %in% non_google_country_codes, ] %>%
  left_join(wb_metadata, by = "ISO")
predicted <- predict.gbm(brt, ACAPs_subset[, c(4:68)], n.trees = brt$gbm.call$best.trees, type = "response")

# Plot Predicted Countries' Google Mobility
# par(mfrow = c(5, 5), mar = c(5, 3, 1, 1))
# ACAPs_subset$pred_mob <- predicted
# predicted_countries <- as.character(unique(ACAPs_subset$country))
# for (i in 1:length(predicted_countries)) {
#   country <- predicted_countries[i]
#   country_index <- which(ACAPs_subset$country == country)
#   plot(ACAPs_subset$date[country_index], ACAPs_subset$pred_mob[country_index], 
#        type = "l", lwd = 2, ylim = c(-100, 10), ylab = "", xlab = country, las = 1)
# }

### Supplementary Code for Cross-Validation ###

# Manual Cross Validation of Individual Countries to Check Overfitting Handled By gmb.step
# country_fold_generator <- function(country_data, number_folds) {
#   fold_indices <- vector("list", number_folds)
#   samples <- length(unique(country_data))/number_folds
#   excluded_indices <- c()
#   for (i in 1:number_folds) {
#     if (length(excluded_indices) == 0) {
#       temp_fold <- sample(1:length(unique(country_data)), samples)
#     } else {
#       temp_fold <- sample(c(1:length(unique(country_data)))[-excluded_indices], samples)
#     }
#     excluded_indices <- c(excluded_indices, temp_fold)
#     fold_indices[[i]] <- temp_fold
#     #print(i)
#   }
#   return(fold_indices)
# }
# 
# number_folds <- 5
# folds <- country_fold_generator(overall$country, number_folds)
# fold_choices <- combinations(n = number_folds, r = number_folds - 1, v = 1:number_folds, repeats.allowed = FALSE)
# cv_outputs <- vector("list", number_folds) 
# countries <- unique(overall$country)
# brts <- list()
# sample_indices_list <- list()
# for (i in 1:number_folds) {
#   folds_to_use <- fold_choices[i, ]
#   sample_indices <- c()
#   for (j in 1:(number_folds - 1)) {
#     sample_indices <- c(sample_indices, folds[[folds_to_use[j]]])
#   }
#   sample_indices_list[[i]] <- sample_indices
#   countries_to_use <- countries[sample_indices]
#   sample_data <- x[overall$country %in% countries_to_use, ]
#   brt_cv <- gbm.step(data = x, gbm.x = c(2:66), gbm.y = 1, family = "gaussian", 
#                      tree.complexity = tree_complexity, learning.rate = learning_rate, 
#                      bag.fraction = bag_fraction, max.trees = max_trees, n.folds = 5)
#   brts[[i]] <- brt_cv
#   observed_train <- x$overall[overall$country %in% countries_to_use]
#   predicted_train <- predict.gbm(brt_cv, 
#                                  x[overall$country %in% countries_to_use, c(2:66)],
#                                  n.trees = brt_cv$gbm.call$best.trees, type = "response")
#   observed_test <- x$overall[!(overall$country %in% countries_to_use)]
#   predicted_test <- predict.gbm(brt_cv, 
#                                 x[!(overall$country %in% countries_to_use),  c(2:66)],
#                                 n.trees = brt_cv$gbm.call$best.trees, type = "response")
#   plot(observed_train, predicted_train, pch = 20)
#   points(observed_test, predicted_test, pch = 20, col = "red")
#   
#   cv_outputs[[i]]$observed_train <- observed_train
#   cv_outputs[[i]]$predicted_train <- predicted_train
#   cv_outputs[[i]]$train_dist <- abs((observed_train - predicted_train))
#   cv_outputs[[i]]$train_ss <- sum((observed_train - predicted_train)^2)
#   cv_outputs[[i]]$observed_test <- observed_test
#   cv_outputs[[i]]$predicted_test <- predicted_test
#   cv_outputs[[i]]$test_dist <- abs((observed_test - predicted_test))
#   cv_outputs[[i]]$test_ss <- sum((observed_test - predicted_test)^2)
# }
# 
# par(mfrow = c(2, 3))
# for (i in 1:number_folds) {
#   observed_train <- cv_outputs[[i]]$observed_train
#   predicted_train <- cv_outputs[[i]]$predicted_train
#   observed_test <- cv_outputs[[i]]$observed_test
#   predicted_test <- cv_outputs[[i]]$predicted_test
#   train_ss <- cv_outputs[[i]]$train_ss
#   test_ss <- cv_outputs[[i]]$test_ss
#   train_dist <- cv_outputs[[i]]$train_dist
#   test_dist <- cv_outputs[[i]]$test_dist
#   print(c(cor(observed_train, predicted_train), cor(observed_test, predicted_test),
#           mean(train_dist), mean(test_dist)))
#   plot(observed_train, predicted_train, pch = 20, cex = 1.5)
#   points(observed_test, predicted_test, pch = 20, col = "red", cex = 1.5)
# }
# 
# par(mfrow = c(5, 5), mar = c(5, 1, 1, 1))
# colours <- c("red", "blue", "orange", "purple", "green")
# countries <- unique(overall$country)
# for (j in 1:length(brts)) {
#   brt <- brts[[j]]
#   sample_indices <- sample_indices_list[[j]]
#   countries_to_use <- countries[-sample_indices]
#   sample_data <- x[overall$country %in% countries_to_use, ]
#   predicted_test <- predict.gbm(brt_cv, 
#                                 x[overall$country %in% countries_to_use,  c(2:59)],
#                                 n.trees = brt_cv$gbm.call$best.trees, type = "response")
#   for (i in 1:length(countries_to_use)) {
#     index <- which(overall$country[overall$country %in% countries_to_use] == countries_to_use[i])
#     predicted <- predicted_test[index]
#     actual_index <- which(overall$country == countries_to_use[i])
#     actual <- overall[actual_index, ]
#     plot(actual$date, actual$overall, ylim = c(-100, 10), pch = 20, cex = 2, ylab = "", xlab = countries_to_use[i])
#     lines(actual$date, predicted, type = "l", ylim = c(-100, 10), lwd = 2, col = colours[j])
#   }
# }