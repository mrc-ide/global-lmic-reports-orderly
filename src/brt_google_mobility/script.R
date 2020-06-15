## -----------------------------------------------------------------------------
## Step 0: Function Prep
## -----------------------------------------------------------------------------

# Loading Required Libraries
library(gbm); library(dismo); library(conflicted); library(gtools); library(lubridate)
conflict_prefer("select", "dplyr")
conflict_prefer("filter", "dplyr")
conflict_prefer("area", "patchwork")

download_url <- function(url) {
  tryCatch({
    tf <- tempfile()
    code <- download.file(url, tf, mode = "wb")
    if (code != 0) {
      stop("Error downloading file")
    }
  },
  error = function(e) {
    stop(sprintf("Error downloading file '%s': %s, please check %s",
                 url, e$message))
  })
  return(tf)
}

match_clean <- function(a,b, quiet=TRUE){
  a <- gsub("[[:punct:][:space:]]","",tolower(stringi::stri_trans_general(a, "latin-ascii")))
  b <- gsub("[[:punct:][:space:]]","",tolower(stringi::stri_trans_general(b, "latin-ascii")))
  ret <- match(a,b)
  if(sum(is.na(ret)>0)){
    dists <- stringdist::seq_distmatrix(lapply(a,utf8ToInt),lapply(b,utf8ToInt))
    ret[is.na(ret)] <- apply(dists[which(is.na(ret)),,drop=FALSE],1,which.min)
    if(!quiet){
      return(unique(cbind(a,b[ret])))
    }
  }
  return(ret)
}

## -----------------------------------------------------------------------------
## Step 1: Data Loading
## -----------------------------------------------------------------------------

date <- as.Date(date)

# Loading Google Mobility Data
goog_tf <- download_url("https://www.gstatic.com/covid19/mobility/Global_Mobility_Report.csv")
goog <- read.csv(goog_tf, stringsAsFactors = FALSE)
goog$iso3c <- countrycode::countrycode(goog$country_region, "country.name", "iso3c")


# N.B. The date format changed recently in Google and may change again. Look out for errors related to this
mob <- goog %>%
  filter(sub_region_1 == "") %>%
  mutate(overall = 1/4 * retail_and_recreation_percent_change_from_baseline +
           1/4 * grocery_and_pharmacy_percent_change_from_baseline + 
           1/4 * transit_stations_percent_change_from_baseline + 
           1/4 * workplaces_percent_change_from_baseline) %>%
  mutate(date = as.Date(date, format="%Y-%m-%d")) %>%
  select(country_region, iso3c, date, overall)

# Loading World Bank Metadata (Charlie can update to url as needed)
wb_metadata <- read.csv("World_Bank_Country_Metadata.csv", fileEncoding="UTF-8-BOM") %>%
  rename(ISO = country_code) %>%
  select(ISO, income_group, region) %>%
  filter(region != "")

# Loading in ACAPs Data


acap_site <- "http://www.acaps.org/covid19-government-measures-dataset"
xml <- xml2::read_html(acap_site)
url <- rvest::html_attr(rvest::html_nodes(xml, ".file a"), "href")
acap_tf <- download_url(url)

# HTPPS ISSUE - TEMPORARY FIX TO MANUALLY PLACE acaps download into directory
#acap_tf <- "acaps.xlsx"

acap <- readxl::read_excel(acap_tf, progress = FALSE, sheet = "Database")

# country name fixes. We want to use the ISO3C eventually but there are typos...
acap$ISO <- countrycode::countrycode(acap$COUNTRY, "country.name", "iso3c",
                                     custom_match = c("Eswatini"="SWZ", "Micronesia"="FSM"))

ACAPs_measure <- acap %>%
  rename(country = COUNTRY, measure = MEASURE, type = LOG_TYPE, date = DATE_IMPLEMENTED) %>%
  select(ISO, measure, type, date) %>%
  filter(measure != "", type != "") %>%
  mutate(measure = as.numeric(factor(measure)), type = as.numeric(factor(type))) %>%
  mutate(combined = as.factor(paste0("m_",type, "_", measure))) %>%
  mutate(date = as.Date(date)) %>%
  select(ISO, date, combined) %>%
  pivot_wider(names_from = combined, values_from = combined) %>%
  filter(!is.na(date))

## -----------------------------------------------------------------------------
## Step 2: Data Cleaning
## -----------------------------------------------------------------------------

# Renaming Columns of ACAPs Data and Inputting In Number of Each Measure Put In On Each Date
measures <- colnames(ACAPs_measure)[-(1:2)]
for (i in 1:length(measures)) {
  index <- measures[i]
  ACAPs_measure[, index] <- unlist(lapply(ACAPs_measure[[index]], length))
}

# Tracking Cumulative Number of Each Type of Measure Implemented and Creating New Rows for Dates 
# Not Present in Each Country
new_ACAPs_measure <- ACAPs_measure %>%
  group_by(ISO) %>%
  arrange(date) %>%
  mutate_at(vars(-ISO, -date), funs(cumsum(.))) %>%    
  complete(date = seq.Date(min(ACAPs_measure$date), max(ACAPs_measure$date), by = "days")) %>%
  mutate_at(vars( -ISO, -date), funs(replace(., row_number() == 1, 0))) 

# Filling Those Newly Created Rows With the Value In the Row Above (As We're Tracking Cumulative)
column_names <- colnames(new_ACAPs_measure)[-c(1:2)]
new_ACAPs_cat <- fill(new_ACAPs_measure, column_names, .direction = c("down"))

# Joining This Data to the World Bank and Google Mobility Data
overall <- new_ACAPs_cat %>%
  left_join(mob, by = c("ISO" = "iso3c", "date" = "date")) %>% 
  filter(!is.na(overall)) %>%
  left_join(wb_metadata, by = "ISO") %>% 
  select(-country_region)

overall_test <- overall %>%
  ungroup(ISO) %>%
  select(overall, everything(), -ISO, -date)

## -----------------------------------------------------------------------------
## Step 3: BRT
## -----------------------------------------------------------------------------

# Running the BRT 
tree_complexity <- 8
bag_fraction <- 0.5
if (short_run) {
  max_trees <- 30
} else {
max_trees <- 3000
}
learning_rate <- 0.05
x <- as.data.frame(overall_test)
brt <- gbm.step(data = x, 
                gbm.x = 2:ncol(x),
                gbm.y = 1,
                family = "gaussian", 
                tree.complexity = tree_complexity, 
                learning.rate = learning_rate, 
                bag.fraction = bag_fraction, 
                max.trees = max_trees, 
                n.folds = 5, 
                plot.main = FALSE, 
                plot.folds = FALSE)

# predicted <- predict.gbm(brt, x[, c(2:ncol(x))], n.trees = brt$gbm.call$best.trees, type = "response")
# plot(overall$overall, predicted, ylim = c(-100, 20), xlim = c(-100, 20), pch = 20, cex = 2, ylab = "")

# let's create the complete data set and predict where data is missing
output_data <- new_ACAPs_cat %>%
  left_join(mob, by = c("ISO" = "iso3c", "date" = "date")) %>% 
  left_join(wb_metadata, by = "ISO") %>%
  ungroup(ISO) %>%
  select(overall, everything(), -country_region)
predicted <- predict.gbm(brt, output_data[which(is.na(output_data$overall)), c(4:(ncol(output_data)))], n.trees = brt$gbm.call$best.trees, type = "response")
output_data$observed <- !is.na(output_data$overall)
output_data$overall[which(is.na(output_data$overall))] <- predicted
output_data$all_overall <- predict.gbm(brt, output_data[, c(4:(ncol(output_data)))], n.trees = brt$gbm.call$best.trees, type = "response")

## -----------------------------------------------------------------------------
## Step 4: Data Formatting
## -----------------------------------------------------------------------------

## Format the data into the correct form the orderly task
res <- select(output_data, ISO, date, overall, all_overall, observed, income_group) %>% 
  mutate(overall = (overall+100)/100,
         all_overall = (all_overall+100)/100) %>% 
  rename(C = overall,
         C_predict = all_overall,
         iso3c = ISO)

res <- split.data.frame(res, res$iso3c)

# Add in missing countries
nms <- unique(squire::population$iso3c)
res_no <- lapply(nms[!nms %in% names(res)], function(x) {
  return(data.frame())  
})
names(res_no) <- nms[!nms %in% names(res)]

# group them and make sure they are all data frames
res <- append(res, res_no)
res <- lapply(res, function(x) {
  
  if(nrow(x) > 0) {
    return(as.data.frame(x))
  } else {
    return(x)
  }
})

# make sure they all have a pre first switch date
res <- lapply(res,function(x){
  
  if(length(unique(x$C)) == 1) {
    x <- rbind(x[1,],x)
    x$date[1] <- x$date[2]-1
    x[1,"C"] <- 1
  }
  return(x)
  
})

# adjust the unobserved data in the future to be the mean of the last 2 weeks
for(r in seq_along(res)) {
  
  if(any(res[[r]]$observed)) {
    
    # first for prior to the first mobility date use the mean of the first week
    lw <- res[[r]]$C[which(res[[r]]$observed)][1:7]
    res[[r]]$C[1:(which(res[[r]]$observed)[1]-1)] <- mean(lw)
    
    # for the end use the mean of the final 2 weeks
    rw <- tail(res[[r]]$C[which(res[[r]]$observed)], 14)
    rw_end <- tail(which(res[[r]]$observed),1)
    rw_7 <- res[[r]]$C[(rw_end+1):(rw_end+7)] 
    res[[r]]$C[(rw_end+1):nrow(res[[r]])] <- mean(rw)
    
  }
  
}



saveRDS(res, "google_brt.rds")
saveRDS(brt, "google_brt_model.rds")
saveRDS(overall, "overall.rds")

## -----------------------------------------------------------------------------
## Step 5: Data Validation / Plotting
## -----------------------------------------------------------------------------


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