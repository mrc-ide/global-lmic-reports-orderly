#A function to pre-process the data for model fitting
preprocess_for_fitting <- function(data, iso3c, date_0){

  # get the raw data correct
  data <- data[,c("date", "deaths", "cases")]
  data <- data[order(data$date),]
  data$date <- as.Date(data$date)
  data <- data[data$date <= as.Date(date_0), ]

  deaths_removed <- 0

  # Handle for countries that have eliminated and had reintroduction events
  reintroduction_iso3cs <- c("MMR", "BLZ", "TTO", "BHS", "HKG", "ABW", "GUM", "ISL", "BRB", "MUS", "BRN")
  if (iso3c %in% reintroduction_iso3cs) {
    if(iso3c == "BRN"){
      deaths_removed <- deaths_removed + sum(data$deaths[data$date < as.Date("2020-07-01")])
      data$deaths[data$date < as.Date("2020-07-01")] <- 0
    } else {
      deaths_removed <- deaths_removed + sum(data$deaths[data$date < as.Date("2020-06-01")])
      data$deaths[data$date < as.Date("2020-06-01")] <- 0
    }
  }

  # and remove the rows with no data up to the first date that a death was reported
  first_report <- which(data$deaths>0)[1]
  missing <- which(data$deaths == 0 | is.na(data$deaths))
  to_remove <- missing[missing<first_report]
  if(length(to_remove) > 0) {
    if(length(to_remove) == (nrow(data)-1)) {
      data <- data[-head(to_remove,-1),]
    } else {
      data <- data[-to_remove,]
    }
  }

  #setup week start/week_end to make it compatible with excess mortality likelihood
  data <- data %>%
    transmute(
      week_start = date - 1,
      week_end = date,
      deaths = deaths,
      cases = cases
    )
  return(list(data, deaths_removed))
}
