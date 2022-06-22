date_0 <- as.Date(date, "%Y-%m-%d")

##Get most recent fit from github repo
excess_deaths_raw <- read.csv(
  "https://raw.githubusercontent.com/TheEconomist/covid-19-the-economist-global-excess-deaths-model/main/output-data/export_country.csv"
  ) %>%
  mutate(date = as_date(date)) %>%
  filter(date <= date_0)

if(lubridate::week(max(excess_deaths_raw$date)) < lubridate::week(date_0)){
  warning(paste0("Excess mortality data ", lubridate::week(date_0) - lubridate::week(max(excess_deaths_raw$date)), " weeks out of date."))
}

#modify to get a single value
excess_deaths <- excess_deaths_raw %>%
  mutate(date = as.Date(date)) %>%
  arrange(iso3c, date) %>%
  #use real data where poss
  mutate(
    deaths = if_else(
      is.na(daily_excess_deaths),
      estimated_daily_excess_deaths,
      daily_excess_deaths
    ),
    #ensure date is first of week then move to mid week
    date_start = lubridate::floor_date(date, unit = "week"),
    date_end = lubridate::ceiling_date(date, unit = "week")
  ) %>% #summarise incase multiple entries a week
  group_by(iso3c, date_start, date_end) %>%
  summarise(deaths = mean(deaths), .groups = "drop") %>%
  mutate(
    #make it weekly
    deaths = as.integer(deaths*as.numeric(date_end - date_start))
  ) %>%
  #if deaths are negative set to 0
  mutate(
    deaths = if_else(
      deaths < 0,
      as.integer(0),
      deaths
    )
  )

saveRDS(excess_deaths, "excess_deaths.Rds")
