date_0 <- as.Date(date, "%Y-%m-%d")

##Get most recent fit from github repo
excess_deaths_raw <- read.csv(
  "https://raw.githubusercontent.com/TheEconomist/covid-19-the-economist-global-excess-deaths-model/main/output-data/export_country.csv"
  ) %>%
  mutate(date = as.Date(date, origin = "1970-01-01")) %>%
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
    ), #if less than reported covid deaths then replace with that
    deaths = if_else(
      !is.na(daily_covid_deaths) & deaths < daily_covid_deaths,
      daily_covid_deaths,
      deaths
    ),#ensure date is first of week then move to mid week
    week_start = lubridate::floor_date(date, unit = "week"),
    week_end = lubridate::ceiling_date(date, unit = "week")
  ) %>% #summarise incase multiple entries a week
    group_by(iso3c, week_start, week_end) %>%
    summarise(deaths = mean(deaths)) %>%
    ungroup() %>%
  mutate(
    #make it weekly
    deaths = as.integer(deaths*7)
  ) %>%
  #if deaths are negative set to 0
  mutate(
    deaths = if_else(
      deaths < 0,
      as.integer(0),
      deaths
    )
  )

#checks, we print these to a pdf so that we can check them
dir.create("calibration", showWarnings = TRUE)
pdf("calibration/plots.pdf")
#plot known covid deaths against excess mortality
for(country in unique(excess_deaths_raw %>% filter(!is.na(daily_covid_deaths)) %>% pull(iso3c))){
  print(ggplot(excess_deaths_raw %>%
           filter(iso3c == country),
         aes(x = date)
         ) + geom_line(
           aes(y = estimated_daily_excess_deaths),
           colour = "red",
           linetype = "dashed"
         ) + geom_line(
           data = excess_deaths_raw %>%
             filter(iso3c == country, !is.na(daily_excess_deaths)),
           aes(y = daily_excess_deaths),
           colour = "red"
         )+  geom_ribbon(
           aes(ymin = estimated_daily_excess_deaths_ci_50_bot,
               ymax = estimated_daily_excess_deaths_ci_50_top),
           fill = "red",
           alpha = 0.1
         )  +  geom_ribbon(
           aes(ymin = estimated_daily_excess_deaths_ci_90_bot,
               ymax = estimated_daily_excess_deaths_ci_90_top),
           fill = "red",
           alpha = 0.1
         )  +  geom_ribbon(
           aes(ymin = estimated_daily_excess_deaths_ci_95_bot,
               ymax = estimated_daily_excess_deaths_ci_95_top),
           fill = "red",
           alpha = 0.1
         ) + geom_line(
           aes(y = daily_covid_deaths),
           colour = "black"
         ) + labs(
           x = "Week",
           y = "Daily average for the week",
           title = country
         ))
}
dev.off()

saveRDS(excess_deaths, "excess_deaths.rds")
