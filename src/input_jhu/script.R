date_0 <- as_date(date)

## -----------------------------------------------------------------------------
## JHU
## -----------------------------------------------------------------------------

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

format_jhu <- function(data, parameter){
  data %>%
    filter(!`Country/Region` %in% c("Diamond Princess", "MS Zaandam", "Summer Olympics 2020", "Winter Olympics 2022")) %>%
    mutate(iso3c = countrycode::countrycode(`Country/Region`, "country.name.en", "iso3c",
                                            custom_match = c(Kosovo = "KSV",
                                                             Micronesia = "FSM"))) %>%
    pivot_longer(cols = matches("\\d+\\/\\d+\\/\\d+"), names_to = "date", values_to = parameter) %>%
    select(iso3c, date, all_of(parameter)) %>%
    mutate(date = mdy(date)) %>%
    #merge provinces into national deaths
    group_by(iso3c, date) %>%
    summarise(
      across(all_of(parameter), ~sum(.x, na.rm = TRUE)),
      .groups = "drop"
    ) %>%
    mutate(
      across(all_of(parameter), ~diff(c(0, .x), na.rm = TRUE))
    )
}

## Get the worldometers data from JHU
jhu_url <- "https://raw.githubusercontent.com/CSSEGISandData/COVID-19/master/csse_covid_19_data/csse_covid_19_time_series/time_series_covid19_deaths_global.csv"
jhu_tf <- download_url(jhu_url)
deaths <- read_csv(jhu_tf) %>%
  format_jhu("deaths")

# now the same for cases
jhu_url <- "https://raw.githubusercontent.com/CSSEGISandData/COVID-19/master/csse_covid_19_data/csse_covid_19_time_series/time_series_covid19_confirmed_global.csv"
jhu_tf <- download_url(jhu_url)
cases <- read_csv(jhu_tf) %>%
  format_jhu("cases")

jhu_df <- full_join(
  deaths, cases,
  by = c("iso3c", "date")
) %>%
  group_by(iso3c) %>%
  complete(date = seq(min(date), max(date), by = 1),
           fill = list(deaths = 0, cases = 0))
#finalise formatting
jhu_df <- jhu_df %>%
  #if less than 0 just set to 0
  mutate(
    across(c(deaths, cases),
           ~if_else(.x < 0, 0, .x))
  ) %>%
  ungroup()

## -----------------------------------------------------------------------------
## Worldometers
## -----------------------------------------------------------------------------
#sometime better than JHU
# function to get data from worldometers
get_country_data <- function(link, iso3c) {
  url <- paste0("https://www.worldometers.info/coronavirus/country/", link)
  html <- xml2::read_html(url)
  scrs <- rvest::html_nodes(html, "script")
  hcs <- grep("Highcharts.chart", unlist(lapply(scrs, as.character)))

  text <- scrs[hcs]

  date_func <- function(dates_d){

    #simplify things by using stringr to reformat dates
    dates_d <- gsub(",", "", str_match_all(dates_d, '\\"\\s*(.*?)\\s*\\"')[[1]][,2])
    #standardise date format
    dates_d <- as.character(mdy(dates_d))

    return(dates_d)

  }

  death_dat <- text[grep("coronavirus-deaths-linear", text)]
  if(length(death_dat) == 0) {

    dates_d <- NA
    deaths <- NA

  } else {

    txt <- xml2::xml_text(death_dat)
    spl <- strsplit(txt, "\n")[[1]]

    dates_d <- spl[grep("categories", spl)][1]
    dates_d <- date_func(dates_d)
    deaths <- spl[grep("data", spl)[1]]
    deaths <- tail(head(strsplit(deaths, ",|\\[|\\]")[[1]],-1),-1)
    deaths <- suppressWarnings(as.numeric(deaths))
    deaths <- c(0,diff(deaths))

  }
  cases_dat <- text[grep("graph-cases-daily", text)]
  txt <- xml2::xml_text(cases_dat)
  spl <- strsplit(txt, "\n")[[1]]

  dates_c <- spl[grep("categories", spl)]
  dates_c <- date_func(dates_c)

  cases <- spl[grep("data", spl)[1]]
  cases <- tail(head(strsplit(cases, ",|\\[|\\]")[[1]],-1),-1)
  cases <- suppressWarnings(as.numeric(cases))

  df <- data.frame("date" = dates_c, "cases" = cases)
  df$deaths <- deaths[match(df$date, dates_d)]
  df$deaths[is.na(df$deaths)] <- 0

  df$iso3c <- iso3c
  df <- df[order(df$date, decreasing = TRUE),]

  return(df)

}

# country names from worldometers
wo <- "https://www.worldometers.info/coronavirus/#countries" %>%
  xml2::read_html() %>%
  rvest::html_nodes(".mt_a")

# create country names and links
worldometer_df <-
  tibble(countries = xml2::xml_text(wo),
         link = gsub("country/|/","",rvest::html_attr(wo, "href"))) %>%
  filter(!countries %in% c("Channel Islands", "Saint Martin", "St. Barth")) %>%
  mutate(iso3c = countrycode::countrycode(countries, "country.name.en", "iso3c",
                                          custom_match = c("CAR" = "CAF",
                                                           Micronesia = "FSM"))
  ) %>%
  select(!countries) %>%
  unique() %>%
  as.list() %>%
  purrr::pmap_dfr(get_country_data) %>%
  mutate(
    across(c(cases, deaths),
           ~if_else(is.na(.x), 0, .x)),
    date = as_date(date)
  )

# OJ's manual cleaning from before, shouldn't matter too much these are no for fitting
worldometer_df <- worldometer_df %>%
  mutate(
    deaths = case_when(
      # FRA
      # -217 deaths day
      date == "2020-05-20" & iso3c == "FRA" ~ 125,
      date == "2020-05-19" & iso3c == "FRA" ~ 186,
      date == "2020-05-18" & iso3c == "FRA" ~  68,
      date == "2020-05-17" & iso3c == "FRA" ~  88,
      date == "2020-05-16" & iso3c == "FRA" ~ 130,
      # CYP
      # -2 deaths day
      date == "2020-04-05" & iso3c == "CYP" ~   0,
      date == "2020-04-04" & iso3c == "CYP" ~   0,
      date == "2020-04-03" & iso3c == "CYP" ~   0,
      # CZE
      # 2 x -1 deaths day
      date == "2020-06-14" & iso3c == "CZE" ~   0,
      date == "2020-06-15" & iso3c == "CZE" ~   0,
      date == "2020-05-19" & iso3c == "CZE" ~   0,
      date == "2020-05-18" & iso3c == "CZE" ~   1,
      # FIN
      # -1 deaths day
      date == "2020-04-07" & iso3c == "FIN" ~   0,
      date == "2020-04-08" & iso3c == "FIN" ~   2,
      # IRL
      # 2 x -deaths day
      date == "2020-06-01" & iso3c == "IRL" ~   0,
      date == "2020-05-31" & iso3c == "IRL" ~   1,
      date == "2020-05-26" & iso3c == "IRL" ~   0,
      date == "2020-05-25" & iso3c == "IRL" ~   2,
      date == "2020-11-24" & iso3c == "IRL" ~   0,
      date == "2020-11-23" & iso3c == "IRL" ~   0,
      # LUX
      # -2 deaths day
      date == "2020-04-15" & iso3c == "LUX" ~   0,
      date == "2020-04-14" & iso3c == "LUX" ~   1,
      # COG
      # few off days
      date == "2020-09-10" & iso3c == "COG" ~   0,
      date == "2020-09-09" & iso3c == "COG" ~   0,
      date == "2020-09-08" & iso3c == "COG" ~   1,
      date == "2020-09-04" & iso3c == "COG" ~   0,
      date == "2020-09-03" & iso3c == "COG" ~   4,
      TRUE ~ deaths
    )
  )

# worldometers is a day ahead of ECDC - so to keep it all aligned
worldometer_df <-
  mutate(
    worldometer_df,
    date = date - 1
  )

## Create a merged data set using WO for certain countries, this way it is easier
#to keep consistent across tasks
combined_df <- jhu_df %>% filter(!(iso3c %in% c("BOL", "ITA", "FRA", "ECU", "CHL", "COD", "ESP", "IRN",
                                                                   "JPN", "GUF","KGZ", "PER", "HKG", "MAC", "TWN",
                                                                   "SDN", "IRL", "TUR", "NPL", "AZE", "BIH", "CRI", "HND",
                                                                   "HTI", "MEX", "SOM", "VEN", "MNG", "LUX", "SAU",
                                                                   "SWE", "USA")))

#now add every country in the worldometer that's not already in JHU
combined_df <- combined_df %>%
  rbind(
    worldometer_df %>% filter(!(iso3c %in% unique(combined_df$iso3c)))
  ) %>%
  arrange(iso3c, date) %>%
  filter(date < date_0)

saveRDS(combined_df, "reported_covid.Rds")
