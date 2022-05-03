
date <- as.Date(date, "%Y-%m-%d")
if(date <= as.Date("2020-12-12")) {

if (is.na(date)) {
  stop("Date must be provided in ISO format (i.e., YYYY-MM-DD)")
}

fmt <- "https://www.ecdc.europa.eu/sites/default/files/documents/COVID-19-geographic-disbtribution-worldwide-%s.xlsx"
date_iso <- as.character(date)
url <- sprintf(fmt, date_iso)

url_page <- "https://www.ecdc.europa.eu/en/publications-data/download-todays-data-geographic-distribution-covid-19-cases-worldwide"
tryCatch({
  code <- download.file(url, "ecdc.xlsx", mode = "wb")
  if (code != 0) {
    stop("Error downloading file")
  }
},
error = function(e) {
  stop(sprintf("Error downloading file '%s': %s, please check %s",
               url, e$message, url_page))
})

## -----------------------------------------------------------------------------
## ECDC
## -----------------------------------------------------------------------------

d <- readxl::read_excel("ecdc.xlsx", progress = FALSE)
names(d)[1] <- "dateRep"

names(d)[names(d) %in% c("Countries and territories", "countriesAndTerritories")] <- "Region"

# sort out the date malarkey that happened on the 16th May
if(is.na(as.Date(d$dateRep[1], "%Y-%m-%d"))) {
  d$dateRep <- as.Date(d$dateRep, format = "%d/%m/%Y")
}

# Data corrections

# fix negative deaths
d$deaths[which(d$deaths<0)] <- 0
d$deaths[which(is.na(d$deaths))] <- 0

# date issues
suppressWarnings(md <- lapply(unique(d$countryterritoryCode), function(x){
  max(as.Date(d$dateRep[d$countryterritoryCode==x]),na.rm = TRUE)
  }))

to_fix <- na.omit(unique(d$countryterritoryCode)[which(unlist(md) != as.numeric(date))])

# spain seems to have stopped reporting now too...
if(length(to_fix) > 0) {

  for(i in seq_along(to_fix)) {
  df_esp <- d[which(d$countryterritoryCode==to_fix[i]),][1,]
  df_esp$dateRep = max(d$dateRep)
  df_esp$day = as.numeric(format(as.Date(date), "%d"))
  df_esp$month = as.numeric(format(as.Date(date), "%m"))
  df_esp$year = as.numeric(format(as.Date(date), "%Y"))
  df_esp$cases = 0
  df_esp$deaths = 0
  d <- rbind(df_esp, d)
  }

}
ecdc <- d

# save
saveRDS(ecdc, "ecdc_all.rds")

} else {
  saveRDS(data.frame(), "ecdc_all.rds")
  file.create("ecdc.xlsx")
}

## -----------------------------------------------------------------------------
## JHU
## -----------------------------------------------------------------------------

#### AND let's get the JHU as well/instead as looks liek it is less susceptible to blips

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


## Get the worldometers data from JHU
jhu_url <- "https://raw.githubusercontent.com/CSSEGISandData/COVID-19/master/csse_covid_19_data/csse_covid_19_time_series/time_series_covid19_deaths_global.csv"
jhu_tf <- download_url(jhu_url)
data <- read.csv(jhu_tf)

# format into the same style as ecdc so easy to swap back and forth
data$countryterritoryCode <- suppressWarnings(countrycode::countrycode(data$Country.Region, "country.name.en", "iso3c",
                                                      custom_match = c(Kosovo = "KSV")))
data <- data %>% tidyr::pivot_longer(matches("X\\d"))
names(data) <- c("", "Region","lat","lon","countryterritoryCode","date","deaths")

data <- data[,c("date","deaths","countryterritoryCode","Region")]
data$date <- as.Date(data$date, format = "X%m.%d.%y")

# and into daily deaths nationally
data <- group_by(data, date, countryterritoryCode, Region) %>%
  summarise(deaths = sum(deaths, na.rm = TRUE))
data <- group_by(data, countryterritoryCode, Region) %>%
  mutate(deaths = c(0, diff(deaths)))

# now the same for cases
jhu_url <- "https://raw.githubusercontent.com/CSSEGISandData/COVID-19/master/csse_covid_19_data/csse_covid_19_time_series/time_series_covid19_confirmed_global.csv"
jhu_tf <- download_url(jhu_url)
cases <- read.csv(jhu_tf)

# format into the same style as ecdc so easy to swap back and forth
cases$countryterritoryCode <- suppressWarnings(countrycode::countrycode(cases$Country.Region, "country.name.en", "iso3c",
                                                                       custom_match = c(Kosovo = "KSV")))
cases <- cases %>% tidyr::pivot_longer(matches("X\\d"))
names(cases) <- c("", "Region","lat","lon","countryterritoryCode","date","cases")

cases <- cases[,c("date","cases","countryterritoryCode","Region")]
cases$date <- as.Date(cases$date, format = "X%m.%d.%y")

# and into daily cases nationally
cases <- group_by(cases, date, countryterritoryCode, Region) %>%
  summarise(cases = sum(cases, na.rm = TRUE))

cases <- group_by(cases, countryterritoryCode, Region) %>%
  mutate(cases = c(0, diff(cases)))

jhu_data <- left_join(data, cases, by = c("date", "countryterritoryCode", "Region"))
jhu_data$dateRep <- jhu_data$date

if(sum(jhu_data$deaths)>0){
  jhu_data$deaths[jhu_data$deaths < 0] <- 0
}

#spread out tanzania spike
if(jhu_data %>%
   filter(countryterritoryCode == "TZA", dateRep == "2021-10-01") %>%
   pull(deaths) == 669){
  #split up over previous and have it linearly increase
  new_values <- round(seq(0, 2/14, length.out = 14)*669)
  jhu_data <- jhu_data %>%
    mutate(
      deaths = if_else(
        countryterritoryCode == "TZA" & dateRep == "2021-10-01",
        0,
        deaths
      )
    )
  jhu_data[jhu_data$countryterritoryCode == "TZA" & (jhu_data$dateRep %in%
             (as.Date("2021-10-01") - 0:13)) &
             !is.na(jhu_data$countryterritoryCode),
           "deaths"] <-jhu_data[jhu_data$countryterritoryCode == "TZA" &
                                  (jhu_data$dateRep %in% (as.Date("2021-10-01") - 0:13)) &
                                  !is.na(jhu_data$countryterritoryCode),"deaths"] +
    new_values
}

# save
saveRDS(jhu_data, "jhu_all.rds")

## -----------------------------------------------------------------------------
## Worldometers
## -----------------------------------------------------------------------------

# function to get data from worldometers
get_country_data <- function(link, iso3c, name) {

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

  df <- data.frame("dateRep" = dates_c, "cases" = cases)
  df$deaths <- deaths[match(df$dateRep, dates_d)]
  df$deaths[is.na(df$deaths)] <- 0

  df$countryterritoryCode <- iso3c
  df$Region <- name
  df <- df[order(df$dateRep, decreasing = TRUE),]

  return(df)

}

# country names from worldometers
wo <- "https://www.worldometers.info/coronavirus/#countries"
wo <- xml2::read_html(wo) %>% rvest::html_nodes(".mt_a")

# create country names and links
countries <- xml2::xml_text(wo)
links <- gsub("country/|/","",rvest::html_attr(wo, "href"))
iso3cs <- countrycode::countrycode(countries, "country.name.en", "iso3c",
                                   custom_match = c("CAR" = "CAF"))

# df of args to run
df_args <- data.frame(countries = countries, links = links, iso3cs = iso3cs)
df_args <- unique(df_args)
df_args <- na.omit(df_args)

# loop over data needs
dats <- lapply(seq_along(df_args$countries), function(i) {
  get_country_data(df_args$links[i], df_args$iso3cs[i], df_args$countries[i])
})
df <- do.call(rbind, dats)

# and fill in leading NAs
df$cases[is.na(df$cases)] <- 0
df$deaths[is.na(df$deaths)] <- 0

# AND handling their peculiar negative deaths based on comparison against ECDC and manually cleaning :)

# FRA
# -217 deaths day
df$deaths[df$dateRep == as.Date("2020-05-20") & df$countryterritoryCode == "FRA"] <- 125
df$deaths[df$dateRep == as.Date("2020-05-19") & df$countryterritoryCode == "FRA"] <- 186
df$deaths[df$dateRep == as.Date("2020-05-18") & df$countryterritoryCode == "FRA"] <- 68
df$deaths[df$dateRep == as.Date("2020-05-17") & df$countryterritoryCode == "FRA"] <- 88
df$deaths[df$dateRep == as.Date("2020-05-16") & df$countryterritoryCode == "FRA"] <- 130

# CYP
# -2 deaths day
df$deaths[df$dateRep == as.Date("2020-04-05") & df$countryterritoryCode == "CYP"] <- 0
df$deaths[df$dateRep == as.Date("2020-04-04") & df$countryterritoryCode == "CYP"] <- 0
df$deaths[df$dateRep == as.Date("2020-04-03") & df$countryterritoryCode == "CYP"] <- 0

# CZE
# 2 x -1 deaths day
df$deaths[df$dateRep == as.Date("2020-06-14") & df$countryterritoryCode == "CZE"] <- 0
df$deaths[df$dateRep == as.Date("2020-06-15") & df$countryterritoryCode == "CZE"] <- 0
df$deaths[df$dateRep == as.Date("2020-05-19") & df$countryterritoryCode == "CZE"] <- 0
df$deaths[df$dateRep == as.Date("2020-05-18") & df$countryterritoryCode == "CZE"] <- 1

# FIN
# -1 deaths day
df$deaths[df$dateRep == as.Date("2020-04-07") & df$countryterritoryCode == "FIN"] <- 0
df$deaths[df$dateRep == as.Date("2020-04-08") & df$countryterritoryCode == "FIN"] <- 2

# IRL
# 2 x -deaths day
df$deaths[df$dateRep == as.Date("2020-06-01") & df$countryterritoryCode == "IRL"] <- 0
df$deaths[df$dateRep == as.Date("2020-05-31") & df$countryterritoryCode == "IRL"] <- 1
df$deaths[df$dateRep == as.Date("2020-05-26") & df$countryterritoryCode == "IRL"] <- 0
df$deaths[df$dateRep == as.Date("2020-05-25") & df$countryterritoryCode == "IRL"] <- 2
df$deaths[df$dateRep == as.Date("2020-11-24") & df$countryterritoryCode == "IRL"] <- 0
df$deaths[df$dateRep == as.Date("2020-11-23") & df$countryterritoryCode == "IRL"] <- 0

# LUX
# -2 deaths day
df$deaths[df$dateRep == as.Date("2020-04-15") & df$countryterritoryCode == "LUX"] <- 0
df$deaths[df$dateRep == as.Date("2020-04-14") & df$countryterritoryCode == "LUX"] <- 1

# COG
# few off days
df$deaths[df$dateRep == as.Date("2020-09-10") & df$countryterritoryCode == "COG"] <- 0
df$deaths[df$dateRep == as.Date("2020-09-09") & df$countryterritoryCode == "COG"] <- 0
df$deaths[df$dateRep == as.Date("2020-09-08") & df$countryterritoryCode == "COG"] <- 1
df$deaths[df$dateRep == as.Date("2020-09-04") & df$countryterritoryCode == "COG"] <- 0
df$deaths[df$dateRep == as.Date("2020-09-03") & df$countryterritoryCode == "COG"] <- 4

#adjust the strange 0 deaths issue with Azerbaijan
if(df$deaths[df$dateRep == as.Date("2021-09-28") & df$countryterritoryCode == "AZE"] == 0){
  #check if its actually an issue
  df$deaths[df$dateRep == as.Date("2021-09-28") & df$countryterritoryCode == "AZE"] <-
    ceiling((df$deaths[df$dateRep == as.Date("2021-09-29") & df$countryterritoryCode == "AZE"])/2)
  df$deaths[df$dateRep == as.Date("2021-09-29") & df$countryterritoryCode == "AZE"] <-
    df$deaths[df$dateRep == as.Date("2021-09-29") & df$countryterritoryCode == "AZE"] -
    df$deaths[df$dateRep == as.Date("2021-09-28") & df$countryterritoryCode == "AZE"]
}

# worldometers is a day ahead of ECDC - so to keep it all aligned
df$dateRep <- as.Date(df$dateRep) - 1
saveRDS(df, "worldometers_all.rds")

## Create a merged data set using WO for certain countries, this way it is easier
#to keep consistent across tasks
countries_not_in_jhu <- setdiff(unique(df$countryterritoryCode), unique(df$jhu_data))
combined_data <- rbind(
  jhu_data %>% filter(!(countryterritoryCode %in% c("BOL", "ITA", "FRA", "ECU", "CHL", "COD", "ESP", "IRN",
                                                "JPN", "GUF","KGZ", "PER", "HKG", "MAC", "TWN",
                                                "SDN", "IRL", "TUR", "NPL", "AZE", "BIH", "CRI", "HND",
                                                "HTI", "MEX", "SOM", "VEN", "MNG", "LUX", "SAU",
                                                "SWE", "USA"))) %>%
    select(!date),
  df %>% filter((countryterritoryCode %in% c("BOL", "ITA", "FRA", "ECU", "CHL", "COD", "ESP", "IRN",
                                            "JPN", "GUF","KGZ", "PER", "HKG", "MAC", "TWN",
                                            "SDN", "IRL", "TUR", "NPL", "AZE", "BIH", "CRI", "HND",
                                            "HTI", "MEX", "SOM", "VEN", "MNG",  "LUX", "SAU",
                                            "SWE", "USA")) | countryterritoryCode %in% countries_not_in_jhu) %>%
    mutate(
      dateRep = as.Date(dateRep)
    )
)

saveRDS(combined_data, "combined_data.Rds")

## -----------------------------------------------------------------------------
## OWID
## -----------------------------------------------------------------------------

# And let's add owid here
owid_url <- "https://covid.ourworldindata.org/data/owid-covid-data.csv"
owid_tf <- download_url(owid_url)
data <- read.csv(owid_tf)

# just to align
names(data)[which(names(data) == "iso_code")] <- "countryterritoryCode"

# save data out
saveRDS(data, "owid.rds")

## -----------------------------------------------------------------------------
## OWID Vaccines by manufacturer and agreed sale
## -----------------------------------------------------------------------------

vacc_loc_url <- "https://raw.githubusercontent.com/owid/covid-19-data/master/public/data/vaccinations/locations.csv"
vacc_loc <- download_url(vacc_loc_url)
vacc_loc <- read.csv(vacc_loc)
vacc_loc$vaccine_types <- strsplit(vacc_loc$vaccines, ", ")
names(data)[which(names(data) == "iso_code")] <- "countryterritoryCode"
saveRDS(vacc_loc, "vaccine_agreements.rds")

vacc_by_type_url <- "https://raw.githubusercontent.com/owid/covid-19-data/master/public/data/vaccinations/vaccinations-by-manufacturer.csv"
vacc_by_type <- download_url(vacc_by_type_url)
vacc_by_type <- read.csv(vacc_by_type)
vacc_by_type$countryterritoryCode <- countrycode::countrycode(
  vacc_by_type$location, "country.name.en", "iso3c"
)
saveRDS(vacc_by_type, "vaccine_doses_by_manufacturer.rds")

who_vacc_url <- "https://covid19.who.int/who-data/vaccination-data.csv"
who_vacc <- download_url(who_vacc_url)
who_vacc <- read.csv(who_vacc)
who_vacc$countryterritoryCode <- who_vacc$ISO3
who_vacc$DATE_UPDATED <- as.Date(who_vacc$DATE_UPDATED)

saveRDS(who_vacc, "who_vacc.rds")


who_vacc_meta_url <- "https://covid19.who.int/who-data/vaccination-metadata.csv"
who_vacc_meta <- download_url(who_vacc_meta_url)
who_vacc_meta <- read.csv(who_vacc_meta)
who_vacc_meta$countryterritoryCode <- who_vacc_meta$ISO3
who_vacc_meta$START_DATE[who_vacc_meta$START_DATE == ""] <- NA
who_vacc_meta$START_DATE <- as.Date(who_vacc_meta$START_DATE)
saveRDS(who_vacc_meta, "who_vacc_meta.rds")
