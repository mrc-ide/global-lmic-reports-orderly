date <- as.Date(date, "%Y-%m-%d")
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

# save 
saveRDS(d, "ecdc_all.rds")


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

# save 
saveRDS(d, "jhu_all.rds")
