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

# remove imported early death in Philippines
d[which(d$countryterritoryCode=="PHL" & as.Date(d$dateRep) == as.Date("2020-02-02")),]$deaths <- 0

# fix panama's negative deaths
d[which(d$countryterritoryCode=="PAN" & as.Date(d$dateRep)=="2020-06-04"),]$deaths <-
  d[which(d$countryterritoryCode=="PAN" & as.Date(d$dateRep)=="2020-06-04"),]$deaths +
  d[which(d$countryterritoryCode=="PAN" & as.Date(d$dateRep)=="2020-06-03"),]$deaths 
d[which(d$countryterritoryCode=="PAN" & as.Date(d$dateRep)=="2020-06-03"),]$deaths <- 0


# save 
saveRDS(d, "ecdc_all.rds")
