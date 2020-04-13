if (date == "today") {
  date <- Sys.Date()
} else if (date == "yesterday") {
  date <- Sys.Date() - 1
} else {
  date <- as.Date(date, "%Y-%m-%d")
  if (is.na(date)) {
    stop("Date must be provided in ISO format (i.e., YYYY-MM-DD)")
  }
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


keep <- c("Italy", "Germany", "France", "United_Kingdom",
          "United_States_of_America", "South_Korea")

d <- readxl::read_excel("ecdc.xlsx", progress = FALSE)

names(d)[names(d) %in% c("Countries and territories", "countriesAndTerritories")] <- "Region"
saveRDS(d, "ecdc_all.rds")

# original countries to keep
keep <- c("Italy", "Germany", "France", "United_Kingdom",
          "United_States_of_America", "South_Korea")

d <- d[d$Region %in% keep, ]

stopifnot(length(unique(d$Region)) == length(keep))

d$t <- lubridate::decimal_date(d$dateRep)
d <- d[order(d$Region, d$t), ]

saveRDS(d, "ecdc.rds")
