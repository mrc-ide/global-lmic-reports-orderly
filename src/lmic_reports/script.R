orderly_id <- tryCatch(orderly::orderly_run_info()$id,
                       error = function(e) "<id>") # bury this in the html, docx

# prepare fitting first
start <- 10
replicates <- 20
library(dplyr)
library(ggplot2)
library(squire)

## Simulations
ecdc <- readRDS("ecdc_all.rds")
country <- squire::population$country[match(iso3c, squire::population$iso3c)[1]]
df <- ecdc[which(ecdc$countryterritoryCode == iso3c),]

data <- df[,c("dateRep", "deaths", "cases")]
names(data)[1] <- "date"

data$deaths <- rev(cumsum(rev(data$deaths)))
data$cases <- rev(cumsum(rev(data$cases)))
data$date <- as.Date(data$date)

if(data$deaths[1] >= 10) {
  show_sentence <- "" 
} else  {
  show_sentence <-  paste("**N.B. ", country, " is not shown in the following plot as only ", data$deaths[1]," have been reported to date**")
}

pop <- get_population(country)

out <- squire::calibrate(country = country, deaths = max(data$deaths),
                         replicates = replicates, 
                         min_seeding_cases = 1, max_seeding_cases = 5,
                         ICU_bed_capacity = sum(pop$n), 
                         hosp_bed_capacity = sum(pop$n),
                         dt = 0.1)

# prepare reports
rmarkdown::render("index.Rmd", 
                  output_format = c("html_document","pdf_document"), 
                  params = list("out" = out, "replicates" = replicates))


# url_structure: /<iso_date>/<iso_country>/report.html
# url_latest: /latest/<iso_country>/report.html
# get the figures out into a run directory
# figures/
# fig.path or fig.prefix
# pdf/
# can get pdfs to sharepoint easily or latest update by nuking the previous reports
# nightly github release with attached binaries
# rewrite the html output to remove bootstrap
# report.html to index.html