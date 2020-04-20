orderly_id <- tryCatch(orderly::orderly_run_info()$id,
                       error = function(e) "<id>") # bury this in the html, docx

set.seed(123)
date <- as.Date(date)
# saveRDS("unfinished", paste0("/home/oj/GoogleDrive/AcademicWork/covid/githubs/global-lmic-reports-orderly/scripts/",iso3c,".rds"))

# prepare fitting first
start <- 10
replicates <- 50

library(dplyr)
library(ggplot2)
library(squire)

## Simulations
ecdc <- readRDS("ecdc_all.rds")
country <- squire::population$country[match(iso3c, squire::population$iso3c)[1]]
df <- ecdc[which(ecdc$countryterritoryCode == iso3c),]

# get the raw data correct
data <- df[,c("dateRep", "deaths", "cases")]
names(data)[1] <- "date"
data$deaths <- rev(cumsum(rev(data$deaths)))
data$cases <- rev(cumsum(rev(data$cases)))
data$date <- as.Date(data$date)

# what was the date being calibrated to
if(max(data$deaths) < 10) {
date_0 <- date
} else {
date_0 <- data$date[tail(which(data$deaths==max(data$deaths)),1)]
}

# get country data
oxford_grt <- readRDS("oxford_grt.rds")

# conduct unmitigated
pop <- squire::get_population(country)
out <- squire::calibrate(country = country, deaths = max(data$deaths),
                         replicates = replicates,
                         min_seeding_cases = 1, max_seeding_cases = 5,
                         ICU_bed_capacity = sum(pop$n),
                         hosp_bed_capacity = sum(pop$n),
                         tt_R0 = oxford_grt[[iso3c]]$tt_R0,
                         R0 = oxford_grt[[iso3c]]$R0,
                         dt = 0.1)
# conduct scnearios
mit <- squire::projections(out, R0_change = 0.5, tt_R0 = 0)
r_list <- list(out, mit)
o_list <- lapply(r_list, squire::format_output,
                 var_select = c("infections","deaths","hospital_demand","ICU_demand", "D"),
                 date_0 = date_0)




# prepare reports
rmarkdown::render("index.Rmd", 
                  output_format = c("html_document","pdf_document"), 
                  params = list("r_list" = r_list,
                                "o_list" = o_list,
                                "replicates" = replicates, 
                                "data" = data,
                                "date_0" = date_0,
                                "country" = country),
                  output_options = list(pandoc_args = paste0("--metadata=title:\"",country," COVID-19 report\"")))

# saveRDS("finished", paste0("/home/oj/GoogleDrive/AcademicWork/covid/githubs/global-lmic-reports-orderly/scripts/",iso3c,".rds"))

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