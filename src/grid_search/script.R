orderly_id <- tryCatch(orderly::orderly_run_info()$id,
                       error = function(e) "<id>") # bury this in the html, docx

set.seed(123)
date <- as.Date(date)
# saveRDS("unfinished", paste0("/home/oj/GoogleDrive/AcademicWork/covid/githubs/global-lmic-reports-orderly/scripts/",iso3c,".rds"))

# prepare fitting first
start <- 10
replicates <- 50

## Get the ECDC data
ecdc <- readRDS("ecdc_all.rds")
country <- squire::population$country[match(iso3c, squire::population$iso3c)[1]]
df <- ecdc[which(ecdc$countryterritoryCode == iso3c),]

# get the raw data correct
data <- df[,c("dateRep", "deaths", "cases")]
names(data)[1] <- "date"
data <- data[order(data$date),]
data$date <- as.Date(data$date)

# what was the date being calibrated to
if(max(data$deaths) < 10) {
  date_0 <- date
} else {
  #date_0 <- data$date[tail(which(data$deaths==max(data$deaths)),1)]
  date_0 <- date
}

# get country data
oxford_grt <- readRDS("oxford_grt.rds")

# conduct unmitigated
pop <- squire::get_population(country)

# calibration arguments
reporting_fraction = 1
R0_min = 2.4
R0_max = 3.4
R0_step = 0.1
first_start_date = "2020-01-30"
last_start_date = "2020-02-12"
day_step = 2
int_unique <- squire:::interventions_unique(oxford_grt[[iso3c]], "C")
R0_change <- int_unique$change
date_R0_change <- as.Date(as.character(int_unique$dates_change), "%Y%m%d")
date_contact_matrix_set_change <- NULL
squire_model <- explicit_model()
pars_obs <- NULL
n_particles <- 50


future::plan(future::multiprocess())

out <- squire::calibrate_particle(
  data = data,
  R0_min = R0_min,
  R0_max = R0_max,
  R0_step = R0_step,
  first_start_date = first_start_date,
  last_start_date = last_start_date,
  day_step = day_step,
  squire_model = squire_model,
  pars_obs = pars_obs,
  n_particles = n_particles,
  reporting_fraction = reporting_fraction,
  R0_change = R0_change,
  date_R0_change = date_R0_change,
  replicates = replicates,
  country = country,
  forecast = 120
)

saveRDS(out, "grid_out.rds")

## and save the info for the interface
pos <- which(out$scan_results$renorm_mat_LL == max(out$scan_results$renorm_mat_LL), arr.ind = TRUE)

R0 <- out$scan_results$x[pos[1]]
date <- out$scan_results$y[pos[2]]
oxford

