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

# and remove the rows with no data up to the first date that a death was reported
first_report <- which(data$deaths>0)[1]
missing <- which(data$deaths == 0 | is.na(data$deaths))
to_remove <- missing[missing<first_report]
if(length(to_remove) > 0) {
  data <- data[-to_remove,]
}

# dat_0 is just the current date now
date_0 <- date

# get country data
oxford_grt <- readRDS("oxford_grt.rds")

# conduct unmitigated
pop <- squire::get_population(country)

# calibration arguments
reporting_fraction = 1
R0_min = 2.4
R0_max = 4.5
R0_step = 0.2
day_step = 2
int_unique <- squire:::interventions_unique(oxford_grt[[iso3c]], "C")
R0_change <- int_unique$change
date_R0_change <- int_unique$dates_change
date_contact_matrix_set_change <- NULL
squire_model <- explicit_model()
pars_obs <- NULL
n_particles <- 50

null_na <- function(x) {if(is.null(x)) {NA} else {x}}

min_death_date <- data$date[which(data$deaths>0)][1]
last_start_date <- min(as.Date(null_na(date_R0_change[1]))-2, as.Date(null_na(min_death_date))-10, na.rm = TRUE)
first_start_date <- max(as.Date("2020-01-04"),last_start_date - 30, na.rm = TRUE)

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
pos <- which(out$scan_results$mat_log_ll == max(out$scan_results$mat_log_ll), arr.ind = TRUE)

# get tthe R0, betas and times into a data frame
R0 <- out$scan_results$x[pos[1]]
date <- out$scan_results$y[pos[2]]

if(!is.null(date_R0_change)) {
  tt_beta <- c(0, squire:::intervention_dates_for_odin(dates = date_R0_change,
                                                       start_date = date,
                                                       steps_per_day = 1))
} else {
  tt_beta <- 0
}

if(!is.null(R0_change)) {
  R0 <- c(R0, R0 * R0_change)
} else {
  R0 <- R0
}
beta_set <- squire:::beta_est(squire_model = squire_model,
                              model_params = out$scan_results$inputs$model_params,
                              R0 = R0)

df <- data.frame(tt_beta = tt_beta, beta_set = beta_set, date = date + tt_beta)
jsonlite::write_json(jsonlite::toJSON(df), "input_params.json")

