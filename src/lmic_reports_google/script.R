orderly_id <- tryCatch(orderly::orderly_run_info()$id,
                       error = function(e) "<id>") # bury this in the html, docx

version_min <- "0.4.8"
if(packageVersion("squire") < version_min) {
  stop("squire needs to be updated to at least", version_min)
}

## -----------------------------------------------------------------------------
## Step 1: Incoming Date
## -----------------------------------------------------------------------------
system(paste0("echo ",iso3c))
set.seed(123)
date <- as.Date(date)
reports <- reports_day(as.character(date-1))

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
  if(length(to_remove) == (nrow(data)-1)) {
    data <- data[-head(to_remove,-1),]
  } else {
    data <- data[-to_remove,]
  }
}

# dat_0 is just the current date now
date_0 <- date

# get country data
# interventions <- readRDS("oxford_grt.rds")
interventions <- readRDS("google_brt.rds")

# conduct unmitigated
pop <- squire::get_population(country)

## -----------------------------------------------------------------------------
## Step 2: Times for the interace based on the deterministic
## -----------------------------------------------------------------------------

## -----------------------------------------------------------------------------
## Step 2a: Sourcing previous fits to limit grid search
## -----------------------------------------------------------------------------

# what is the date of first death
null_na <- function(x) {if(is.null(x)) {NA} else {x}}
min_death_date <- data$date[which(data$deaths>0)][1]

# calibration arguments
reporting_fraction = 1
int_unique <- squire:::interventions_unique(interventions[[iso3c]], "C")
R0_change <- int_unique$change
date_R0_change <- int_unique$dates_change
date_contact_matrix_set_change <- NULL
squire_model <- explicit_model()
pars_obs <- NULL
day_step = 1
R0_step = 0.2
n_particles <- 100
replicates <- 100

# 1. Do we have a previous report for this country
if (!is.null(reports) && (iso3c %in% reports$country)) {
  
  out <- file.path(here::here(), "archive", "lmic_reports", reports$id[which(reports$country==iso3c)], "grid_out.rds")
  out <- readRDS(out)
  
  if (!is.nulll(out$scan_results$z)) {
    
    # recreate the grids
    x_grid <- array(out$scan_results$x, dim(out$scan_results$renorm_mat_LL))
    y_grid <- array(mapply(rep, out$scan_results$y, length(out$scan_results$x)), dim(out$scan_results$renorm_mat_LL))
    if (!is.null(out$scan_results$z)) {
      z_grid <- array(mapply(rep, out$scan_results$z, length(out$scan_results$x)*length(out$scan_results$y)), dim(out$scan_results$renorm_mat_LL))
    }
    
    # first get the sorted density
    ord <- order(out$scan_results$renorm_mat_LL, decreasing = TRUE)
    cum <- cumsum(out$scan_results$renorm_mat_LL[ord])
    ninety <- which(cum > 0.9)[1]
    
    # get the range from this for R0 and grow it by 0.2
    R0_max <- max(x_grid[ord[seq_len(ninety)]]) + 0.2
    R0_min <- min(x_grid[ord[seq_len(ninety)]]) - 0.2
    
    # get the range for dates and grow it by 3 days
    last_start_date <- max(y_grid[ord[seq_len(ninety)]])
    first_start_date <- min(y_grid[ord[seq_len(ninety)]])
    last_start_date <- as.Date(out$scan_results$y[match(last_start_date, as.numeric(out$scan_results$y))]) + 3
    first_start_date <- as.Date(out$scan_results$y[match(first_start_date, as.numeric(out$scan_results$y))]) -3
    
    # adust the dates so they are compliant with the data
    last_start_date <- min(c(last_start_date, as.Date(null_na(min_death_date))-10), na.rm = TRUE)
    first_start_date <- max(as.Date("2020-01-04"), first_start_date, na.rm = TRUE)
    
    # get the range for Meff
    Meff_max <- max(z_grid[ord[seq_len(ninety)]]) + 0.2
    Meff_min <- min(z_grid[ord[seq_len(ninety)]]) - 0.2
    Meff_step <- 0.1
    
  } else {
    
    # Defualts if no previous data
    R0_min = 2.0
    R0_max = 5.6
    Meff_min = 0.5
    Meff_max = 2
    Meff_step = 0.3
    last_start_date <- as.Date(null_na(min_death_date))-10
    first_start_date <- max(as.Date("2020-01-04"),last_start_date - 30, na.rm = TRUE)
    
  }
  
} else {
  
  # Defualts if no previous data
  R0_min = 2.0
  R0_max = 5.6
  Meff_min = 0.5
  Meff_max = 2
  Meff_step = 0.3
  last_start_date <- as.Date(null_na(min_death_date))-20
  first_start_date <- max(as.Date("2020-01-04"),last_start_date - 30, na.rm = TRUE)
  
}

out_det <- squire::calibrate(
  data = data,
  R0_min = R0_min,
  R0_max = R0_max,
  R0_step = 0.05,
  R0_prior = list("func" = dnorm, args = list("mean"= 3.5, "sd"= 3, "log = TRUE")),
  Meff_min = Meff_min,
  Meff_max = 1,
  Meff_step = 0.01,
  first_start_date = first_start_date,
  last_start_date = last_start_date,
  day_step = 1,
  squire_model = squire:::deterministic_model(),
  pars_obs = pars_obs,
  n_particles = 2,
  reporting_fraction = reporting_fraction,
  R0_change = R0_change,
  date_R0_change = date_R0_change,
  replicates = replicates,
  country = country,
  forecast = 28,
)

## and save the info for the interface
pos <- which(out_det$scan_results$mat_log_ll == max(out_det$scan_results$mat_log_ll), arr.ind = TRUE)

# get tthe R0, betas and times into a data frame
R0 <- out_det$scan_results$x[pos[1]]
start_date <- out_det$scan_results$y[pos[2]]
Meff <- out_det$scan_results$z[pos[3]]

if(!is.null(date_R0_change)) {
  start_date <- min(start_date, date_R0_change-1)
  tt_beta <- squire:::intervention_dates_for_odin(dates = date_R0_change,
                                                  change = R0_change,
                                                  start_date = start_date,
                                                  steps_per_day = 1)
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
                              R0 = R0*Meff)

df <- data.frame(tt_beta = c(0,tt_beta$tt), beta_set = beta_set, date = start_date + c(0,tt_beta$tt), Meff = Meff)
writeLines(jsonlite::toJSON(df,pretty = TRUE), "input_params.json")


## -----------------------------------------------------------------------------
## Step 3: Particle Filter
## -----------------------------------------------------------------------------


# future::plan(future::multiprocess())
out <- squire::calibrate(
  data = data,
  R0_min = R0_min,
  R0_max = R0_max,
  R0_step = R0_step,
  Meff_min = Meff_min,
  Meff_max = Meff_max,
  Meff_step = Meff_step,
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
  forecast = 28
)

saveRDS(out, "grid_out.rds")

## summarise what we have
prob1 <- plot(out$scan_results, what="probability", log = FALSE, show = c(1,2))
prob2 <- plot(out$scan_results, what="probability", log = FALSE, show = c(1,3))
prob3 <- plot(out$scan_results, what="probability", log = FALSE, show = c(2,3))

index <- squire:::odin_index(out$model)
forecast <- 14

d <- plot(out, "deaths", date_0 = date, x_var = "date")
ymax <- max(out$scan_results$inputs$data$deaths, d$layers[[1]]$data$ymax)
d <- d + geom_point(data = out$scan_results$inputs$data, 
                    mapping = aes(x=date,y=deaths), inherit.aes = FALSE) + 
  scale_x_date(limits = c(min(data$date),date+forecast)) +
  scale_y_continuous(limits = c(0,ymax)) + 
  geom_vline(xintercept = date, linetype = "dashed") +
  ylab("Deaths") + 
  xlab("") +
  theme(legend.position = "none")

#intervention <- intervention_plot(interventions[[iso3c]], date)
intervention <- intervention_plot_google(interventions[[iso3c]], date)

title <- cowplot::ggdraw() + 
  cowplot::draw_label(
    country,
    fontface = 'bold',
    x = 0.5
  )

line <- ggplot() + cowplot::draw_line(x = 0:10,y=1) + 
  theme(panel.background = element_blank(),
        axis.title = element_blank(), 
        axis.text = element_blank(), 
        axis.ticks = element_blank())

top_row <- cowplot::plot_grid(prob1, prob2, prob3, ncol=3)

pdf("fitting.pdf",width = 6,height = 10)
print(cowplot::plot_grid(title,line,top_row,intervention,d,ncol=1,rel_heights = c(0.1,0.1,0.8,0.6,1)))
dev.off()



## -----------------------------------------------------------------------------
## Step 4: Scenarios
## -----------------------------------------------------------------------------

# conduct scnearios
mit <- squire::projections(out, R0_change = 0.5, tt_R0 = 0)

if(!is.null(out$interventions$R0_change)) {
  rev_change <- 1-((1-tail(out$interventions$R0_change,1))/2)
} else {
  rev_change <- 1
}
rev <- squire::projections(out, R0 = mean(out$replicate_parameters$R0)*rev_change, tt_R0 = 0)

r_list <- list(out, mit, rev)
o_list <- lapply(r_list, squire::format_output,
                 var_select = c("infections","deaths","hospital_demand","ICU_demand", "D"),
                 date_0 = date_0)

## -----------------------------------------------------------------------------
## Step 5: Report
## -----------------------------------------------------------------------------

# get data in correct format for plotting
df <- ecdc[which(ecdc$countryterritoryCode == iso3c),]

# get the raw data correct
data <- df[,c("dateRep", "deaths", "cases")]
names(data)[1] <- "date"
data$daily_deaths <- data$deaths
data$daily_cases <- data$cases
data$deaths <- rev(cumsum(rev(data$deaths)))
data$cases <- rev(cumsum(rev(data$cases)))
data$date <- as.Date(data$date)

# prepare reports
rmarkdown::render("index.Rmd", 
                  output_format = c("html_document","pdf_document"), 
                  params = list("r_list" = r_list,
                                "o_list" = o_list,
                                "replicates" = replicates, 
                                "data" = data,
                                "date_0" = date_0,
                                "country" = country),
                  output_options = list(pandoc_args = paste0("--metadata=title:",country," COVID-19 report")))

data_sum <- lapply(o_list, function(pd){
  
  # remove any NA rows (due to different start dates)
  if(sum(is.na(pd$t) | is.na(pd$y))>0) {
    pd <- pd[-which(is.na(pd$t) | is.na(pd$y)),]
  }
  
  # Format summary data
  pds <- pd %>%
    dplyr::group_by(.data$date, .data$compartment) %>%
    dplyr::summarise(y_025 = stats::quantile(.data$y, 0.025),
                     y_25 = stats::quantile(.data$y, 0.25),
                     y_median = median(.data$y),
                     y_mean = mean(.data$y),
                     y_75 = stats::quantile(.data$y, 0.75),
                     y_975 = stats::quantile(.data$y, 0.975))
  
  return(as.data.frame(pds, stringsAsFactors = FALSE))
})
data_sum[[1]]$scenario <- "Maintain Status Quo"
data_sum[[2]]$scenario <- "Additional 50% Reduction"
data_sum[[3]]$scenario <- "Relax Interventions 50%"
data_sum <- do.call(rbind, data_sum)
data_sum$country <- country
data_sum$iso3c <- iso3c
data_sum$report_date <- date
data_sum <- data_sum[data_sum$compartment != "D",]
write.csv(data_sum, "projections.csv", row.names = FALSE, quote = FALSE)

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