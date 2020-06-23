orderly_id <- tryCatch(orderly::orderly_run_info()$id,
                       error = function(e) "<id>") # bury this in the html, docx

version_min <- "0.4.22"
if(packageVersion("squire") < version_min) {
  stop("squire needs to be updated to at least ", version_min)
}

## -----------------------------------------------------------------------------
## Step 1: Incoming Date
## -----------------------------------------------------------------------------
system(paste0("echo LMIC Reports Google Mobility for  ",iso3c, ". Short Run = ", short_run, ". Parallel = ", parallel))
set.seed(123)
date <- as.Date(date)
short_run <- as.logical(short_run)
parallel <- as.logical(parallel)
full_scenarios <- as.logical(full_scenarios)

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
## Step 2a: Set up args for pmcmc/grid
## -----------------------------------------------------------------------------

# what is the date of first death
null_na <- function(x) {if(is.null(x)) {NA} else {x}}
min_death_date <- data$date[which(data$deaths>0)][1]

# calibration arguments
reporting_fraction = 1
R0_change <- interventions[[iso3c]]$C
date_R0_change <- interventions[[iso3c]]$date
R0_change <- R0_change[as.Date(date_R0_change) <= date]
date_R0_change <- date_R0_change[as.Date(date_R0_change) <= date]

date_contact_matrix_set_change <- NULL

squire_model <- squire::explicit_model()
pars_obs <- NULL
R0_prior <- list("func" = dnorm, args = list("mean"= 3, "sd"= 1, "log" = TRUE))
Rt_func <- function(R0_change, R0, Meff) {
  R0 * (2 * plogis(-(R0_change-1) * -Meff))
}

if(short_run) {
  n_particles <- 2
  replicates <- 2
  n_mcmc <- 100
  n_chains <- 3
  grid_spread <- 2
  sleep <- 2
} else {
  n_particles <- 50
  replicates <- 200
  n_mcmc <- 10000
  n_chains <- 3
  grid_spread <- 11
  sleep <- 120
}

if (parallel) {
  suppressWarnings(future::plan(future::multiprocess()))
}

# Defualt edges
R0_min <- 1.6
R0_max <- 5.6
Meff_min <- 0.5
Meff_max <- 10
#Meff_pl_min <- -3
Meff_pl_min <- 0
#Meff_pl_max <- 5
Meff_pl_max <- 1
last_start_date <- as.Date(null_na(min_death_date))-10
first_start_date <- as.Date(null_na(min_death_date))-55
R0_pl_shift_min <- -2
R0_pl_shift_max <- 5

## -----------------------------------------------------------------------------
## Step 2b: Sourcing previous fits to start pmcmc nearby
## -----------------------------------------------------------------------------

# 1. Do we have a previous run for this country
# 1. Do we have a previous run for this country
json <- NULL
json <- tryCatch({
  json_path <- file.path("https://raw.githubusercontent.com/mrc-ide/global-lmic-reports/master/",iso3c,"input_params_dashboard.json")
  suppressWarnings(jsonlite::read_json(json_path))
}, error = function(e){NULL})

if (!is.null(json) && !is.null(json$Meff) && json$Meff_pl[1]<=1) {
  
  R0_start <- json[[1]]$Rt
  date_start <- json[[1]]$date
  Meff_start <- json[[1]]$Meff
  Meff_pl_start <- json[[1]]$Meff_pl
  
} else {
  
  # have at least a week span for start date
  span_date_currently <- seq.Date(first_start_date, last_start_date, 1)
  day_step <- as.numeric(round((last_start_date - first_start_date + 1)/12))
  
  Sys.setenv("SQUIRE_PARALLEL_DEBUG" = "TRUE")
  
  # do coarse grid search to get in the right ball park
  out_det <- squire::calibrate(
    data = data,
    R0_min = R0_min,
    R0_max = R0_max,
    R0_step = (R0_max - R0_min)/grid_spread,
    R0_prior = R0_prior,
    Meff_min = Meff_min,
    Meff_max = Meff_max,
    Meff_step = (Meff_max - Meff_min)/grid_spread,
    Rt_func = Rt_func,
    first_start_date = first_start_date,
    last_start_date = last_start_date,
    day_step = day_step,
    squire_model = squire:::deterministic_model(),
    pars_obs = pars_obs,
    n_particles = 1,
    seeding_cases = 5,
    reporting_fraction = reporting_fraction,
    R0_change = R0_change,
    date_R0_change = date_R0_change,
    replicates = replicates,
    country = country,
    forecast = 0
  )
  
  Sys.setenv("SQUIRE_PARALLEL_DEBUG" = FALSE)
  
  ## and get the best  position
  pos <- which(out_det$scan_results$mat_log_ll == max(out_det$scan_results$mat_log_ll), arr.ind = TRUE)
  
  # get tthe R0, betas and times into a data frame
  R0_start <- out_det$scan_results$x[pos[1]]
  date_start <- as.Date(out_det$scan_results$y[pos[2]])
  Meff_start <- out_det$scan_results$z[pos[3]]
  Meff_pl_start <- 0.2
  R0_pl_shift_start <- 0
  
}

R0_start <- min(max(R0_start, R0_min), R0_max)
date_start <- min(max(as.Date(date_start), as.Date(first_start_date)), as.Date(last_start_date))
Meff_start <- min(max(Meff_start, Meff_min), Meff_max)
Meff_pl_start <- min(max(Meff_pl_start,Meff_pl_min), Meff_pl_max)
R0_pl_shift_start <- min(max(R0_pl_shift_start,R0_pl_shift_min), R0_pl_shift_max)



## -----------------------------------------------------------------------------
## Step 2b: PMCMC parameter set up
## -----------------------------------------------------------------------------

# PMCMC Parameters
pars_init = list('start_date' = date_start, 
                 'R0' = R0_start, 
                 'Meff' = Meff_start, 
                 'Meff_pl' = Meff_pl_start,
                 'R0_pl_shift' = R0_pl_shift_start,
                 )
pars_min = list('start_date' = first_start_date, 'R0' = R0_min, 'Meff' = Meff_min, 'Meff_pl' = Meff_pl_min, 'R0_pl_shift' = R0_pl_shift_min)
pars_max = list('start_date' = last_start_date, 'R0' = R0_max, 'Meff' = Meff_max, 'Meff_pl' = Meff_pl_max, 'R0_pl_shift' = R0_pl_shift_max)
pars_discrete = list('start_date' = TRUE, 'R0' = FALSE, 'Meff' = FALSE, 'Meff_pl' = FALSE, 'R0_pl_shift' = FALSE)
pars_obs = list(phi_cases = 1, k_cases = 2, phi_death = 1, k_death = 2, exp_noise = 1e6)

# Covriance Matrix
proposal_kernel <- diag(length(names(pars_init))) * 0.3
rownames(proposal_kernel) <- colnames(proposal_kernel) <- names(pars_init)
proposal_kernel["start_date", "start_date"] <- 1.5

# MCMC Functions - Prior and Likelihood Calculation
# logprior <- function(pars){
#   squire:::assert_in(names(pars), c("start_date", "R0", "Meff", "Meff_pl")) # good sanity check
#   ret <- dunif(x = pars[["start_date"]], min = -55, max = -10, log = TRUE) +
#     dnorm(x = pars[["R0"]], mean = 3, sd = 1, log = TRUE) +
#     dnorm(x = pars[["Meff"]], mean = 3, sd = 3, log = TRUE) +
#     dnorm(x = pars[["Meff_pl"]], mean = 0, sd = 3, log = TRUE)
#   return(ret)
# }

logprior <- function(pars){
  squire:::assert_in(names(pars), c("start_date", "R0", "Meff", "Meff_pl", "R0_pl_shift")) # good sanity check
  ret <- dunif(x = pars[["start_date"]], min = -55, max = -10, log = TRUE) +
    dnorm(x = pars[["R0"]], mean = 3, sd = 1, log = TRUE) +
    dnorm(x = pars[["Meff"]], mean = 3, sd = 3, log = TRUE) +
    dunif(x = pars[["Meff_pl"]], min = 0, max = 1, log = TRUE) +
    dnorm(x = pars[["R0_pl_shift"]], mean = 0, sd = 0.5, log = TRUE) 
    
  return(ret)
}

# Meff_date_change. look at when mobility has increased by 20%
above <- 1.2

# These countries have peculiar weekend effec ts that are slghtly messing with calculating this
# so have to switch the point at which we calculate their lockdown date
if (iso3c %in% c("BRA", "OMA", "USA")){
  above <- 1.1
} else if(iso3c %in% c("SWE")) {
  above <- 1.05
} else if(iso3c %in% c("MEX")) {
  above <- 1.025
}

# N.B. look at strucchange and segmented for maybe a better way to do this
pld <- post_lockdown_date(interventions[[iso3c]], 1, 
                          max_date = as.Date("2020-06-04"),
                          min_date = as.Date("2020-02-01"))
                        
# sleep so parallel is chill
Sys.sleep(time = runif(1, 0, sleep))
out_det <- squire::pmcmc(data = data, 
                         n_mcmc = n_mcmc,
                         log_prior = logprior,
                         n_particles = 1,
                         steps_per_day = 4,
                         log_likelihood = NULL,
                         squire_model = squire:::deterministic_model(),
                         output_proposals = FALSE,
                         n_chains = n_chains,
                         pars_obs = pars_obs,
                         pars_init = pars_init,
                         pars_min = pars_min,
                         pars_max = pars_max,
                         pars_discrete = pars_discrete,
                         proposal_kernel = proposal_kernel,
                         country = country, 
                         R0_change = R0_change,
                         date_R0_change = date_R0_change,
                         date_Meff_change = pld, 
                         burnin = ceiling(n_mcmc/10),
                         seeding_cases = 5,
                         replicates = replicates,
                         required_acceptance_ratio = 0.20,
                         start_adaptation = 1000,
                         roll = 7)

## -----------------------------------------------------------------------------
## Step 2b: Summarise Fits for Interface
## -----------------------------------------------------------------------------

## and save the info for the interface
all_chains <- do.call(rbind,lapply(out_det$pmcmc_results$chains, "[[", "results"))
best <- all_chains[which.max(all_chains$log_posterior), ]

# get the R0, betas and times into a data frame
R0 <- best$R0
start_date <- squire:::offset_to_start_date(data$date[1],round(best$start_date))
Meff <- best$Meff
Meff_pl <- best$Meff_pl

df_dash <- data.frame("start_date" = start_date, "R0" = R0, "Meff" = Meff, "Meff_pl" = Meff_pl, 
                      "iso3c" = iso3c, "run_date" = date, "pld" = pld)

if(!is.null(date_R0_change)) {
  tt_beta <- squire:::intervention_dates_for_odin(dates = date_R0_change,
                                                  change = R0_change,
                                                  start_date = start_date,
                                                  steps_per_day = 1)
} else {
  tt_beta <- 0
}

if(!is.null(R0_change)) {
  R0 <- squire:::evaluate_Rt_pmcmc(R0_change = tt_beta$change, R0 = R0, Meff = Meff, 
                             Meff_pl = Meff_pl, start_date = start_date,
                             date_R0_change = date_R0_change[date_R0_change>start_date], 
                             date_Meff_change = out_det$pmcmc_results$inputs$interventions$date_Meff_change,
                             roll = 7)
} else {
  R0 <- R0
}
beta_set <- squire:::beta_est(squire_model = squire_model,
                              model_params = out_det$pmcmc_results$inputs$model_params,
                              R0 = R0)

df <- data.frame(tt_beta = c(0,tt_beta$tt), beta_set = beta_set, 
                 date = start_date + c(0,tt_beta$tt), Rt = R0, 
                 grey_bar_start = FALSE)

# add in grey bar start for interface
ox_interventions <- readRDS("oxford_grt.rds")
ox_interventions_unique <- squire:::interventions_unique(ox_interventions[[iso3c]], "C")
df$grey_bar_start[which.min(abs(as.numeric(df$date - ox_interventions_unique$dates_change[1])))] <- TRUE

writeLines(jsonlite::toJSON(df,pretty = TRUE), "input_params.json")
writeLines(jsonlite::toJSON(df_dash,pretty = TRUE), "input_params_dashboard.json")


## -----------------------------------------------------------------------------
## Step 3: Particle Filter
## -----------------------------------------------------------------------------

## -----------------------------------------------------------------------------
## Step 3a: Fit Stochastic Model
## -----------------------------------------------------------------------------

Sys.setenv("SQUIRE_PARALLEL_DEBUG" = "TRUE")

out <- out_det
out$pmcmc_results$inputs$squire_model <- explicit_model()
out$pmcmc_results$inputs$model_params$dt <- 0.05
pmcmc <- out$pmcmc_results
out <- generate_draws_pmcmc(pmcmc = pmcmc,
                            burnin = ceiling(n_mcmc/10),
                            n_chains = n_chains,
                            squire_model = out$pmcmc_results$inputs$squire_model,
                            replicates = replicates,
                            n_particles = n_particles,
                            forecast = 0,
                            country = country,
                            population = squire::get_population(iso3c = iso3c)$n,
                            interventions = out$interventions,
                            data = out$pmcmc_results$inputs$data)

# Add the prior
out$pmcmc_results$inputs$prior <- as.function(c(formals(logprior), 
                                                  body(logprior)), 
                                                envir = new.env(parent = environment(stats::acf)))

## -----------------------------------------------------------------------------
## Step 3d: Summarise Fits
## -----------------------------------------------------------------------------

## summarise what we have
png("top_row.png", height = 6, width = 8, units = "in", res = 300)
plot(out$pmcmc_results, thin = 0.25)
dev.off()

title <- cowplot::ggdraw() + 
  cowplot::draw_label(
    country,
    fontface = 'bold',
    x = 0.5
  )

line <- ggplot() + cowplot::draw_line(x = 0:10, y=1) + 
  theme(panel.background = element_blank(),
        axis.title = element_blank(), 
        axis.text = element_blank(), 
        axis.ticks = element_blank())

header <- cowplot::plot_grid(title, line, ncol = 1)

png("header.png", height = 0.5, width = 8, units = "in", res = 300)
header
dev.off()


index <- squire:::odin_index(out$model)
forecast <- 0

suppressWarnings(d <- deaths_plot_single(out, data, date = date,date_0 = date_0, forecast = forecast, single = TRUE) + 
                   theme(legend.position = "none"))

intervention <- intervention_plot_google(interventions[[iso3c]], date, data, forecast)


bottom <- cowplot::plot_grid(intervention, d,
                             ncol=1,
                             rel_heights = c(0.4,0.6))
cowplot::save_plot("bottom.png", bottom, base_height = 6, base_width = 8)

plots <- list() 

nms <- c("header.png", "top_row.png", "bottom.png")
for(i in 1:3) {
  x <- nms[i]
  img <- png::readPNG(x)
  plots[[i]] <- grid::rasterGrob(img, interpolate = FALSE)
}

ggsave("fitting.pdf",width=7, height=11, 
       gridExtra::marrangeGrob(grobs = c(list(cowplot::as_grob(header)),
                                         plots[2],
                                         list(cowplot::as_grob(bottom))), 
                               nrow=3, ncol=1,top=NULL, heights = c(1.5, 6, 6)))
dev.off()
file.remove(nms)
file.remove("Rplots.pdf")

## Save the grid out object

# remove states to keep object memory save down
for(i in seq_along(out$pmcmc_results$chains)) {
  out$pmcmc_results$chains[[i]]$states <- NULL
}

saveRDS(out, "grid_out.rds")

## -----------------------------------------------------------------------------
## Step 4: Scenarios
## -----------------------------------------------------------------------------

## -----------------------------------------------------------------------------
## 4.1. Conduct our projections based on no surging
## -----------------------------------------------------------------------------

## Functions for working out the relative changes in R0 for given scenarios
fr0 <- tail(out$interventions$R0_change,1)
time_period <- 365
rel_R0 <- function(rel = 0.5, Meff_mult = 1) {
  R0_ch <- 1-((1-fr0)*rel)
  
  current <-  vapply(seq_along(out$replicate_parameters$R0), function(y){
    
    if(!is.null(date_R0_change)) {
      tt_beta <- squire:::intervention_dates_for_odin(dates = out$interventions$date_R0_change,
                                                      change = out$interventions$R0_change,
                                                      start_date = out$replicate_parameters$start_date[y],
                                                      steps_per_day = 1/out$parameters$dt)
    } else {
      tt_beta <- 0
    }
    
    tail(squire:::evaluate_Rt_pmcmc(
      R0_change = c(out$interventions$R0_change[out$interventions$date_R0_change>out$replicate_parameters$start_date[y]],fr0), 
      R0 = out$replicate_parameters$R0[y], 
      Meff = out$replicate_parameters$Meff[y], 
      Meff_pl = out$replicate_parameters$Meff_pl[y],
      date_R0_change = c(out$interventions$date_R0_change[out$interventions$date_R0_change>out$replicate_parameters$start_date[y]],as.Date(date)+1),
      date_Meff_change = out$interventions$date_Meff_change, 
      roll = out$pmcmc_results$inputs$roll,
      start_date = out$replicate_parameters$start_date[y]),1)
  }, numeric(1))
  
  wanted <-  vapply(seq_along(out$replicate_parameters$R0), function(y){
    
    if(!is.null(date_R0_change)) {
      tt_beta <- squire:::intervention_dates_for_odin(dates = out$interventions$date_R0_change,
                                                      change = out$interventions$R0_change,
                                                      start_date = out$replicate_parameters$start_date[y],
                                                      steps_per_day = 1/out$parameters$dt)
    } else {
      tt_beta <- 0
    }
    
    tail(squire:::evaluate_Rt_pmcmc(
      R0_change = c(out$interventions$R0_change[out$interventions$date_R0_change>out$replicate_parameters$start_date[y]],R0_ch), 
      R0 = out$replicate_parameters$R0[y], 
      Meff = out$replicate_parameters$Meff[y], 
      Meff_pl = out$replicate_parameters$Meff_pl[y],
      date_R0_change = c(out$interventions$date_R0_change[out$interventions$date_R0_change>out$replicate_parameters$start_date[y]],as.Date(date)+1),
      date_Meff_change = out$interventions$date_Meff_change, 
      roll = out$pmcmc_results$inputs$roll,
      start_date = out$replicate_parameters$start_date[y]),1)
  }, numeric(1))
  
  mean(wanted/current)
}

# Maintaining the current set of measures for a further 3 months following which contacts return to pre-intervention levels  
maintain_3months_lift <- squire::projections(out, 
                                             R0_change = c(1, rel_R0(0)), 
                                             tt_R0 = c(0,90), 
                                             time_period = time_period)

# Enhancing movement restrictions for 3 months (50% further reduction in contacts) which then return to pre-intervention levels (Mitigation) 
mitigation_3months_lift <- squire::projections(out, 
                                               R0_change = c(0.5, rel_R0(0)), 
                                               tt_R0 = c(0,90), 
                                               time_period = time_period)

# Relax by 50% for 3 months and then return to pre-intervention levels 
reverse_3_months_lift <- squire::projections(out, 
                                             R0_change = c(rel_R0(0.5), rel_R0(0)), 
                                             tt_R0 = c(0, 90), 
                                             time_period = time_period)

## -----------------------------------------------------------------------------
## 4.2. Investigating a capacity surge
## -----------------------------------------------------------------------------

## Is capacity passed prior to today
icu_cap <- squire:::get_ICU_bed_capacity(country)
hosp_cap <- squire:::get_hosp_bed_capacity(country)

t_test_safe <- function(x, ...) {
  out <- try(t.test(x, ...), silent = TRUE)
  if (inherits(out, "try-error"))
  {
    out <- list("conf.int"=c(mean(x),mean(x)))
  }
  return(out)
}

# update for surged healthcare
pmcmc <- out$pmcmc_results
pmcmc$inputs$model_params$hosp_beds <- 1e10
pmcmc$inputs$model_params$ICU_beds <- 1e10

out_surged <- generate_draws_pmcmc(pmcmc = pmcmc,
                                   burnin = ceiling(n_mcmc/10),
                                   n_chains = n_chains,
                                   squire_model = out$pmcmc_results$inputs$squire_model, 
                                   replicates = replicates, 
                                   n_particles = n_particles, 
                                   forecast = 0, 
                                   country = country, 
                                   population = squire::get_population(iso3c = iso3c)$n, 
                                   interventions = out$interventions, 
                                   data = out$pmcmc_results$inputs$data)

out_surged$model$set_user(ICU_beds = 1e10)
out_surged$model$set_user(hosp_beds = 1e10)

out_surged$parameters$hosp_bed_capacity <- 1e10
out_surged$parameters$ICU_bed_capacity <- 1e10

# will surging be required within the next 28 day forecast
icu <- format_output(maintain_3months_lift, "ICU_demand", date_0 = date_0)
icu <- icu[icu$compartment == "ICU_demand",]
icu_28 <- group_by(icu[icu$t==28,], replicate) %>% 
  summarise(tot = sum(y, na.rm = TRUE)) %>% 
  summarise(i_tot = mean(tot, na.rm = TRUE), 
            i_min = t_test_safe(tot)$conf.int[1],
            i_max = t_test_safe(tot)$conf.int[2])

hosp <- format_output(maintain_3months_lift, "hospital_demand", date_0 = date_0)
hosp <- hosp[hosp$compartment == "hospital_demand",]
hosp_28 <- group_by(hosp[hosp$t==28,], replicate) %>% 
  summarise(tot = sum(y, na.rm = TRUE)) %>% 
  summarise(i_tot = mean(tot, na.rm = TRUE), 
            i_min = t_test_safe(tot)$conf.int[1],
            i_max = t_test_safe(tot)$conf.int[2])

if(icu_28$i_tot > icu_cap || hosp_28$i_tot > hosp_cap) {
  
  surging <- TRUE
  
} else {
  
  surging <- FALSE
  
}


## Now simulate the surged scenarios
out_surged$pmcmc_results$chains <- NULL

maintain_3months_lift_surged <- squire::projections(out_surged, 
                                                    R0_change = c(1, rel_R0(0)), 
                                                    tt_R0 = c(0,90), 
                                                    time_period = time_period)

mitigation_3months_lift_surged <- squire::projections(out_surged, 
                                                      R0_change = c(0.5, rel_R0(0)), 
                                                      tt_R0 = c(0,90), 
                                                      time_period = time_period)

reverse_3_months_lift_surged <- squire::projections(out_surged, 
                                                    R0_change = c(rel_R0(0.5), rel_R0(0)), 
                                                    tt_R0 = c(0, 90), 
                                                    time_period = time_period)

## -----------------------------------------------------------------------------
## 4.3. Extra Scenarios
## N.B. Not ready yet  ---------------------------------------------------------
## Need to think about these more and think about putting in linear mobility time
## -----------------------------------------------------------------------------

if (full_scenarios) {
  
  # # Enhancing movement restrictions until the end of the year (50% further reduction in contacts) which then return to pre-intervention levels   
  # mitigation_rest_year_lift <- squire::projections(out, R0_change = c(0.5, rel_R0(0)), 
  #                                                  tt_R0 = c(0,as.Date("2020-12-31")-Sys.Date()), time_period = time_period)
  # 
  # # Enhance movement restrictions for 3 months (50% further reduction in contacts), ease restrictions for a further 3 months (current contacts), then return to pre-intervention levels  
  # mitigation_3months_ease_3months_lift <- squire::projections(out, R0_change = c(0.5, 1, rel_R0(0)), tt_R0 = c(0, 90, 180), time_period = time_period)
  # 
  # # Suppression for 3 months (75% reduction in contacts) which then return to pre-intervention levels 
  # suppress_3months_lift <- squire::projections(out, R0_change = c((rel_R0(0))*0.25, rel_R0(0)), tt_R0 = c(0, 90), time_period = time_period)
  # 
  # # Long-term sustained suppression (75% reduction in contacts)  
  # suppress_full <- squire::projections(out, R0_change = c((rel_R0(0))*0.25), tt_R0 = c(0), time_period = time_period)
  # 
  # # Full lifting of emergency measures in a week – contact rates are assumed to return to pre-intervention levels  
  # lift_week <- squire::projections(out, R0_change = c(1,rel_R0(0)), tt_R0 = c(0, 7), time_period = time_period)
  # 
  # ## -----------------------------------------------------------------------------
  # ## 4.2. Assuming that this will uncoupel further, with Meff being 50% lower
  # ## -----------------------------------------------------------------------------
  # 
  # # Maintaining the current set of measures for a further 3 months following which contacts return to pre-intervention levels  
  # maintain_3months_lift_meff_50 <- squire::projections(out, 
  #                                                      R0_change = c(1, rel_R0(0)*0.5), 
  #                                                      tt_R0 = c(0,90), 
  #                                                      time_period = time_period)
  # 
  # # Enhancing movement restrictions for 3 months (50% further reduction in contacts) which then return to pre-intervention levels (Mitigation) 
  # mitigation_3months_lift_meff_50 <- squire::projections(out, R0_change = c(0.5,rel_R0(0)*0.5), tt_R0 = c(0,90), time_period = time_period)
  # 
  # # Relax by 50% for 3 months and then return to pre-intervention levels 
  # reverse_3_months_lift_meff_50 <- squire::projections(out, R0_change = c((rel_R0(0))*(1-(fr0/2)), rel_R0(0)*0.5), tt_R0 = c(0, 90), time_period = time_period)
  # 
  # # Enhancing movement restrictions until the end of the year (50% further reduction in contacts) which then return to pre-intervention levels   
  # mitigation_rest_year_lift_meff_50 <- squire::projections(out, R0_change = c(0.5, rel_R0(0)*0.5), tt_R0 = c(0,as.Date("2020-12-31")-Sys.Date()), time_period = time_period)
  # 
  # # Enhance movement restrictions for 3 months (50% further reduction in contacts), ease restrictions for a further 3 months (current contacts), then return to pre-intervention levels  
  # mitigation_3months_ease_3months_lift_meff_50 <- squire::projections(out, R0_change = c(0.5, 1*0.5, rel_R0(0)*0.5), tt_R0 = c(0, 90, 180), time_period = time_period)
  # 
  # # Suppression for 3 months (75% reduction in contacts) which then return to pre-intervention levels 
  # suppress_3months_lift_meff_50 <- squire::projections(out, R0_change = c((rel_R0(0))*0.25, rel_R0(0)*0.5), tt_R0 = c(0, 90), time_period = time_period)
  # 
  # # Long-term sustained suppression (75% reduction in contacts)  
  # suppress_full_meff_50 <- squire::projections(out, R0_change = c((rel_R0(0))*0.25), tt_R0 = c(0), time_period = time_period)
  # 
  # # Full lifting of emergency measures in a week – contact rates are assumed to return to pre-intervention levels  
  # lift_week_meff_50 <- squire::projections(out, R0_change = c(1,rel_R0(0)*0.5), tt_R0 = c(0, 7), time_period = time_period)
  # 
  # ## -----------------------------------------------------------------------------
  # ## 4.3. Assuming that this will uncoupel further, with Meff being 25% lower
  # ## -----------------------------------------------------------------------------
  # 
  # # Maintaining the current set of measures for a further 3 months following which contacts return to pre-intervention levels  
  # maintain_3months_lift_meff_75 <- squire::projections(out, R0_change = c(1, rel_R0(0)*0.75), tt_R0 = c(0,90), time_period = time_period)
  # 
  # # Enhancing movement restrictions for 3 months (50% further reduction in contacts) which then return to pre-intervention levels (Mitigation) 
  # mitigation_3months_lift_meff_75 <- squire::projections(out, R0_change = c(0.5,rel_R0(0)*0.75), tt_R0 = c(0,90), time_period = time_period)
  # 
  # # Relax by 50% for 3 months and then return to pre-intervention levels 
  # reverse_3_months_lift_meff_75 <- squire::projections(out, R0_change = c((rel_R0(0))*(1-(fr0/2)), rel_R0(0)*0.75), tt_R0 = c(0, 90), time_period = time_period)
  # 
  # # Enhancing movement restrictions until the end of the year (50% further reduction in contacts) which then return to pre-intervention levels   
  # mitigation_rest_year_lift_meff_75 <- squire::projections(out, R0_change = c(0.5, rel_R0(0)*0.75), tt_R0 = c(0,as.Date("2020-12-31")-Sys.Date()), time_period = time_period)
  # 
  # # Enhance movement restrictions for 3 months (50% further reduction in contacts), ease restrictions for a further 3 months (current contacts), then return to pre-intervention levels  
  # mitigation_3months_ease_3months_lift_meff_75 <- squire::projections(out, R0_change = c(0.5, 1*0.75, rel_R0(0)*0.75), tt_R0 = c(0, 90, 180), time_period = time_period)
  # 
  # # Suppression for 3 months (75% reduction in contacts) which then return to pre-intervention levels 
  # suppress_3months_lift_meff_75 <- squire::projections(out, R0_change = c((rel_R0(0))*0.25, rel_R0(0)*0.75), tt_R0 = c(0, 90), time_period = time_period)
  # 
  # # Long-term sustained suppression (75% reduction in contacts)  
  # suppress_full_meff_75 <- squire::projections(out, R0_change = c((rel_R0(0))*0.25), tt_R0 = c(0), time_period = time_period)
  # 
  # # Full lifting of emergency measures in a week – contact rates are assumed to return to pre-intervention levels  
  # lift_week_meff_75 <- squire::projections(out, R0_change = c(1,rel_R0(0)*0.75), tt_R0 = c(0, 7), time_period = time_period)
  # 
  # 
  # ## BIND THESE TOGETHER
  # 
  # r_list <- named_list(maintain_3months_lift, mitigation_3months_lift,  reverse_3_months_lift, 
  #                      mitigation_rest_year_lift, mitigation_3months_ease_3months_lift, 
  #                      suppress_3months_lift, suppress_full, lift_week,
  #                      
  #                      maintain_3months_lift_meff_50, mitigation_3months_lift_meff_50, reverse_3_months_lift_meff_50, 
  #                      mitigation_3months_lift_meff_50, mitigation_3months_ease_3months_lift_meff_50, 
  #                      suppress_3months_lift_meff_50, suppress_full_meff_50, lift_week_meff_50,
  #                      
  #                      maintain_3months_lift_meff_75, mitigation_3months_lift_meff_75, reverse_3_months_lift_meff_75, 
  #                      mitigation_rest_year_lift_meff_75, mitigation_3months_ease_3months_lift_meff_75, 
  #                      suppress_3months_lift_meff_75, suppress_full_meff_75, lift_week_meff_75)
  
} else {
  
  r_list <-
    named_list(
      maintain_3months_lift,
      mitigation_3months_lift,
      reverse_3_months_lift,
      maintain_3months_lift_surged,
      mitigation_3months_lift_surged,
      reverse_3_months_lift_surged
    )
  
}

r_list_pass <- r_list

o_list <- lapply(r_list_pass, squire::format_output,
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
                  params = list("r_list" = r_list_pass,
                                "o_list" = o_list,
                                "replicates" = replicates, 
                                "data" = data,
                                "date_0" = date_0,
                                "country" = country,
                                "surging" = surging),
                  output_options = list(pandoc_args = c(paste0("--metadata=title:",country," COVID-19 report "))))

# summarise the projections

data_sum <- lapply(o_list, function(pd){
  
  # remove any NA rows (due to different start dates)
  if(sum(is.na(pd$t) | is.na(pd$y))>0) {
    pd <- pd[-which(is.na(pd$t) | is.na(pd$y)),]
  }
  
  # add in cumulative infections
  cum_i <- pd %>% 
    filter(compartment == "infections") %>% 
    group_by(replicate) %>% 
    mutate(y = cumsum(y),
           compartment = "cumulative_infections") %>% 
    ungroup
  
  pd <- rbind(pd, cum_i)
  
  # Format summary data
  pds <- pd %>%
    dplyr::filter(.data$date < (date_0+90)) %>% 
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
data_sum[[4]]$scenario <- "Surged Maintain Status Quo"
data_sum[[5]]$scenario <- "Surged Additional 50% Reduction"
data_sum[[6]]$scenario <- "Surged Relax Interventions 50%"

# summarise the Rt
rt_sum <- lapply(r_list_pass, rt_creation, date_0, date_0+89)
rt_sum[[1]]$scenario <- "Maintain Status Quo"
rt_sum[[2]]$scenario <- "Additional 50% Reduction"
rt_sum[[3]]$scenario <- "Relax Interventions 50%"
rt_sum[[4]]$scenario <- "Surged Maintain Status Quo"
rt_sum[[5]]$scenario <- "Surged Additional 50% Reduction"
rt_sum[[6]]$scenario <- "Surged Relax Interventions 50%"

# combine and annotate
data_sum <- do.call(rbind, data_sum)
rt_sum <- do.call(rbind, rt_sum)
data_sum <- rbind(data_sum, rt_sum) %>% arrange(date, scenario)
rownames(data_sum) <- NULL

data_sum$country <- country
data_sum$iso3c <- iso3c
data_sum$report_date <- date
data_sum <- data_sum[data_sum$compartment != "D",]
data_sum$version <- "v3"
write.csv(data_sum, "projections.csv", row.names = FALSE, quote = FALSE)

## -----------------------------------------------------------------------------
## Step 6: Full saves
## -----------------------------------------------------------------------------

if (full_scenarios) {
  
  o_list <- lapply(r_list, squire::format_output,
                   var_select = c("infections","deaths","hospital_demand","ICU_demand", "ICase"),
                   date_0 = date_0)
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
  
  for(i in seq_along(r_list)) {
    data_sum[[i]]$scenario <- names(r_list)[i]
  }
  data_sum <- do.call(rbind, data_sum)
  data_sum$country <- country
  data_sum$iso3c <- iso3c
  data_sum$report_date <- date
  write.csv(data_sum, "full_projections.csv", row.names = FALSE, quote = FALSE)
  
}
