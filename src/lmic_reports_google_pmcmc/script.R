orderly_id <- tryCatch(orderly::orderly_run_info()$id,
                       error = function(e) "<id>") # bury this in the html, docx

version_min <- "0.4.18"
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
## Step 2a: Sourcing previous fits to start pmcmc nearby
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
R0_prior <- list("func" = dnorm, args = list("mean"= 3, "sd"= 0.5, "log" = TRUE))
Rt_func <- function(R0_change, R0, Meff) {
  R0 * (2 * plogis(-(R0_change-1) * -Meff))
}

if(short_run) {
  n_particles <- 2
  replicates <- 2
  n_mcmc <- 5
  n_chains <- 3
} else {
  n_particles <- 5
  replicates <- 100
  n_mcmc <- 7500
  n_chains < 3
}

# 1. Do we have a previous run for this country
json <- NULL
json <- tryCatch({
  json_path <- file.path("https://raw.githubusercontent.com/mrc-ide/global-lmic-reports/master/",iso3c,"input_params.json")
  suppressWarnings(jsonlite::read_json(json_path))
}, error = function(e){NULL})

if (!is.null(json) && !is.null(json$Meff_pl)) {
  
  R0 <- json[[1]]$Rt
  date <- json[[1]]$date
  Meff <- json[[1]]$Meff
  
  # get the range from this for R0 and grow it by 0.75
  R0_max <- min(R0 + 0.75, 5.6)
  R0_min <- max(R0 - 0.75, 1.0)
  
  # get the range for dates and grow it by 7 days
  last_start_date <- as.Date(date) + 7
  first_start_date <- as.Date(date) - 7
  
  # adust the dates so they are compliant with the data
  last_start_date <- min(c(last_start_date, as.Date(null_na(min_death_date))-10), na.rm = TRUE)
  first_start_date <- max(as.Date("2020-01-04"), first_start_date, na.rm = TRUE)
  
  # get the range for Meff
  Meff_max <- min(Meff + 1, 7.5)
  Meff_min <- max(Meff - 1, 0.1)
  
  # get the range for Meff_pl
  Meff_pl_max <- min(Meff + 1, 15)
  Meff_pl_min <- max(Meff - 1, 0.1)
  
} else {
  
  # Defualts if no previous data
  R0_min <- 1.6
  R0_max <- 5.6
  Meff_min <- 0.1
  Meff_max <- 7.5
  Meff_pl_min <- 0.1
  Meff_pl_max <- 15
  last_start_date <- as.Date(null_na(min_death_date))-10
  first_start_date <- max(as.Date("2020-01-04"),last_start_date - 45, na.rm = TRUE)
  
}

if (parallel) {
  suppressWarnings(future::plan(future::multiprocess()))
}

## -----------------------------------------------------------------------------
## Step 2b: PMCMC parameter set up
## -----------------------------------------------------------------------------

# PMCMC Parameters
pars_init = list('start_date' = first_start_date + (last_start_date - first_start_date)/2, 
                 'R0' = R0_min + (R0_max - R0_min)/2, 
                 'Meff' = Meff_min + (Meff_max - Meff_min)/2, 
                 'Meff_pl' = (Meff_min + (Meff_max - Meff_min)/2)*1.5)
pars_min = list('start_date' = first_start_date, 'R0' = R0_min, 'Meff' = Meff_min, 'Meff_pl' = Meff_pl_min)
pars_max = list('start_date' = last_start_date, 'R0' = R0_max, 'Meff' = Meff_max, 'Meff_pl' = Meff_pl_max)
pars_discrete = list('start_date' = TRUE, 'R0' = FALSE, 'Meff' = FALSE, 'Meff_pl' = FALSE)
pars_obs = list(phi_cases = 1, k_cases = 2, phi_death = 1, k_death = 2, exp_noise = 1e6)

# Covriance Matrix
proposal_kernel <- diag(length(names(pars_init))) * 0.2
rownames(proposal_kernel) <- colnames(proposal_kernel) <- names(pars_init)
proposal_kernel["start_date", "start_date"] <- 1.5

# MCMC Functions - Prior and Likelihood Calculation
logprior <- function(pars){
  squire:::assert_in(names(pars), c("start_date", "R0", "Meff", "Meff_pl")) # good sanity check
  ret <- dunif(x = pars[["start_date"]], min = -55, max = -10, log = TRUE) +
    dnorm(x = pars[["R0"]], mean = 3, sd = 0.5, log = TRUE) +
    dnorm(x = pars[["Meff"]], mean = 2, sd = 2, log = TRUE) +
    dnorm(x = pars[["Meff_pl"]], mean = 6, sd = 3, log = TRUE)
  return(ret)
}

# Meff_date_change
pld <- post_lockdown_date(interventions[[iso3c]], 1)

out_det <- squire::pmcmc(data = data, 
           n_mcmc = n_mcmc,
           log_prior = logprior,
           n_particles = 1,
           steps_per_day = 4,
           log_likelihood = NULL,
           squire_model = squire:::deterministic_model(),
           output_proposals = FALSE,
           n_chains = n_chains,
           Rt_func = Rt_func,
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
           required_acceptance_ratio = 0.13,
           start_covariance_adaptation = 150,
           start_scaling_factor_adaptation = 10,
           initial_scaling_factor = 0.05
)

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

if(!is.null(date_R0_change)) {
  tt_beta <- squire:::intervention_dates_for_odin(dates = date_R0_change,
                                                  change = R0_change,
                                                  start_date = start_date,
                                                  steps_per_day = 1)
} else {
  tt_beta <- 0
}

if(!is.null(R0_change)) {
  R0 <- squire:::evaluate_Rt(R0_change = tt_beta$change, R0 = R0, Meff = Meff, 
                             Meff_pl = Meff_pl, 
                             date_R0_change = date_R0_change[date_R0_change>start_date], 
                             date_Meff_change = out_det$pmcmc_results$inputs$interventions$date_Meff_change,
                             Rt_func = Rt_func)
} else {
  R0 <- R0
}
beta_set <- squire:::beta_est(squire_model = squire_model,
                              model_params = out_det$pmcmc_results$inputs$model_params,
                              R0 = R0)

df <- data.frame(tt_beta = c(0,tt_beta$tt), beta_set = beta_set, 
                 date = start_date + c(0,tt_beta$tt), Rt = R0, Meff = Meff,
                 grey_bar_start = FALSE)

# add in grey bar start for interface
ox_interventions <- readRDS("oxford_grt.rds")
ox_interventions_unique <- squire:::interventions_unique(ox_interventions[[iso3c]], "C")
df$grey_bar_start[which.min(abs(as.numeric(df$date - ox_interventions_unique$dates_change[1])))] <- TRUE

writeLines(jsonlite::toJSON(df,pretty = TRUE), "input_params.json")


## -----------------------------------------------------------------------------
## Step 3: Particle Filter
## -----------------------------------------------------------------------------

## -----------------------------------------------------------------------------
## Step 3a: Fit Detreministic Model Again
## -----------------------------------------------------------------------------

# First redo the deterministic model fit based on 20 seeds to get the parameter space
out_det <- squire::pmcmc(data = data, 
                         n_mcmc = n_mcmc*2,
                         log_prior = logprior,
                         n_particles = 1,
                         steps_per_day = 4,
                         log_likelihood = NULL,
                         squire_model = squire:::deterministic_model(),
                         output_proposals = FALSE,
                         n_chains = n_chains,
                         Rt_func = Rt_func,
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
                         replicates = replicates,
                         required_acceptance_ratio = 0.13,
                         start_covariance_adaptation = 150,
                         start_scaling_factor_adaptation = 10,
                         initial_scaling_factor = 0.05
)

## take the density from the deterministic to focus the grid
# recreate the grids

## get the density
all_chains <- unique(do.call(rbind,lapply(out_det$pmcmc_results$chains, "[[", "results")))

# first get the sorted density
all_chains <- all_chains[order(all_chains$log_posterior, decreasing = TRUE),]
all_chains$prob <- exp(all_chains$log_posterior)/sum(exp(all_chains$log_posterior))
drop <- 0.9
while(any(is.na(all_chains$prob))) {
  all_chains$prob <- exp(all_chains$log_posterior*drop)/sum(exp(all_chains$log_posterior*drop))
  drop <- drop^2
}
cum <- cumsum(all_chains$prob)
cut <- which(cum > 0.8)[1]

# get the range from this for R0 and grow it by 0.2
R0_max <- min(max(all_chains$R0[seq_len(cut)]) + 0.4, 5.6)
R0_min <- max(min(all_chains$R0[seq_len(cut)]) - 0.4, 1.0)

# get the range for dates and grow it by 1 day
last_start_date <- squire:::offset_to_start_date(data$date[1], round(max(all_chains$start_date))) + 4
first_start_date <- squire:::offset_to_start_date(data$date[1], round(min(all_chains$start_date))) - 4

# adust the dates so they are compliant with the data
last_start_date <- min(c(last_start_date, as.Date(null_na(min_death_date))-10), na.rm = TRUE)
first_start_date <- max(as.Date("2020-01-04"), first_start_date, na.rm = TRUE)

# get the range for Meff
Meff_max <- min(max(all_chains$Meff[seq_len(cut)]) + 0.25, 10)
Meff_min <- max(min(all_chains$Meff[seq_len(cut)]) - 0.25, 0.1)

# get the range for Meff
Meff_pl_max <- min(max(all_chains$Meff_pl[seq_len(cut)]) + 0.25, 20)
Meff_pl_min <- max(min(all_chains$Meff_pl[seq_len(cut)]) - 0.25, 0.1)

# PMCMC Parameters
pars_init = list(
  list('start_date' = first_start_date, 'R0' = R0_min, 'Meff' = Meff_min, 'Meff_pl' = Meff_pl_min),
  list('start_date' = last_start_date, 'R0' = R0_max, 'Meff' = Meff_max, 'Meff_pl' = Meff_pl_max),
  list('start_date' = first_start_date + (last_start_date - first_start_date)/2, 
       'R0' = R0_min + (R0_max - R0_min)/2, 
       'Meff' = Meff_min + (Meff_max - Meff_min)/2, 
       'Meff_pl' = (Meff_pl_min + (Meff_pl_max - Meff_pl_min)/2))
)
pars_min = list('start_date' = first_start_date, 'R0' = R0_min, 'Meff' = Meff_min, 'Meff_pl' = Meff_pl_min)
pars_max = list('start_date' = last_start_date, 'R0' = R0_max, 'Meff' = Meff_max, 'Meff_pl' = Meff_pl_max)
pars_discrete = list('start_date' = TRUE, 'R0' = FALSE, 'Meff' = FALSE, 'Meff_pl' = FALSE)
pars_obs = list(phi_cases = 1, k_cases = 2, phi_death = 1, k_death = 2, exp_noise = 1e6)

## -----------------------------------------------------------------------------
## Step 3b: Fit Stochastic Model
## -----------------------------------------------------------------------------

out <- squire::pmcmc(data = data, 
                         n_mcmc = n_mcmc,
                         log_prior = logprior,
                         n_particles = n_particles,
                         steps_per_day = 4,
                         log_likelihood = NULL,
                         squire_model = squire:::explicit_model(),
                         output_proposals = FALSE,
                         n_chains = n_chains,
                         Rt_func = Rt_func,
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
                         replicates = replicates,
                         required_acceptance_ratio = 0.13,
                         start_covariance_adaptation = 150,
                         start_scaling_factor_adaptation = 10,
                         initial_scaling_factor = 0.05
)

# Reassign the Rt_func_replace with stats environment as for some reason this is grabbing the environment
# and adding like 50Mb plus to the object being saved!
Rt_func_replace <- function(R0_change, R0, Meff) {
  R0 * (2 * plogis(-(R0_change-1) * -Meff))
}
out$pmcmc_results$inputs$Rt_func <- as.function(c(formals(Rt_func_replace), 
                                                 body(Rt_func_replace)), 
                                               envir = new.env(parent = environment(stats::acf)))
saveRDS(out, "grid_out.rds")

## -----------------------------------------------------------------------------
## Step 3b: Summarise Fits
## -----------------------------------------------------------------------------

## summarise what we have
top_row <- plot(out$pmcmc_results)
top_row <- recordPlot()

index <- squire:::odin_index(out$model)
forecast <- 0

d <- deaths_plot_single(out, data, date = date,date_0 = date_0, forecast = forecast, single = TRUE) + 
  theme(legend.position = "none")

#intervention <- intervention_plot(interventions[[iso3c]], date)
intervention <- intervention_plot_google(interventions[[iso3c]], date, data, forecast)

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

pdf("fitting.pdf",width = 8.5,height = 12)
suppressWarnings(print(cowplot::plot_grid(title,line,top_row,intervention,d,ncol=1,rel_heights = c(0.1,0.1,1,0.4,0.6))))
dev.off()

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
  wanted <- vapply(seq_along(out$replicate_parameters$R0), function(i){
    out$pmcmc_results$inputs$Rt_func(R0_ch, 
                                    R0 = out$replicate_parameters$R0[i], 
                                    Meff = out$replicate_parameters$Meff[i]*Meff_mult)
  }, numeric(1))
  current <- vapply(seq_along(out$replicate_parameters$R0), function(i){
    out$pmcmc_results$inputs$Rt_func(fr0, 
                                    R0 = out$replicate_parameters$R0[i], 
                                    Meff = out$replicate_parameters$Meff[i])
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

icu <- format_output(maintain_3months_lift, "ICU_demand", date_0 = date_0)
icu <- icu[icu$compartment == "ICU_demand",]
icu_0 <- group_by(icu[icu$t==0,], replicate) %>% 
  summarise(tot = sum(y, na.rm = TRUE)) %>% 
  summarise(i_tot = mean(tot, na.rm = TRUE), 
            i_min = t_test_safe(tot)$conf.int[1],
            i_max = t_test_safe(tot)$conf.int[2])

hosp <- format_output(maintain_3months_lift, "hospital_demand", date_0 = date_0)
hosp <- hosp[hosp$compartment == "hospital_demand",]
hosp_0 <- group_by(hosp[hosp$t==0,], replicate) %>% 
  summarise(tot = sum(y, na.rm = TRUE)) %>% 
  summarise(i_tot = mean(tot, na.rm = TRUE), 
            i_min = t_test_safe(tot)$conf.int[1],
            i_max = t_test_safe(tot)$conf.int[2])

# if it is then we need to redo the fit without it for our surged predictions
if(icu_0$i_tot > icu_cap || hosp_0$i_tot > hosp_cap) {
  
  out_surged <- squire::pmcmc(data = data, 
                       n_mcmc = n_mcmc,
                       log_prior = logprior,
                       n_particles = n_particles,
                       steps_per_day = 4,
                       log_likelihood = NULL,
                       squire_model = out$pmcmc_results$inputs$squire_model,
                       output_proposals = FALSE,
                       n_chains = n_chains,
                       Rt_func = Rt_func,
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
                       replicates = replicates,
                       required_acceptance_ratio = 0.13,
                       start_covariance_adaptation = 150,
                       start_scaling_factor_adaptation = 10,
                       initial_scaling_factor = 0.05,
                       baseline_hosp_bed_capacity = 1e10, 
                       baseline_ICU_bed_capacity = 1e10
  )
 
  surging <- TRUE
  
} else {
  
  # if not then we can use the current fit with altered capacity:
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
  
}

## Now simulate the surged scenarios

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
