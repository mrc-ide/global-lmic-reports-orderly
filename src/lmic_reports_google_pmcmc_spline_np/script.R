orderly_id <- tryCatch(orderly::orderly_run_info()$id,
                       error = function(e) "<id>") 

print(sessionInfo())
RhpcBLASctl::blas_set_num_threads(1L)
RhpcBLASctl::omp_set_num_threads(1L)

version_min <- "0.4.30"
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

# Remove any deaths at beginning that were followed by 21 days of no deaths as we have no information in these situations
if(sum(df$deaths>0)>1) {
  if(tail(diff(which(df$deaths>0)),1) > 21) {
    df$deaths[tail(which(df$deaths>0),1)] <- 0
  }
}

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
  n_mcmc <- 250
  n_chains <- 3
  grid_spread <- 2
  sleep <- 2
  start_adaptation <- 50
} else {
  n_particles <- 10
  replicates <- 100
  n_mcmc <- 10000
  n_chains <- 3
  grid_spread <- 11
  sleep <- 120
  start_adaptation <- 1000
}

if (parallel) {
  suppressWarnings(future::plan(future::multiprocess()))
}

# Defualt edges
R0_min <- 1.6
R0_max <- 5.6
Meff_min <- 0.5
Meff_max <- 10
Meff_pl_min <- 0
Meff_pl_max <- 1
Rt_shift_min <- 0
Rt_shift_max <- 0.001
Rt_shift_scale_min <- 0.1
Rt_shift_scale_max <- 10


last_start_date <- as.Date(null_na(min_death_date))-10
first_start_date <- as.Date(null_na(min_death_date))-55

## -----------------------------------------------------------------------------
## Step 2b: Sourcing previous fits to start pmcmc nearby
## -----------------------------------------------------------------------------

# 1. Do we have a previous run for this country
pars_former <- readRDS("pars_init.rds")
pars_former <- pars_former[pars_former$iso3c == iso3c,]

if (nrow(pars_former) == 1) {
  
  R0_start <- pars_former$R0
  date_start <- pars_former$start_date
  Meff_start <- pars_former$Meff
  Meff_pl_start <- pars_former$Meff_pl
  Rt_shift_start <- pars_former$Rt_shift
  Rt_shift_duration <- pars_former$Rt_shift_duration
  Rt_shift_scale_start <- pars_former$Rt_shift_scale
  date_Meff_change <- pars_former$date_Meff_change
  
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
  Rt_shift_start <- 0.5
  Rt_shift_scale_start <- 2
  Rt_shift_duration <- 30
  
  date_Meff_change <- post_lockdown_date_relative(interventions[[iso3c]], 1.05, 
                                     max_date = as.Date("2020-06-02"),
                                     min_date = as.Date("2020-02-01"))
  
}

R0_start <- min(max(R0_start, R0_min), R0_max)
date_start <- min(max(as.Date(date_start), as.Date(first_start_date)), as.Date(last_start_date))
Meff_start <- min(max(Meff_start, Meff_min), Meff_max)
Meff_pl_start <- min(max(Meff_pl_start, Meff_pl_min), Meff_pl_max)
Rt_shift_start <- min(max(Rt_shift_start, Rt_shift_min), Rt_shift_max)
Rt_shift_scale_start <- min(max(Rt_shift_scale_start, Rt_shift_scale_min), Rt_shift_scale_max)


## -----------------------------------------------------------------------------
## Step 2ab: Spline set up
## -----------------------------------------------------------------------------

last_shift_date <- as.Date(date_Meff_change) + 7
remaining_days <- as.Date(date_0) - last_shift_date - 21 # reporting delay in place

# how many spline pars do we need
Rt_rw_duration <- 7
rw_needed <- as.numeric(round(remaining_days/Rt_rw_duration))

# set up rw pars
pars_init_rw <- as.list(rep(0, rw_needed))
pars_min_rw <- as.list(rep(-5, rw_needed))
pars_max_rw <- as.list(rep(5, rw_needed))
pars_discrete_rw <- as.list(rep(FALSE, rw_needed))
names(pars_init_rw) <- names(pars_min_rw) <- names(pars_max_rw) <- names(pars_discrete_rw) <- paste0("Rt_rw_", seq_len(rw_needed))

## -----------------------------------------------------------------------------
## Step 2b: PMCMC parameter set up
## -----------------------------------------------------------------------------

# PMCMC Parameters
pars_init = list('start_date' = date_start, 
                 'R0' = R0_start, 
                 'Meff' = Meff_start, 
                 'Meff_pl' = Meff_pl_start,
                 "Rt_shift" = 0,
                 "Rt_shift_scale" = Rt_shift_scale_start)
pars_min = list('start_date' = first_start_date, 
                'R0' = R0_min, 
                'Meff' = Meff_min, 
                'Meff_pl' = Meff_pl_min,
                "Rt_shift" = Rt_shift_min,
                "Rt_shift_scale" = Rt_shift_scale_min)
pars_max = list('start_date' = last_start_date, 
                'R0' = R0_max, 
                'Meff' = Meff_max, 
                'Meff_pl' = Meff_pl_max,
                "Rt_shift" = Rt_shift_max,
                "Rt_shift_scale" = Rt_shift_scale_max)
pars_discrete = list('start_date' = TRUE, 'R0' = FALSE, 'Meff' = FALSE, 
                     'Meff_pl' = FALSE, "Rt_shift" = FALSE, "Rt_shift_scale" = FALSE)
pars_obs = list(phi_cases = 1, k_cases = 2, phi_death = 1, k_death = 2, exp_noise = 1e6)

# add in the spline list
pars_init <- append(pars_init, pars_init_rw)
pars_min <- append(pars_min, pars_min_rw)
pars_max <- append(pars_max, pars_max_rw)
pars_discrete <- append(pars_discrete, pars_discrete_rw)

# Covriance Matrix
proposal_kernel <- diag(length(names(pars_init))) * 0.3
rownames(proposal_kernel) <- colnames(proposal_kernel) <- names(pars_init)
proposal_kernel["start_date", "start_date"] <- 1.5

# MCMC Functions - Prior and Likelihood Calculation
logprior <- function(pars){
  ret <- dunif(x = pars[["start_date"]], min = -55, max = -10, log = TRUE) +
    dnorm(x = pars[["R0"]], mean = 3, sd = 1, log = TRUE) +
    dnorm(x = pars[["Meff"]], mean = 3, sd = 3, log = TRUE) +
    dunif(x = pars[["Meff_pl"]], min = 0, max = 1, log = TRUE) +
    dnorm(x = pars[["Rt_shift"]], mean = 0, sd = 1, log = TRUE) +
    dunif(x = pars[["Rt_shift_scale"]], min = 0.1, max = 10, log = TRUE)
  
  # get rw spline parameters
  if(any(grepl("Rt_rw", names(pars)))) {
    Rt_rws <- pars[grepl("Rt_rw", names(pars))]
    for (i in seq_along(Rt_rws)) {
      ret <- ret + dnorm(x = Rt_rws[[i]], mean = 0, sd = 0.2, log = TRUE)
      }
  }
  return(ret)
}


# input params
hosp_beds <- squire:::get_hosp_bed_capacity(country)
icu_beds <- squire:::get_ICU_bed_capacity(country)

# Increase ICU beds where known to be too low:
if (iso3c == "BRA") {
  # https://g1.globo.com/bemestar/coronavirus/noticia/2020/06/08/casos-de-coronavirus-e-numero-de-mortes-no-brasil-em-8-de-junho.ghtml - date we predicted ICU to be at capacity and reported to be at 70% 
  icu_beds <- icu_beds / 0.7
} 

# sleep so parallel is chill
Sys.sleep(time = runif(1, 0, sleep))
out_det <- squire::pmcmc(data = data, 
                         n_mcmc = n_mcmc,
                         log_prior = logprior,
                         n_particles = 1,
                         steps_per_day = 1,
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
                         Rt_args = squire:::Rt_args_list(
                           date_Meff_change = date_Meff_change,
                           scale_Meff_pl = TRUE,
                           Rt_shift_duration = 7,
                           Rt_rw_duration = Rt_rw_duration), 
                         burnin = ceiling(n_mcmc/10),
                         seeding_cases = 5,
                         replicates = replicates,
                         required_acceptance_ratio = 0.20,
                         start_adaptation = start_adaptation,
                         baseline_hosp_bed_capacity = hosp_beds, 
                         baseline_ICU_bed_capacity = icu_beds)

## -----------------------------------------------------------------------------
## Step 2b: Summarise Fits for Interface
## -----------------------------------------------------------------------------

## and save the info for the interface
all_chains <- do.call(rbind,lapply(out_det$pmcmc_results$chains, "[[", "results"))
best <- all_chains[which.max(all_chains$log_posterior), ]

## BEST
## -----------------------------------------------------------------------------

# get the R0, betas and times into a data frame
R0 <- best$R0
start_date <- squire:::offset_to_start_date(data$date[1],round(best$start_date))
Meff <- best$Meff
Meff_pl <- best$Meff_pl
Rt_shift <- best$Rt_shift
Rt_shift_scale <- best$Rt_shift_scale

if(!is.null(date_R0_change)) {
  tt_beta <- squire:::intervention_dates_for_odin(dates = date_R0_change,
                                                  change = R0_change,
                                                  start_date = start_date,
                                                  steps_per_day = 1)
} else {
  tt_beta <- 0
}

if(!is.null(R0_change)) {
  R0 <- squire:::evaluate_Rt_pmcmc(R0_change = tt_beta$change, 
                                   date_R0_change = tt_beta$dates, 
                                   R0 = R0, 
                                   pars = as.list(best[1,-(1:2)]),
                                   Rt_args = out_det$pmcmc_results$inputs$Rt_args)
} else {
  R0 <- R0
}
beta_set <- squire:::beta_est(squire_model = squire_model,
                              model_params = out_det$pmcmc_results$inputs$model_params,
                              R0 = R0)


# ## UPPER
# ## -----------------------------------------------------------------------------
# 
# # do the same for the best and worst case
# ft_d <- format_output(out_det, "deaths")
# ft_d <- group_by(ft_d, replicate) %>% summarise(d = sum(y,na.rm=TRUE))
# max_pos <- which.max(ft_d$d)
# min_pos <- which.min(ft_d$d)
# 
# # get the R0, betas and times into a data frame
# R0 <- out_det$replicate_parameters$R0[max_pos]
# start_date <- squire:::offset_to_start_date(data$date[1],round(out_det$replicate_parameters$start_date[max_pos]))
# Meff <- out_det$replicate_parameters$Meff[max_pos]
# Meff_pl <- out_det$replicate_parameters$Meff_pl[max_pos]
# Rt_shift <- out_det$replicate_parameters$Rt_shift[max_pos]
# Rt_shift_scale <- out_det$replicate_parameters$Rt_shift_scale[max_pos]
# 
# if(!is.null(date_R0_change)) {
#   tt_beta_upper <- squire:::intervention_dates_for_odin(dates = date_R0_change,
#                                                   change = R0_change,
#                                                   start_date = start_date,
#                                                   steps_per_day = 1)
# } else {
#   tt_beta_upper <- 0
# }
# 
# if(!is.null(R0_change)) {
#   R0_upper <- squire:::evaluate_Rt_pmcmc(R0_change = tt_beta_upper$change, 
#                                    date_R0_change = tt_beta_upper$dates, 
#                                    R0 = R0, 
#                                    pars = list(
#                                      Meff = Meff,
#                                      Meff_pl = Meff_pl,
#                                      Rt_shift = Rt_shift,
#                                      Rt_shift_scale = Rt_shift_scale
#                                    ),
#                                    Rt_args = out_det$pmcmc_results$inputs$Rt_args)
# } else {
#   R0_upper <- R0
# }
# 
# beta_set_upper <- squire:::beta_est(squire_model = squire_model,
#                               model_params = out_det$pmcmc_results$inputs$model_params,
#                               R0 = R0_upper)
# 
# ## LOWER
# ## -----------------------------------------------------------------------------
# 
# # get the R0, betas and times into a data frame
# R0 <- out_det$replicate_parameters$R0[min_pos]
# start_date <- squire:::offset_to_start_date(data$date[1],round(out_det$replicate_parameters$start_date[min_pos]))
# Meff <- out_det$replicate_parameters$Meff[min_pos]
# Meff_pl <- out_det$replicate_parameters$Meff_pl[min_pos]
# Rt_shift <- out_det$replicate_parameters$Rt_shift[min_pos]
# Rt_shift_scale <- out_det$replicate_parameters$Rt_shift_scale[min_pos]
# 
# if(!is.null(date_R0_change)) {
#   tt_beta_lower <- squire:::intervention_dates_for_odin(dates = date_R0_change,
#                                                         change = R0_change,
#                                                         start_date = start_date,
#                                                         steps_per_day = 1)
# } else {
#   tt_beta_lower <- 0
# }
# 
# if(!is.null(R0_change)) {
#   R0_lower <- squire:::evaluate_Rt_pmcmc(R0_change = tt_beta_lower$change, 
#                                          date_R0_change = tt_beta_lower$dates, 
#                                          R0 = R0, 
#                                          pars = list(
#                                            Meff = Meff,
#                                            Meff_pl = Meff_pl,
#                                            Rt_shift = Rt_shift,
#                                            Rt_shift_scale = Rt_shift_scale
#                                          ),
#                                          Rt_args = out_det$pmcmc_results$inputs$Rt_args)
# } else {
#   R0_lower <- R0
# }
# 
# beta_set_lower <- squire:::beta_est(squire_model = squire_model,
#                                     model_params = out_det$pmcmc_results$inputs$model_params,
#                                     R0 = R0_lower)

## -----------------------------------------------------------------------------


df <- data.frame(tt_beta = tt_beta$tt, beta_set = beta_set, 
                 date = start_date + tt_beta$tt, Rt = R0, 
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
## Step 3a: Fit Stochastic Model
## -----------------------------------------------------------------------------

Sys.setenv("SQUIRE_PARALLEL_DEBUG" = "TRUE")

out <- out_det
out$pmcmc_results$inputs$squire_model <- explicit_model()
out$pmcmc_results$inputs$model_params$dt <- 0.02
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
g1 <- simple_pmcmc_plot(out) + theme(text = element_text(size = 2))

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

line_v <- ggplot() + cowplot::draw_line(y = 0:10, x=1) + 
  theme(panel.background = element_blank(),
        axis.title = element_blank(), 
        axis.text = element_blank(), 
        axis.ticks = element_blank())


header <- cowplot::plot_grid(title, line, ncol = 1)

index <- squire:::odin_index(out$model)
forecast <- 0

suppressMessages(suppressWarnings(
  d <- deaths_plot_single(
    out, data, date = date,
    date_0 = date_0, forecast = forecast, 
    single = TRUE, full = TRUE) + 
    theme(legend.position = "none")))

intervention <- intervention_plot_google(
  interventions[[iso3c]], date, data, 
  forecast, 
  start_date = min(as.Date(out$replicate_parameters$start_date))) + 
  geom_vline(xintercept = as.Date(date_Meff_change) + seq(Rt_rw_duration, Rt_rw_duration*rw_needed, by = Rt_rw_duration),
             linetype = "dashed") + xlab("") + theme(axis.text.x = element_text(angle=45, vjust = 0.5))

rtp <- rt_plot(out)$plot

date_range <- as.Date(c(min(as.Date(out$replicate_parameters$start_date)),date_0))
suppressMessages(suppressWarnings(
  bottom <- cowplot::plot_grid(
    intervention + scale_x_date(date_breaks = "1 month", date_labels = "%b" ,limits = date_range), 
    d + scale_x_date(date_breaks = "1 month", date_labels = "%b" ,limits = date_range),
    rtp + scale_x_date(date_breaks = "1 month", date_labels = "%b" ,limits = date_range),
    ncol=1,
    rel_heights = c(0.4,0.6,0.4))
))

combined <- cowplot::plot_grid(header, 
                               cowplot::plot_grid(g1, line_v, bottom, ncol = 3, rel_widths = c(3,0.1,1)), 
                               ncol=1, rel_heights = c(1,15))



ggsave("fitting.pdf",width=24, height=12, 
       combined)


## Save the grid out object

# remove states to keep object memory save down
for(i in seq_along(out$pmcmc_results$chains)) {
  out$pmcmc_results$chains[[i]]$states <- NULL
  out$pmcmc_results$chains[[i]]$covariance_matrix <- NULL
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
  current <- squire:::t0_variables(out)
  current <- unlist(lapply(current, "[[", "R0"))
  
  wanted <-  vapply(seq_along(out$replicate_parameters$R0), function(y){
    
    if(!is.null(out$interventions$date_R0_change)) {
      tt_beta <- squire:::intervention_dates_for_odin(dates = out$interventions$date_R0_change,
                                                      change = out$interventions$R0_change,
                                                      start_date = out$replicate_parameters$start_date[y],
                                                      steps_per_day = 1/out$parameters$dt)
    } else {
      tt_beta <- 0
    }
    
    Rt <- tail(squire:::evaluate_Rt_pmcmc(
      R0_change = c(tt_beta$change,R0_ch), 
      date_R0_change = c(
        tt_beta$dates, 
        tail(out$pmcmc_results$inputs$data$date,1)+1),
      R0 = out$replicate_parameters$R0[y], 
      pars = as.list(out$replicate_parameters[1,-(1:2)]),
      Rt_args = out$pmcmc_results$inputs$Rt_args) ,1)
    
    
  }, numeric(1))
  
  mean(wanted/current)
}

# Maintaining the current set of measures for a further 3 months  
maintain_3months_lift <- squire::projections(out, 
                                             R0_change = c(1), 
                                             tt_R0 = c(0), 
                                             time_period = time_period)

# Enhancing movement restrictions for 3 months (50% further reduction in contacts) (Mitigation) 
mitigation_3months_lift <- squire::projections(out, 
                                               R0_change = c(0.5), 
                                               tt_R0 = c(0), 
                                               time_period = time_period)

# Relax by 50% for 3 months 
reverse_3_months_lift <- squire::projections(out, 
                                             R0_change = 1.5, 
                                             tt_R0 = c(0), 
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
                                                    R0_change = c(1), 
                                                    tt_R0 = c(0), 
                                                    time_period = time_period)

mitigation_3months_lift_surged <- squire::projections(out_surged, 
                                                      R0_change = c(0.5), 
                                                      tt_R0 = c(0), 
                                                      time_period = time_period)

reverse_3_months_lift_surged <- squire::projections(out_surged, 
                                                    R0_change = c(1.5), 
                                                    tt_R0 = c(0), 
                                                    time_period = time_period)

## -----------------------------------------------------------------------------
## 4.3. Extra Scenarios
## N.B. Not implemented but leave here if we want to add more ------------------
## -----------------------------------------------------------------------------

if (full_scenarios) {
  
  
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
                 var_select = c("infections","deaths","hospital_demand",
                                "ICU_demand", "D", "hospital_incidence","ICU_incidence"),
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
  
  # add in cumulative deaths
  cum_i <- pd %>%
    filter(compartment == "deaths") %>%
    group_by(replicate) %>%
    mutate(y = cumsum(y),
           compartment = "cumulative_deaths") %>%
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
data_sum$version <- "v4"
data_sum <- dplyr::mutate(data_sum, across(dplyr::starts_with("y_"), ~round(.x,digits = 2)))
data_sum <- data_sum[]
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
