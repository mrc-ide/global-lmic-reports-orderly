orderly_id <- tryCatch(orderly::orderly_run_info()$id,
                       error = function(e) "<id>")

print(sessionInfo())
RhpcBLASctl::blas_set_num_threads(1L)
RhpcBLASctl::omp_set_num_threads(1L)

## -----------------------------------------------------------------------------
## Step 1: Incoming Date
## -----------------------------------------------------------------------------
system(paste0("echo Vaccine Reports for  ",iso3c))
if(!identical(seed, FALSE)){
  set.seed(seed)
}

# format user provided arguments correctly
date <- as.Date(date)

#set real end date to 2022-01-01
date <- as.Date("2022-01-01")

country <- squire::population$country[match(iso3c, squire::population$iso3c)[1]]

## Get the excess mortality estimates from the economist
excess_deaths <- readRDS("excess_deaths.Rds") %>%
  rename(iso = iso3c) %>%
  filter(iso == iso3c) %>%
  arrange(date_start) %>%
  filter(date_end <= date)

if(nrow(excess_deaths) == 0){
  excess_deaths <-
    data.frame(
      iso = iso3c,
      date_start = seq(-7, 7, by = 7) + date,
      date_end = seq(0, 14, by = 7) + date,
      deaths = 0
    )
}

death_limit <- 10

fit_excess <- sum(excess_deaths$deaths) > death_limit

#setup model and parameters

#get model start dates
if(fit_excess){
  if(iso3c == "TKM"){
    #remove starting deaths before major waves
    if(sum(cumsum(excess_deaths$deaths) < 40) > 15){
      excess_deaths <- excess_deaths %>%
        filter(cumsum(deaths) > 40)
    }
  }

  first_report <- which(excess_deaths$deaths>0)[1]
  missing <- which(excess_deaths$deaths == 0 | is.na(excess_deaths$deaths))
  to_remove <- missing[missing<first_report]
  if(length(to_remove) > 0) {
    if(length(to_remove) == (nrow(excess_deaths)-1)) {
      excess_deaths <- excess_deaths[-head(to_remove,-1),]
    } else {
      excess_deaths <- excess_deaths[-to_remove,]
    }
  }
  excess_start_date <-  excess_deaths$date_start[1] - 30
  first_start_date <- excess_start_date
}

pop <- squire::get_population(country)
squire_model <- squire.page:::nimue_booster_min_model()
default_parameters_func <- squire_model$parameter_func
define_parameters_func <- function(default_parameters_func){
  function(...){
    args <- list(...)
    S_0 <- args$S_0
    args$S_0 <- NULL
    odin_pars <- do.call(default_parameters_func, args)
    if(!is.null(S_0)){
      odin_pars$E1_0 <- array(0, dim = dim(odin_pars$E1_0))
      odin_pars$S_0 <- S_0
    }
    odin_pars
  }
}
end_date <- date
#Invariant parameters
parameters <- list(
  country = country,
  hosp_bed_capacity = squire:::get_hosp_bed_capacity(country),
  ICU_bed_capacity = squire:::get_ICU_bed_capacity(country)
)
# Increase ICU beds where known to be too low:
if (iso3c == "BRA") {
  # https://g1.globo.com/bemestar/coronavirus/noticia/2020/06/08/casos-de-coronavirus-e-numero-de-mortes-no-brasil-em-8-de-junho.ghtml - date we predicted ICU to be at capacity and reported to be at 70%
  parameters$hosp_bed_capacity <- parameters$hosp_bed_capacity / 0.7
}
#load inputs
vacc_inputs <- get_vaccine_inputs(iso3c)
tt_vaccine <- as.numeric(c(first_start_date, vacc_inputs$date_vaccine_change) - first_start_date)
parameters$primary_doses <- c(0, vacc_inputs$primary_doses)
parameters$tt_primary_doses <- tt_vaccine
parameters$second_dose_delay <- vacc_inputs$second_dose_delay
parameters$booster_doses <- c(0, vacc_inputs$booster_doses)
parameters$tt_booster_doses <- tt_vaccine
rm(tt_vaccine)
parameters$vaccine_coverage_mat <- vacc_inputs$vaccine_coverage_mat
parameters$rel_infectiousness_vaccinated <- c(0.5)
parameters$vaccine_booster_follow_up_coverage <- c(rep(0, 12), rep(1, 5))
parameters$protection_delay_rate <- 1/7
parameters$protection_delay_shape <- 2
parameters$protection_delay_time <- as.numeric(date - first_start_date)
parameters$time_period <- as.numeric(date - first_start_date) + 1
#Variant dependant parameters
variants_to_model <- c("Delta")
#load inputs
sample_vaccine_efficacies <- readRDS("vaccine_params.Rds")$sample_vaccine_efficacies
variant_timings <- readRDS("variant_timings.Rds")[[iso3c]] %>%
  filter(variant %in% variants_to_model & !is.na(start_date)) %>%
  select(!mid_date) %>%
  arrange(.data$start_date) %>%
  mutate(end_date = if_else(
    .data$end_date > lead(.data$start_date, n = 1, default = as_date(max(.data$end_date) + 1)),
    lead(.data$start_date, 1, default = as_date(max(.data$end_date) + 1)),
    .data$end_date
  ))
variants_to_model <- variant_timings$variant
parameters <-
  append(parameters, estimate_healthcare_durations(variant_timings, first_start_date))
#Sample from random parameter
#generate samples
dur_R <- sample_duration_natural_immunity(samples)
ifr <- sample_ifr(samples, iso3c)
immune_escape <- sample_variant_immune_escape(samples, variants_to_model)
prob_hosp_multiplier <- sample_variant_prob_hosp(samples, variants_to_model)
prob_severe_multiplier <- sample_variant_prob_severe(samples, variants_to_model)
variant_ve <- sample_vaccine_efficacies(samples, names(vacc_inputs$platforms)[as.logical(vacc_inputs$platforms[1,])])
#format into correct setup
distribution <- map(seq_len(samples), function(x){
  pars <- list()
  ## disease parameters:
  pars$prob_hosp <- ifr$prob_hosp[[x]]
  pars$prob_severe <- ifr$prob_severe[[x]]
  pars$prob_non_severe_death_treatment <- ifr$prob_non_severe_death_treatment[[x]]
  pars$prob_severe_death_treatment <- ifr$prob_severe_death_treatment[[x]]
  pars$prob_severe_death_no_treatment <- ifr$prob_severe_death_no_treatment[[x]]
  pars$prob_non_severe_death_no_treatment <- ifr$prob_non_severe_death_no_treatment[[x]]
  ## variant parameters:
  dur_R_change <- variant_immune_escape(variant_timings, immune_escape, x,
                                        dur_R, first_start_date)
  pars$dur_R <- dur_R_change$var
  pars$tt_dur_R <- dur_R_change$tt
  prob_hosp_multiplier_change <- multiplier_changes_over_time(variant_timings,
                                                              prob_hosp_multiplier, x, first_start_date)
  pars$prob_hosp_multiplier <- prob_hosp_multiplier_change$var
  pars$tt_prob_hosp_multiplier <- prob_hosp_multiplier_change$tt
  prob_severe_multiplier_change <- multiplier_changes_over_time(variant_timings,
                                                                prob_severe_multiplier, x, first_start_date)
  pars$prob_severe_multiplier <- prob_severe_multiplier_change$var
  pars$tt_prob_severe_multiplier <- prob_severe_multiplier_change$tt
  #Vaccine Efficacies
  vacc_inputs <- vaccine_eff_over_time(variant_timings, variant_ve, x, first_start_date)
  pars$dur_V <- vacc_inputs$dur_V
  pars$vaccine_efficacy_infection <- vacc_inputs$vaccine_efficacy_infection
  pars$vaccine_efficacy_disease <- vacc_inputs$vaccine_efficacy_disease
  pars$tt_dur_V <- pars$tt_vaccine_efficacy_infection <-
    pars$tt_vaccine_efficacy_disease <- vacc_inputs$tt
  pars
})


## Excess deaths loop
if(fit_excess){
  #ensure parameters are in correct format
  difference_in_t <- as.numeric(excess_start_date - first_start_date)
  excess_parameters <- update_parameters(parameters, difference_in_t)
  excess_distribution <- update_distribution(distribution, difference_in_t)
  excess_squire_model <- squire_model
  if(!is.null(excess_parameters$prefit_vaccines)){
    #if vaccinations occur before epidemic, we run the model with just vaccinations
    #to get the current state of vaccination status at the epidemic start
    excess_distribution <- prefit_vaccines(excess_parameters, excess_distribution, excess_squire_model)
    excess_parameters$prefit_vaccines <- NULL
    excess_squire_model$parameter_func <- define_parameters_func(default_parameters_func)
  }
  #fitting parameters
  #load in country specific default parameters
  fitting_params <- readRDS("fitting_params.Rds")[[iso3c]]$excess
  #check if we need to update parameters
  if(identical(initial_infections_interval, "NULL")){
    initial_infections_interval <- fitting_params$initial_infections_interval
  }
  if(identical(n_particles, "NULL")){
    n_particles <- fitting_params$n_particles
  }
  if(identical(k, "NULL")){
    k <- fitting_params$k
  }
  if(identical(rt_interval, "NULL")){
    rt_interval <- fitting_params$rt_interval
  }
  trimming <- fitting_params$trimming

  #run model
  #use parallel if asked for
  if(parallel){
    future::plan(future::multisession())
  }

  excess_out <- rt_optimise(
    data = excess_deaths,
    distribution = excess_distribution,
    squire_model = excess_squire_model,
    parameters = excess_parameters,
    start_date = excess_start_date,
    parallel = parallel,
    rt_spacing = 14,
    initial_infections_interval = initial_infections_interval,
    n_particles = n_particles,
    k = k,
    rt_interval = rt_interval
  )

  excess_out <- trim_output(excess_out, trimming)

  #Stop using parallel, furrr doesn't like something (maybe model object)
  Sys.setenv(SQUIRE_PARALLEL_DEBUG = "TRUE")
  if(parallel){
    future::plan(future::sequential())
  }

  #save fitting plot
  summarise_fit("excess_fitting.pdf", excess_out, country, iso3c, end_date, excess_start_date)

  save_output(excess_out, "excess_out.Rds")
} else {
  file.create(c("excess_fitting.pdf", "excess_out.Rds"))
}
