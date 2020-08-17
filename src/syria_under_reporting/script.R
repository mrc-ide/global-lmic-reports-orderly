orderly_id <- tryCatch(orderly::orderly_run_info()$id,
                       error = function(e) "<id>") 

print(sessionInfo())
RhpcBLASctl::blas_set_num_threads(1L)
RhpcBLASctl::omp_set_num_threads(1L)

version_min <- "0.4.31"
if(packageVersion("squire") < version_min) {
  stop("squire needs to be updated to at least ", version_min)
}

system(paste0("echo Syria Reporting. Reporting Fraction = ",reporting_fraction, ". Date = ", date, 
              ". Urban = ", urban, ". Poorer Health Outcomes = ", poorer_health_outcomes ,
              ". Younger Cities = ", younger_cities))

iso3c <- "SYR"
country <- "Syria"

## -----------------------------------------------------------------------------
## Step 1: Incoming Date
## -----------------------------------------------------------------------------
set.seed(123)
date <- as.Date(date)

## Get the worldometers data from JHU
data <- as.data.frame(data.table::fread(
  "https://raw.githubusercontent.com/CSSEGISandData/COVID-19/master/csse_covid_19_data/csse_covid_19_time_series/time_series_covid19_deaths_global.csv"
))

# format into the correct style for squire
data <- data %>% filter(`Country/Region` == "Syria") %>% tidyr::pivot_longer(matches("^\\d"))
names(data) <- c("", "country","lat","lon","date","deaths")
data <- data[,c("date","deaths")]
data$date <- as.Date(data$date, format = "%m/%d/%y")

# and into daily deaths
data$deaths <- c(0, diff(data$deaths))
data <- filter(data, date <= as.Date(date))

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

# Maintaining daily sheet for deaths
sheet <- readxl::read_xlsx("death_sheet.xlsx")
sheet$Dam[is.na(sheet$Dam)] <- 0
if (sheet$Deaths[1] != sum(data$deaths)) {
  stop("Death Sheet out of sync with JHU data stream")
}

# use this to remove non damascus deaths
non_dam <- sheet$Daily - sheet$Dam
non_dam_match <- non_dam[match(data$date, as.Date(sheet$Date))]
non_dam_match[is.na(non_dam_match)] <- 0
data$deaths <- data$deaths - non_dam_match

# Then there is 1 death in Aleppo and 1 in Homs not known as to what date so remove from day with most deaths for which we don't know where they cam
data$deaths[which.max(data$deaths)] <- data$deaths[which.max(data$deaths)]-1
data$deaths[which.max(data$deaths)] <- data$deaths[which.max(data$deaths)]-1

# dat_0 is just the current date now
date_0 <- date

# get country data
interventions <- readRDS("google_brt.rds")

## -----------------------------------------------------------------------------
## Step 2: Times for the interace
## -----------------------------------------------------------------------------

## -----------------------------------------------------------------------------
## Step 2a: Set up args for pmcmc/grid
## -----------------------------------------------------------------------------

# what is the date of first death
null_na <- function(x) {if(is.null(x)) {NA} else {x}}
min_death_date <- data$date[which(data$deaths>0)][1]

# calibration arguments
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

# pmcmc arguments 
n_particles <- 10
replicates <- 100
n_mcmc <- 2000
n_chains <- 3
start_adaptation <- 500

# this should be in parallel
suppressWarnings(future::plan(future::multiprocess()))

# Defualt edges to seatch within
R0_min <- 1.6
R0_max <- 5.6
Meff_min <- 0.1
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

# 1. Get previous run for this country
pars_former <- readRDS("pars_init.rds")
pars_former <- pars_former$SYR

# 2. Get starting parameters 
R0_start <- pars_former$R0
date_start <- pars_former$start_date
Meff_start <- pars_former$Meff
Meff_pl_start <- pars_former$Meff_pl
Rt_shift_start <- pars_former$Rt_shift
Rt_shift_duration <- pars_former$Rt_shift_duration
Rt_shift_scale_start <- pars_former$Rt_shift_scale
Rt_rw_duration <- pars_former$Rt_rw_duration
date_Meff_change <- pars_former$date_Meff_change

# 3. Range base them
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
rw_needed <- as.numeric(round(remaining_days/Rt_rw_duration))

# set up rw pars
if (is.null(pars_former)) {
  pars_init_rw <- as.list(rep(0, rw_needed))
} else {
  pars_init_rw <- as.list(pars_former[grep("Rt_rw_\\d",names(pars_former))])
  if(length(pars_init_rw) < rw_needed) {
    pars_init_rw[[rw_needed]] <- 0 
  } else if(length(pars_init_rw) > rw_needed) {
    pars_init_rw <- head(pars_init_rw, rw_needed)
  }
}

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
pars_obs = list(phi_cases = 1, k_cases = 2, phi_death = 0.1, k_death = 2, exp_noise = 1e6)

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

## -----------------------------------------------------------------------------
## Step 3: Syria specific adjustments and IFR/Pop adjustments
## -----------------------------------------------------------------------------

# input params
hosp_beds <- squire:::get_hosp_bed_capacity("Syria")
icu_beds <- squire:::get_ICU_bed_capacity("Syria")

# population for Damascus with demographics of Syria
# dam_pop_city <- 1569394
dam_pop_urban <- 2396587 # https://populationstat.com/syria/damascus
# https://en.wikipedia.org/wiki/Rif_Dimashq_Governorate (2011) - but population growth for Damscus urban has stayed same since 2011
dam_pop_rural <- 2836000 

pop <- squire::get_population("Syria")

# matrix
mix_mat <- get_mixing_matrix("Syria")

# https://eprints.lse.ac.uk/103841/1/CRP_covid_19_in_Syria_policy_memo_published.pdf
icu_beds <- 96
hosp_beds <- 1920

# scan across a range of undereporting and IFR and assumptions about the excess mortality

#reporting_fraction <- 0.01
reporting_fraction <- as.numeric(reporting_fraction)

# poorer_health_outcomes <- TRUE
poorer_health_outcomes <- as.logical(poorer_health_outcomes)
prob_nsdt <- squire:::probs$prob_non_severe_death_treatment
if(poorer_health_outcomes) {
  prob_nsdt[seq_len(length(prob_nsdt)-1)] <- prob_nsdt[length(prob_nsdt)-1]
}

# urban <- TRUE
urban <- as.logical(urban)
if (urban) {
  dam_pop <- dam_pop_urban
} else {
  dam_pop <- dam_pop_rural
}

# younger_cities <- TRUE
younger_cities <- as.logical(younger_cities)
if (younger_cities) {
  altered_pop <- pop$n * (dexp(1:17, 1/20)/mean(dexp(1:17, 1/20)))
  altered_pop <- altered_pop / sum(altered_pop)
  population <- round(altered_pop*dam_pop)
} else {
  population <- round((pop$n/sum(pop$n))*dam_pop)
}

## -----------------------------------------------------------------------------
## Step 3: Run pmcmc
## -----------------------------------------------------------------------------

res <- squire::pmcmc(data = data, 
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
                     population = population,
                     baseline_contact_matrix = mix_mat,
                     R0_change = R0_change,
                     date_R0_change = date_R0_change,
                     Rt_args = squire:::Rt_args_list(
                       date_Meff_change = date_Meff_change,
                       scale_Meff_pl = TRUE,
                       Rt_shift_duration = 7,
                       Rt_rw_duration = Rt_rw_duration), 
                     prob_non_severe_death_treatment = prob_nsdt,
                     burnin = ceiling(n_mcmc/10),
                     reporting_fraction = reporting_fraction,
                     seeding_cases = 5,
                     replicates = replicates,
                     required_acceptance_ratio = 0.20,
                     start_adaptation = start_adaptation,
                     baseline_hosp_bed_capacity = hosp_beds, 
                     baseline_ICU_bed_capacity = icu_beds)

# redraw using stochastic
res$pmcmc_results$inputs$squire_model <- explicit_model()
res$pmcmc_results$inputs$model_params$dt <- 0.02
pmcmc <- res$pmcmc_results
res <- generate_draws_pmcmc(pmcmc = pmcmc,
                            burnin = ceiling(n_mcmc/10),
                            n_chains = n_chains,
                            squire_model = res$pmcmc_results$inputs$squire_model,
                            replicates = replicates,
                            n_particles = n_particles,
                            forecast = 0,
                            country = country,
                            population = res$parameters$population,
                            interventions = res$interventions,
                            data = res$pmcmc_results$inputs$data)

# Add the prior
res$pmcmc_results$inputs$prior <- as.function(c(formals(logprior), 
                                                body(logprior)), 
                                              envir = new.env(parent = environment(stats::acf)))

# remove states to keep object memory save down
for(i in seq_along(res$pmcmc_results$chains)) {
  res$pmcmc_results$chains[[i]]$states <- NULL
  res$pmcmc_results$chains[[i]]$covariance_matrix <- NULL
}

## -----------------------------------------------------------------------------
## Step 4: Comparisons to the data
## -----------------------------------------------------------------------------

## Damascus death sources
## https://www.facebook.com/damascusgovrnorat/posts/191862302286989
reported <- data.frame("date" = seq.Date(as.Date("2020-07-25"), as.Date("2020-08-01"),1),
                       "deaths" = c(78, 88, 87, 108, 133, 127, 104, 107))

# https://aliqtisadi.com/790464-%D9%8A%D9%88%D8%AC%D8%AF-%D9%81%D9%8A-%D8%AF%D9%85%D8%B4%D9%82-%D8%B3%D8%AA-%D9%85%D9%82%D8%A7%D8%A8%D8%B1/
reported$deaths_low <- reported$deaths - 50
reported$deaths_high <- reported$deaths - 15

# https://www.facebook.com/MEENALMASOOL/photos/a.1277434595713288/3039607002829363/?type=3
reported$deaths <- reported$deaths -25

## Let's also compare to an assumption that there is even more excess mortlaity than historic due to secondary pressure on health systems
reported$extra_deaths <- reported$deaths - 20
reported$extra_deaths_low <- reported$deaths_low - 20
reported$extra_deaths_high <- reported$deaths_high - 20

# get the model run deaths for these dates:
deaths <- squire::format_output(res, "deaths", date_0 = date)
model_deaths <- deaths$y[deaths$date %in% reported$date]
all_model_deaths <- deaths$y[deaths$date %in% data$date]

# get likelihoods against the excess deaths
ll_reported <- squire:::ll_nbinom(model = model_deaths, 
                         data = rep(reported$deaths, nrow(res$replicate_parameters)), 
                         phi = 1, k = 1, 
                         exp_noise = res$pmcmc_results$inputs$pars_obs$exp_noise)

ll_extra <- squire:::ll_nbinom(model = model_deaths, 
                                  data = rep(reported$extra_deaths, nrow(res$replicate_parameters)), 
                                  phi = 1, k = 1, 
                                  exp_noise = res$pmcmc_results$inputs$pars_obs$exp_noise)

ll_all <- squire:::ll_nbinom(model = all_model_deaths, 
                                  data = rep(data$deaths, nrow(res$replicate_parameters)), 
                                  phi = res$pmcmc_results$inputs$pars_obs$phi_death, 
                                  k = 1, 
                                  exp_noise = res$pmcmc_results$inputs$pars_obs$exp_noise)

# how much do they appear within the "range"
within_reported_range <- model_deaths < reported$deaths_high & model_deaths > reported$deaths_low
within_extra_range <- model_deaths < reported$extra_deaths_high & model_deaths > reported$extra_deaths_low

# does the CI capture the point estimate
ci_deaths <- group_by(deaths[deaths$date %in% reported$date,], date) %>% 
  summarise(low = quantile(y, 0.025), high = quantile(y, 0.975))

# summarise these model fit metrics
model_fit_summary <- data.frame("ll_reported" = mean(ll_reported),
                               "ll_extra" = mean(ll_extra),
                               "ll_all" = mean(ll_all),
                               "within_reported_range" = mean(within_reported_range),
                               "within_extra_range" = mean(within_extra_range),
                               "reporting_fraction" = res$pmcmc_results$inputs$pars_obs$phi_death,
                               "range_includes_reported" = mean(ci_deaths$low < reported$deaths & ci_deaths$high > reported$deaths),
                               "range_includes_extra" = mean(ci_deaths$low < reported$extra_deaths & ci_deaths$high > reported$extra_deaths),
                               "urban" = urban,
                               "poorer_health_outcomes" = poorer_health_outcomes,
                               "younger_cities" = younger_cities)

## -----------------------------------------------------------------------------
## Step 5: Quick forward simulation and then save
## -----------------------------------------------------------------------------
res_long <- squire::projections(res, 
                                        R0_change = c(1), 
                                        tt_R0 = c(0), 
                                        time_period = 120)
res_long$projection_args$r <- NULL


saveRDS(res_long, "res.rds")
saveRDS(model_fit_summary, "likelihood.rds")