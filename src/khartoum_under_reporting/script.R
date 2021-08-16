orderly_id <- tryCatch(orderly::orderly_run_info()$id,
                       error = function(e) "<id>") 

print(sessionInfo())
RhpcBLASctl::blas_set_num_threads(1L)
RhpcBLASctl::omp_set_num_threads(1L)

version_min <- "0.4.31"
if(packageVersion("squire") < version_min) {
  stop("squire needs to be updated to at least ", version_min)
}

system(paste0("echo Khartoum Reporting. Reporting Fraction = ",reporting_fraction, ". Date = ", date, 
              ". Pop = ", pop, ". Poorer Health Outcomes = ", poorer_health_outcomes ,
              ". City Age = ", city_age, ". Hospital Use = ", hospital_normal_use  ))

iso3c <- "SDN"
country <- "Sudan"

## -----------------------------------------------------------------------------
## Step 1: Incoming Date
## -----------------------------------------------------------------------------
set.seed(123)
data <- read.csv("daily_deaths.csv")
data$date <- as.Date(as.character(data$date), "%d/%m/%Y")
names(data)[2:3] <- c("cumu", "deaths") 

## We also have 46 deaths betwwen the 1st and the 15th
additional <- data.frame(
  "date" = seq.Date(as.Date("2020-11-01"), as.Date("2020-11-15"), 1),
  "cumu" = 0, 
  "deaths" = round(exp((1:15)*0.125)))

additional$deaths[additional$date == "2020-11-07"] <- additional$deaths[additional$date == "2020-11-07"] + 1
data <- data[-which(data$date == "2020-11-07"), ]

data <- rbind(data, additional)
data <- data[order(data$date), ]
data$cumu <- cumsum(data$deaths)

# date_0 is max date in data
date_0 <- max(data$date)

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
R0_change <- R0_change/max(R0_change, na.rm = TRUE) 
date_R0_change <- interventions[[iso3c]]$date
R0_change <- R0_change[as.Date(date_R0_change) <= date_0]
date_R0_change <- date_R0_change[as.Date(date_R0_change) <= date_0]
date_contact_matrix_set_change <- NULL
R0_change[is.na(R0_change)] <- tail(R0_change[!is.na(R0_change)], 1)

# pmcmc arguments 
n_particles <- 2 # doesn't do anything because using the deterministic version
replicates <- 100
n_mcmc <- 2000
n_chains <- 3
start_adaptation <- 1000

# this should be in parallel
suppressWarnings(future::plan(future::multiprocess()))

# Defualt edges to seatch within
R0_min <- 1.6
R0_max <- 5.6
Meff_min <- -10
Meff_max <- 10
Meff_pl_min <- 0
Meff_pl_max <- 1
Rt_shift_min <- 0
Rt_shift_max <- 0.001
Rt_shift_scale_min <- 0.1
Rt_shift_scale_max <- 10
last_start_date <- as.Date(null_na(min_death_date))-10
first_start_date <- as.Date(null_na(min_death_date))-60

## -----------------------------------------------------------------------------
## Step 2b: Sourcing previous fits to start pmcmc nearby
## -----------------------------------------------------------------------------

# 1. Get previous run for this country
pars_former <- readRDS("pars_init.rds")
pars_former <- pars_former[[iso3c]]

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
pars_init_rw <- as.list(pars_former[grep("Rt_rw_\\d",names(pars_former))])
if(length(pars_init_rw) < rw_needed) {
  pars_init_rw[[rw_needed]] <- 0 
} else if(length(pars_init_rw) > rw_needed) {
  pars_init_rw <- head(pars_init_rw, rw_needed)
}

# pars_min_rw <- as.list(rep(-5, rw_needed))
# pars_max_rw <- as.list(rep(5, rw_needed))
pars_min_rw <- as.list(rep(-0.0001, rw_needed))
pars_max_rw <- as.list(rep(0.0005, rw_needed))
pars_init_rw <- as.list(rep(0, rw_needed))

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
  ret <- dunif(x = pars[["start_date"]], min = -60, max = -10, log = TRUE) +
    dnorm(x = pars[["R0"]], mean = 3, sd = 1, log = TRUE) +
    dnorm(x = pars[["Meff"]], mean = 0, sd = 3, log = TRUE) +
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
## Step 3: Country specific adjustments and IFR/Pop adjustments
## -----------------------------------------------------------------------------

# input params

# Sudan demgraphy
dem <- squire::get_population("Sudan")

# matrix
mix_mat <- squire::get_mixing_matrix("Sudan")

# ------------------------------
# POPULATION SIZE
# ------------------------------

# https://www.macrotrends.net/cities/22579/khartoum/population - City
# pop <- 5829000 
# https://reports.unocha.org/en/country/sudan/card/ZXsJqLs913/ - State
pop <- 8877147
# pop <- as.numeric(pop)

# ------------------------------
# BEDS
# ------------------------------

# mid April/May Covid treatament in icus were provided in a few gvt facilities. This included icu access.
# Anecdotally though there is a sense that this was at peak in April May June and since have started to come down.

# https://go.ifrc.org/reports/13075 # gives a live upper estimate
# https://go.ifrc.org/reports/13042
# https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6132516/ - "74-bed ICUs are almost fully occupied year-round"
icu_beds <- 75
# icu_beds <- as.numeric(icu_beds)

# Hospital beds is again harder
# https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6132516/ - 1737 across the 4 government tertiary care hospitals
# per capita reference https://sudannextgen.com/wp-content/uploads/2019/10/Statistics-on-Hospital-sector-in-Sudan-2.pdf has much higher

# hosp_beds <- 4640
hosp_beds <- 2000
hosp_beds <- as.numeric(hosp_beds)

# ------------------------------
# NORMAL OCCUPANCY
# ------------------------------

# Then what level of beds are in use for non covid

hospital_normal_use <- 0.5
hospital_normal_use <- as.numeric(hospital_normal_use)

# ------------------------------
# PARAMS TO SCAN
# ------------------------------

# scan across a range of undereporting and IFR and assumptions about the excess mortality
# reporting_fraction <- 0.05
reporting_fraction <- as.numeric(reporting_fraction)

poorer_health_outcomes <- TRUE
poorer_health_outcomes <- as.logical(poorer_health_outcomes)
prob_nsdt <- squire:::probs$prob_non_severe_death_treatment
if(poorer_health_outcomes) {
  prob_nsdt[seq_len(length(prob_nsdt)-1)] <- prob_nsdt[length(prob_nsdt)-1]
}

city_age <- "older"
if (city_age == "younger") {
  altered_pop <- dem$n * (dexp(1:17, 1/50)/mean(dexp(1:17, 1/500)))
  altered_pop <- altered_pop / sum(altered_pop)
  population <- round(altered_pop*pop)
} else if (city_age == "older") {
  altered_pop <- dem$n * rev((dexp(1:17, 1/50))/mean(dexp(1:17, 1/50)))
  altered_pop <- altered_pop / sum(altered_pop)
  population <- round(altered_pop*pop)
} else {
  population <- round((dem$n/sum(dem$n))*pop)
}


## -----------------------------------------------------------------------------
## Step 3: Run pmcmc
## -----------------------------------------------------------------------------

rfs <- c(0.005, 0.01, 0.02, 0.03, 0.05, 0.075, 0.1)
results <- pbapply::pblapply(rfs, function(x) {
  
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
                     reporting_fraction = x, 
                     treated_deaths_only = FALSE,
                     seeding_cases = 5,
                     replicates = replicates,
                     required_acceptance_ratio = 0.20,
                     start_adaptation = start_adaptation,
                     baseline_hosp_bed_capacity = round(hosp_beds*(1-hospital_normal_use)), 
                     baseline_ICU_bed_capacity = icu_beds) # ICU bed figure is already so low that it fills up immediately anyway

# Add the prior
res$pmcmc_results$inputs$prior <- as.function(c(formals(logprior), 
                                                body(logprior)), 
                                              envir = new.env(parent = environment(stats::acf)))

# remove states to keep object memory save down
for(i in seq_along(res$pmcmc_results$chains)) {
  res$pmcmc_results$chains[[i]]$states <- NULL
  res$pmcmc_results$chains[[i]]$covariance_matrix <- NULL
}

return(res)

})
saveRDS(results, file.path(here::here(), "analysis/khartoum/results_11_22.rds"))

## -----------------------------------------------------------------------------
## Step 4: Comparisons to the data
## -----------------------------------------------------------------------------

## -----------------------------------------------------------------------------
## Cumulative symptomatic attack rate
## -----------------------------------------------------------------------------

ll <- pbapply::pblapply(results, function(res) {
  
  ## Symptomatic infection prevalence estimate as of 2nd June which we compare to cum_hosp_inc + cum_icu_inc
  # https://www.researchgate.net/publication/342122988_Estimation_of_Coronavirus_COVID-19_Infections_in_Khartoum_State_Sudan
  p_e_low <- 0.083
  p_e_med <- 0.11
  p_e_high <- 0.137
  
  # first lets get symptomatic totals for individuals sufficiently severe to need treatment beofre 2nd June and older than 15
  hosp_inc <- squire::format_output(res, c("hospital_incidence"), date_0 = date_0, reduce_age = FALSE) %>% 
    filter(age_group > 3)
  
  icu_inc <- squire::format_output(res, c("ICU_incidence"), date_0 = date_0, reduce_age = FALSE) %>% 
    filter(age_group > 3)
  
  # and then also total infs
  infs <- squire::format_output(res, c("infections"), date_0 = date_0, reduce_age = FALSE) %>% 
    filter(age_group > 3)
  
  # asymptomatic ratios from Davies et al
  davies <- read.csv("davies_clin_susc.csv")
  davies <- davies[davies$Parameter == "Clinical fraction",]
  age_matches <- unlist(lapply(sapply(sort(unique(infs$age_group)), grep, davies$Group), "[[", 1))
  ages_matched <- sort(unique(infs$age_group))
  
  # estimate of symptomatic numbers per age
  infs$non_severe <- infs$y - hosp_inc$y - icu_inc$y
  infs$symptomatic <- hosp_inc$y + icu_inc$y + infs$non_severe*davies$Q50[age_matches[match(infs$age_group, ages_matched)]]
  infs$symptomatic_low <- hosp_inc$y + icu_inc$y + infs$non_severe*davies$Q025[age_matches[match(infs$age_group, ages_matched)]]
  infs$symptomatic_high <- hosp_inc$y + icu_inc$y + infs$non_severe*davies$Q975[age_matches[match(infs$age_group, ages_matched)]]
  
  # summarise back up to date
  infs_sum <- infs %>% group_by(replicate, t, date) %>% 
    summarise(symptomatic = sum(symptomatic, na.rm = TRUE),
              symptomatic_low = sum(symptomatic_low, na.rm = TRUE),
              symptomatic_high = sum(symptomatic_high, na.rm = TRUE))
  infs_sum$p_e <- infs_sum$symptomatic / sum(res$parameters$population[-c(1:3)])
  infs_sum$p_e_low <- infs_sum$symptomatic_low / sum(res$parameters$population[-c(1:3)])
  infs_sum$p_e_high <- infs_sum$symptomatic_high / sum(res$parameters$population[-c(1:3)])
  
  infs_sum <- infs_sum %>% group_by(replicate) %>% 
    mutate(cumu_symp_pe_med = cumsum(p_e),
           cumu_symp_pe_low = cumsum(p_e_low),
           cumu_symp_pe_high = cumsum(p_e_high))
  
# how many cases would this be
expected_cases_low <- round(p_e_low * sum(res$parameters$population))
expected_cases_med <- round(p_e_med * sum(res$parameters$population))
expected_cases_high <- round(p_e_high * sum(res$parameters$population))

# how many were model observed
observed_cases_low <- round(infs_sum$cumu_symp_pe_low[which(infs_sum$date == "2020-06-02")] * sum(res$parameters$population))
observed_cases_med <- round(infs_sum$cumu_symp_pe_med[which(infs_sum$date == "2020-06-02")] * sum(res$parameters$population))
observed_cases_high <- round(infs_sum$cumu_symp_pe_high[which(infs_sum$date == "2020-06-02")] * sum(res$parameters$population))

# get likelihoods against cumulative cases described by a poisson 

ll_cases_low <- dpois(observed_cases_low, lambda = expected_cases_low, log = TRUE) 
ll_cases_med <- dpois(observed_cases_med, lambda = expected_cases_med, log = TRUE) 
ll_cases_high <- dpois(observed_cases_high, lambda = expected_cases_high, log = TRUE) 

# and alternative as just the mae
mae_cases_low <- abs(observed_cases_low - expected_cases_low)/expected_cases_low
mae_cases_med <- abs(observed_cases_med - expected_cases_med)/expected_cases_med
mae_cases_high <- abs(observed_cases_high - expected_cases_high)/expected_cases_high

## SEROLOGY / PCR relationship lls

# rolling funtion for summarising historic pcr and sero
roll_func <- function(x, det) {
  l <- length(det)
  c(NA, 
    zoo::rollapply(x, 
              list(seq(-l, -1)),
              function(i) {
                sum(i*tail(det, length(i)), na.rm = TRUE)
              },
              partial = 1
    ))
}

# susceptible population size
inf <- squire::format_output(res, c("S"), date_0 = max(res$pmcmc_results$inputs$data$date)) %>% 
  mutate(S = as.integer(y)) %>% 
  group_by(replicate) %>%  
  mutate(infections = lag(S, 1)-S) %>% 
  select(replicate, t, date, S, infections)

# correctly format
inf <- left_join(inf, 
                 squire::format_output(res, c("infections"), date_0 = max(res$pmcmc_results$inputs$data$date)) %>% 
                   mutate(symptoms = as.integer(y)) %>% 
                   select(replicate, t, date, symptoms), 
                 by = c("replicate", "t", "date"))

# calc pcr and sero
pcr_det <- readRDS("pcr_det.rds")
sero_det <- readRDS("sero_det.rds")
inf <- inf %>% 
  group_by(replicate) %>%
  mutate(pcr_positive = roll_func(infections, pcr_det),
         sero_positive = roll_func(symptoms, sero_det),
         ps_ratio = pcr_positive/sero_positive, 
         sero_perc = sero_positive/max(S,na.rm = TRUE),
         pcr_perc = pcr_positive/max(S,na.rm = TRUE)) %>% 
  ungroup %>% 
  filter(date < as.Date("2020-07-05") & date > as.Date("2020-05-22")) %>% 
  group_by(replicate) %>% 
  summarise(ps = mean(ps_ratio, na.rm = TRUE))

df <- data.frame("ll_low" = ll_cases_low, 
                 "ll_med" = ll_cases_med,
                 "ll_high" = ll_cases_high, 
                 "mae_cases_low" = mae_cases_low,
                 "mae_cases_med" = mae_cases_med,
                 "mae_cases_high" = mae_cases_high,
                 "ps_ratio" = inf$ps,
                 "ps_diff" = inf$ps - 291/155,
                 "observed_cases_low" = observed_cases_low,
                 "observed_cases_med" = observed_cases_med,
                 "observed_cases_high" = observed_cases_high,
                 "replicate" = seq_along(ll_cases_low),
                 "pop" = sum(res$parameters$population),
                 "rf" = res$pmcmc_results$inputs$pars_obs$phi_death,
                 "hosp_beds" = sum(res$parameters$hosp_bed_capacity),
                 "icu_beds" = sum(res$parameters$ICU_bed_capacity),
                 "poorer_health_outcomes" = poorer_health_outcomes,
                 "city_age" = city_age,
                 "hospital_normal_use" = hospital_normal_use
                 )

return(df)

})
saveRDS(ll, file.path(here::here(), "analysis/khartoum/ll_11_22.rds"))

## -----------------------------------------------------------------------------
## Cemetery Information
## -----------------------------------------------------------------------------

cem_df <- data.frame(
  "cemetery" = as.character(mapply(rep, c("Farouk", "Press", "Ahmed Sharafi", "Hamad al-Nile", "Bandray", "Hillat Hamad"), 4)),
  "month" = c("April","May","April","May"),
  "year" = c(2020, 2020, 2019, 2019),
  "burials" = c(87,145,25,35,224,449,45,129,54,332,90,103,500,600,120,150,206,480,192,242,302,210,25,45)
)

## -----------------------------------------------------------------------------
## Step 5: Quick forward simulation and then save
## -----------------------------------------------------------------------------
res_long <- squire::projections(res, 
                                R0_change = c(1), 
                                tt_R0 = c(0), 
                                time_period = 120)

# save space remove args here
res_long$projection_args$r <- NULL

# round to save memory. 
res_long$output <- round(res_long$output)


saveRDS(res_long, "res.rds")
saveRDS(model_fit_summary, "likelihood.rds")