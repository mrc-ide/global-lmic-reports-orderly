RhpcBLASctl::blas_set_num_threads(1L)
RhpcBLASctl::omp_set_num_threads(1L)
## -----------------------------------------------------------------------------
## 0. Checks and Function Definitions
## -----------------------------------------------------------------------------

# check on the state name
if (!(state %in% c(
  "Andaman and Nicobar Islands","Andhra Pradesh","Arunachal Pradesh","Assam",
  "Bihar","Chandigarh","Chhattisgarh","Dadra and Nagar Haveli and Daman and Diu",
  "Delhi","Goa","Gujarat","Haryana","Himachal Pradesh","Jammu and Kashmir",
  "Jharkhand","Karnataka","Kerala","Ladakh","Lakshadweep","Madhya Pradesh",
  "Maharashtra","Manipur","Meghalaya","Mizoram","Nagaland","Odisha","Puducherry",
  "Punjab","Rajasthan","Sikkim","Tamil Nadu","Telangana","Tripura","Uttar Pradesh",
  "Uttarakhand","West Bengal"))) {
  stop("State is not correct")
}

# fit spline only model function
fit_spline_rt <- function(data,
                          country,
                          pop,
                          reporting_fraction,
                          vacc_inputs,
                          n_mcmc = 10000,
                          replicates = 100,
                          rw_duration = 14,
                          hosp_beds = 10000000000,
                          icu_beds = 10000000000) {
  
  
  
  ## -----------------------------------------------------------------------------
  ## Step 1 DATA CLEANING AND ORDERING
  ## -----------------------------------------------------------------------------
  
  # order data
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
  
  ## -----------------------------------------------------------------------------
  ## Step 2a: PMCMC SETUP
  ## -----------------------------------------------------------------------------
  
  # dat_0 is just the current date now
  date_0 <- max(data$date)
  
  # what is the date of first death
  null_na <- function(x) {if(is.null(x)) {NA} else {x}}
  min_death_date <- data$date[which(data$deaths>0)][1]
  
  # We set the R0_change here to be 1 everywhere to effectively turn off mobility
  R0_change <- rep(1, nrow(data))
  date_R0_change <- data$date
  R0_change <- R0_change[as.Date(date_R0_change) <= date_0]
  date_R0_change <- date_R0_change[as.Date(date_R0_change) <= date_0]
  
  # pmcmc args
  n_particles <- 2 # we use the deterministic model now so this does nothing (makes your life quicker and easier too)
  n_chains <- 3 # number of chains
  start_adaptation <- max(2, round(n_mcmc/10)) # how long before adapting
  
  # parallel call
  suppressWarnings(future::plan(future::multiprocess()))
  
  # Defualt parameter edges for pmcmc
  R0_min <- 1.6
  R0_max <- 5.6
  last_start_date <- as.Date(null_na(min_death_date))-10
  first_start_date <- as.Date(null_na(min_death_date))-55
  start_date <- as.Date(null_na(min_death_date))-30
  
  # These 4 parameters do nothign as setting R0_change to 1 
  Meff_min <- -2
  Meff_max <- 2
  Meff_pl_min <- 0
  Meff_pl_max <- 1
  Rt_shift_min <- 0
  Rt_shift_max <- 0.001
  Rt_shift_scale_min <- 0.1
  Rt_shift_scale_max <- 10
  
  
  ## -----------------------------------------------------------------------------
  ## Step 2b: Sourcing suitable starting conditions
  ## -----------------------------------------------------------------------------
  
  date_start <- data$date[which(cumsum(data$deaths)>10)[1]] - 30
  R0_start <- 3
  
  # These are the the initial conditions now loaded from our previous run. 
  R0_start <- min(max(R0_start, R0_min), R0_max)
  date_start <- min(max(as.Date(start_date), as.Date(first_start_date)), as.Date(last_start_date))
  
  # again these all do nothing
  Meff_start <- min(max(0, Meff_min), Meff_max)
  Meff_pl_start <- min(max(0.5, Meff_pl_min), Meff_pl_max)
  Rt_shift_start <- min(max(0.0005, Rt_shift_min), Rt_shift_max)
  Rt_shift_scale_start <- min(max(5, Rt_shift_scale_min), Rt_shift_scale_max)
  
  # Our random walk parameters start after the Meff change
  # Basically just set this suitably far back in the past
  date_Meff_change <- date_start - 1
  
  ## -----------------------------------------------------------------------------
  ## Step 2c: Spline set up
  ## -----------------------------------------------------------------------------
  
  last_shift_date <- as.Date(date_Meff_change) + 1
  remaining_days <- as.Date(date_0) - last_shift_date - 14 # reporting delay in place
  
  # how many spline pars do we need
  Rt_rw_duration <- rw_duration # i.e. we fit with a 2 week duration for our random walks. 
  rw_needed <- as.numeric(ceiling(remaining_days/Rt_rw_duration))
  
  # set up rw pars
  pars_init_rw <- as.list(rep(0, rw_needed))
  pars_min_rw <- as.list(rep(-5, rw_needed))
  pars_max_rw <- as.list(rep(5, rw_needed))
  pars_discrete_rw <- as.list(rep(FALSE, rw_needed))
  names(pars_init_rw) <- names(pars_min_rw) <- names(pars_max_rw) <- names(pars_discrete_rw) <- paste0("Rt_rw_", seq_len(rw_needed))
  
  ## -----------------------------------------------------------------------------
  ## Step 2d: PMCMC parameter set up
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
      dnorm(x = pars[["Meff"]], mean = 0, sd = 1, log = TRUE) +
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
  
  # Defaults for now for vaccines 
  strategy <- "HCW, Elderly and High-Risk"
  available_doses_proportion <- 0.95  
  vaccine_uptake <- 0.8 
  vaccine_coverage_mat <- get_coverage_mat(
    iso3c = "IND",
    pop = pop,
    available_doses_proportion = available_doses_proportion, 
    strategy = strategy,
    vaccine_uptake = vaccine_uptake
  )
  
  ## -----------------------------------------------------------------------------
  ## Step 3: Run PMCMC
  ## -----------------------------------------------------------------------------
  
  # mixing matrix - assume is same as country as whole
  mix_mat <- squire::get_mixing_matrix(country)
  
  # run the pmcmc
  res <- squire::pmcmc(data = data, 
                       gibbs_days = NULL,
                       gibbs_sampling = FALSE,
                       n_mcmc = n_mcmc,
                       log_prior = logprior,
                       n_particles = 1,
                       steps_per_day = 1,
                       log_likelihood = NULL,
                       reporting_fraction = reporting_fraction,
                       squire_model = nimue::nimue_deterministic_model(),
                       output_proposals = FALSE,
                       n_chains = n_chains,
                       pars_init = pars_init,
                       pars_min = pars_min,
                       pars_max = pars_max,
                       pars_discrete = pars_discrete,
                       proposal_kernel = proposal_kernel,
                       population = pop,
                       baseline_contact_matrix = mix_mat,
                       R0_change = R0_change,
                       date_R0_change = date_R0_change,
                       Rt_args = squire:::Rt_args_list(
                         date_Meff_change = date_Meff_change,
                         scale_Meff_pl = TRUE,
                         Rt_shift_duration = 1,
                         Rt_rw_duration = Rt_rw_duration), 
                       burnin = ceiling(n_mcmc/10),
                       seeding_cases = 5,
                       replicates = replicates,
                       required_acceptance_ratio = 0.20,
                       start_adaptation = start_adaptation,
                       baseline_hosp_bed_capacity = hosp_beds, 
                       baseline_ICU_bed_capacity = icu_beds,
                       date_vaccine_change = vacc_inputs$date_vaccine_change,
                       max_vaccine = vacc_inputs$max_vaccine,
                       baseline_max_vaccine = 0,
                       date_vaccine_efficacy_infection_change = vacc_inputs$date_vaccine_change,
                       vaccine_efficacy_infection = vacc_inputs$vaccine_efficacy_infection,
                       baseline_vaccine_efficacy_infection = vacc_inputs$vaccine_efficacy_infection[[1]],
                       date_vaccine_efficacy_disease_change = vacc_inputs$date_vaccine_change,
                       vaccine_efficacy_disease = vacc_inputs$vaccine_efficacy_disease,
                       baseline_vaccine_efficacy_disease = vacc_inputs$vaccine_efficacy_disease[[1]],
                       rel_infectiousness_vaccinated = vacc_inputs$rel_infectiousness_vaccinated, 
                       vaccine_coverage_mat = vaccine_coverage_mat,
                       dur_R = 365,
                       dur_V = 5000) 
  
  
  ## remove things so they don't atke up so much memory when you save them :)
  
  # Add the prior
  res$pmcmc_results$inputs$prior <- as.function(c(formals(logprior), 
                                                  body(logprior)), 
                                                envir = new.env(parent = environment(stats::acf)))
  
  # remove states to keep object memory save down
  if("chains" %in% names(res$pmcmc_results)) {
    for(i in seq_along(res$pmcmc_results$chains)) {
      res$pmcmc_results$chains[[i]]$states <- NULL
      res$pmcmc_results$chains[[i]]$covariance_matrix <- tail(res$pmcmc_results$chains$chain1$covariance_matrix,1)
    }
  } else {
    res$pmcmc_results$states <- NULL
    res$pmcmc_results$covariance_matrix <- tail(res$pmcmc_results$covariance_matrix, 1)
  }
  
  return(res)
  
}

## -----------------------------------------------------------------------------
## 1. GET INPUT DATA
## -----------------------------------------------------------------------------

## a. Get from local files
## -----------------------------------------------------------------------------

# get pop data from file
pop <- readRDS("pop.rds")
demog <- readRDS("demog.rds")
n <- demog$n[demog$state == state]
pop <- round((sum(n)/pop$population[pop$state == state])*n)

# get icu beds from file
icu_beds <- readRDS("icu_beds.rds")
icu_beds <- icu_beds$icu_beds[icu_beds$state == state]

# get hosp beds from file
hosp_beds <- readRDS("hosp_beds.rds")
hosp_beds <- hosp_beds$hosp_beds[hosp_beds$state == state]

## b. Get from remote changing sources
## -----------------------------------------------------------------------------

# get death data
subnat_df <- read.csv("https://api.covid19india.org/csv/latest/states.csv") %>% 
  filter(!(State %in% c("India", "State Unassigned"))) %>% 
  mutate(Date = as.Date(Date)) %>% 
  group_by(State) %>% 
  complete(Date = seq.Date(min(as.Date(Date)), max(as.Date(Date)), 1)) %>% 
  mutate(State = replace_na(State, unique(!is.na(State)))) %>% 
  mutate(cases = replace_na(Confirmed, 0),
         deaths = replace_na(Deceased, 0)) %>% 
  select(-c("Tested", "Recovered", "Other", "Confirmed", "Deceased")) %>% 
  rename(date = Date, state = State)
df <- subnat_df[subnat_df$state == state, ] %>% 
  ungroup %>% 
  select(date, deaths, cases) %>% 
  arrange(date)
df$deaths <- c(df$deaths[1], diff(df$deaths))
df$cases <- c(df$cases[1], diff(df$cases))

# get vaccination data
subnat_vacc <- read.csv("https://raw.githubusercontent.com/sociepy/covid19-vaccination-subnational/main/data/countries/India.csv")
subnat_vacc <- subnat_vacc %>% 
  mutate(region = replace(region, which(region_iso %in% c("IN-DN", "IN-DD")), "Dadra and Nagar Haveli and Daman and Diu"),
         region_iso = replace(region_iso, which(region_iso %in% c("IN-DN", "IN-DD")), "IN-DD")) %>% 
  group_by(region, date, region_iso) %>% 
  summarise(total_vaccinations = sum(total_vaccinations),
            people_vaccinated = sum(people_vaccinated),
            people_fully_vaccinated = sum(people_fully_vaccinated)) %>% 
  rename(state = region)
subnat_vacc <- subnat_vacc[subnat_vacc$state == state, ]
vacc_inputs <- get_vaccine_inputs(max(df$date), subnat_vacc)

# fit model
res <- fit_spline_rt(data = df, 
                     country = as.character("India"), 
                     pop = pop, 
                     reporting_fraction = as.numeric(rf), 
                     vacc_inputs = vacc_inputs,
                     n_mcmc = as.numeric(n_mcmc),
                     replicates = as.numeric(replicates),
                     hosp_beds = as.numeric(hosp_beds),
                     icu_beds = as.numeric(icu_beds)) 

# make a quick plot so we can check fits easily
rtp <- rt_plot_immunity_vaccine(res)
dp <- plot(res, particle_fit = TRUE)
ggsave("fitting.pdf",width=12, height=6, 
       cowplot::plot_grid(rtp$plot + ggtitle(state), dp, ncol = 1))

# add state for ease and remove the output for memory
res$parameters$state <- state
res$output <- NULL

# save output
saveRDS(res, "res.rds")
