orderly_id <- tryCatch(orderly::orderly_run_info()$id,
                       error = function(e) "<id>")

print(sessionInfo())
RhpcBLASctl::blas_set_num_threads(1L)
RhpcBLASctl::omp_set_num_threads(1L)

# pandoc linking
if(file.exists("L:\\OJ\\pandoc")) {
  rmarkdown:::set_pandoc_info("L:\\OJ\\pandoc")
  Sys.setenv(RSTUDIO_PANDOC="L:\\OJ\\pandoc")
  tinytex::use_tinytex("L:\\OJ\\TinyTex")
  tinytex::tlmgr_update()
}

version_min <- "0.6.9"
if(packageVersion("squire") < version_min) {
  stop("squire needs to be updated to at least v", version_min)
}

version_min <- "0.1.8"
if(packageVersion("nimue") < version_min) {
  stop("nimue needs to be updated to at least v", version_min)
}

## -----------------------------------------------------------------------------
## Step 1: Incoming Date
## -----------------------------------------------------------------------------
system(paste0("echo Vaccine Reports for  ",iso3c, ". Short Run = ", short_run, ". Parallel = ", parallel))
set.seed(123)

# format user provided arguments correctly
date <- as.Date(date)
date_0 <- date
short_run <- as.logical(short_run)
parallel <- as.logical(parallel)
full_scenarios <- as.logical(full_scenarios)

## Get the worldometers if JHU is too erratic
ecdc <- readRDS("jhu_all.rds")
# ecdc <- readRDS("ecdc_all.rds")
if (iso3c %in% c("BOL", "ITA", "FRA", "ECU", "CHL", "COD", "ESP", "IRN",
                 "JPN", "GUF","KGZ", "PER", "HKG", "MAC", "TWN",
                 "SDN", "IRL", "TUR", "NPL")) {
  ecdc <- readRDS("worldometers_all.rds")
}

country <- squire::population$country[match(iso3c, squire::population$iso3c)[1]]
ecdc_df <- ecdc[which(ecdc$countryterritoryCode == iso3c),]
ecdc_df <- ecdc_df[ecdc_df$date <= as.Date(date_0),]

## MAIN LOOP IS ONLY FOR THOSE WITH DEATHS
if(sum(ecdc_df$deaths) > 0) {

  # Remove any deaths at beginning that were followed by 21 days of no deaths as we have no information in these situations
  if(sum(ecdc_df$deaths>0)>1) {
    if(tail(diff(which(ecdc_df$deaths>0)),1) > 21) {
      ecdc_df$deaths[tail(which(ecdc_df$deaths>0),1)] <- 0
      deaths_removed <- sum(ecdc_df$deaths[tail(which(ecdc_df$deaths>0),1)])
    } else {
      deaths_removed <- 0
    }
  } else {
    deaths_removed <- 0
  }

  # get the raw data correct
  data <- ecdc_df[,c("dateRep", "deaths", "cases")]
  names(data)[1] <- "date"
  data <- data[order(data$date),]
  data$date <- as.Date(data$date)
  data <- data[data$date <= as.Date(date_0), ]

  # Handle for countries that have eliminated and had reintroduction events
  reintroduction_iso3cs <- c("MMR", "BLZ", "TTO", "BHS", "HKG", "ABW", "GUM", "ISL", "BRB", "MUS")
  if (iso3c %in% reintroduction_iso3cs) {
    deaths_removed <- deaths_removed + sum(data$deaths[data$date < as.Date("2020-06-01")])
    data$deaths[data$date < as.Date("2020-06-01")] <- 0
  }

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
  ## Step 2: Fit model
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

  # catch for missing mobilty data or China which happened too early for our BRT to be helpful and too late in RWA/PNG case
  # as well as others that are just not at all informed by mobility :/
  spline_iso3cs <- c("CHN","MAC","TWN","KOR", "RWA", "PNG", "DZA", "COD", "SYR", "TUN", "UGA","UZB", "BEL","IRN")
  if(is.null(R0_change) || is.null(date_R0_change) || iso3c %in% spline_iso3cs) {
    date_R0_change <- seq.Date(as.Date("2020-01-01"), as.Date(date), 1)
    R0_change <- rep(1, length(date_R0_change))
  }

  # fix starting interventions to extend back for 55 days to fix rounding date issue
  date_R0_change <- c(seq.Date(as.Date(date_R0_change[1])-55, as.Date(date_R0_change[1]-1), 1), date_R0_change)
  R0_change <- c(rep(R0_change[1], 55), R0_change)

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
    n_mcmc <- 20
    n_chains <- 1
    grid_spread <- 2
    sleep <- 2
    start_adaptation <- 10
  } else {
    n_particles <- 50
    replicates <- 50
    n_mcmc <- as.integer(n_mcmc)
    n_chains <- 1
    grid_spread <- 11
    sleep <- 120
    start_adaptation <- 1000
  }

  # can't figure out why it subthreads now...
  if (parallel) {
    options("future.rng.onMisuse" = "ignore")
    suppressWarnings(future::plan(future::multisession()))
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
  pars_former <- pars_former[[iso3c]]

  if (!is.null(pars_former)) {

    R0_start <- pars_former$R0
    date_start <- pars_former$start_date
    Meff_start <- pars_former$Meff
    Meff_pl_start <- pars_former$Meff_pl
    Rt_shift_start <- pars_former$Rt_shift
    Rt_shift_duration <- pars_former$Rt_shift_duration
    Rt_shift_scale_start <- pars_former$Rt_shift_scale
    Rt_rw_duration <- pars_former$Rt_rw_duration
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
    Rt_rw_duration <- 14

    if (is.null(interventions[[iso3c]]$C)) {
      date_Meff_change <- NA
    } else {
      date_Meff_change <- post_lockdown_date_relative(interventions[[iso3c]], 1.05,
                                                      max_date = as.Date("2020-06-02"),
                                                      min_date = as.Date("2020-02-01"))
    }
  }




  R0_start <- min(max(R0_start, R0_min*1.02), R0_max*0.98)
  date_start <- min(max(as.Date(date_start), as.Date(first_start_date)+1), as.Date(last_start_date)-1)
  Meff_start <- min(max(Meff_start, Meff_min), Meff_max)
  Meff_pl_start <- min(max(Meff_pl_start, Meff_pl_min), Meff_pl_max)
  Rt_shift_start <- min(max(Rt_shift_start, Rt_shift_min), Rt_shift_max)
  Rt_shift_scale_start <- min(max(Rt_shift_scale_start, Rt_shift_scale_min), Rt_shift_scale_max)

  # either set to end if mobility dictates
  if(is.null(date_Meff_change) || is.na(date_Meff_change)) {
    date_Meff_change <- as.Date("2020-06-01")
  }

  # however if the mobility coming in is null then let's set it to the start date and that should
  # ensure that the correct number of rws are used
  if (is.null(interventions[[iso3c]]$C) || iso3c %in% spline_iso3cs) {
    date_Meff_change <- date_start
  }


  ## -----------------------------------------------------------------------------
  ## Step 2c: Spline set up
  ## -----------------------------------------------------------------------------

  last_shift_date <- as.Date(date_Meff_change) + 7
  remaining_days <- as.Date(date_0) - last_shift_date - 14 # reporting delay in place with the rounding means this is closer in practice.

  # how many spline pars do we need
  rw_needed <- as.numeric(round(remaining_days/Rt_rw_duration)) + 1 # because the first rw starts at t0

  # set up rw pars
  if (is.null(pars_former)) {
    pars_init_rw <- as.list(rep(0, rw_needed))
  } else {
    pars_init_rw <- as.list(pars_former[grep("Rt_rw_\\d",names(pars_former))])
    if(length(pars_init_rw) > rw_needed) {
      pars_init_rw <- head(pars_init_rw, rw_needed)
    }
    if(length(pars_init_rw) < rw_needed) {
      pars_init_rw <- c(pars_init_rw, vector("list", rw_needed-length(pars_init_rw)))
    }
    pars_init_rw <- lapply(pars_init_rw, function(x) {
      if(is.null(x)){
        return(0)
      } else {
        return(x)
      }})
  }

  # how long is the last spline explaining now
  last_date_rw_starts <- as.Date(last_shift_date) + Rt_rw_duration*(rw_needed-1)
  number_of_last_rw_days <- as.integer(as.Date(date_0) - last_date_rw_starts)

  pars_min_rw <- as.list(rep(-5, rw_needed))
  pars_max_rw <- as.list(rep(5, rw_needed))
  pars_discrete_rw <- as.list(rep(FALSE, rw_needed))
  names(pars_init_rw) <- names(pars_min_rw) <- names(pars_max_rw) <- names(pars_discrete_rw) <- paste0("Rt_rw_", seq_len(rw_needed))

  ## -----------------------------------------------------------------------------
  ## Step 2d: PMCMC parameter initials set up
  ## -----------------------------------------------------------------------------

  # PMCMC Parameters
  pars_init <- list(
    list('start_date' = date_start-1,
         'R0' = R0_start*0.99,
         'Meff' = Meff_start,
         'Meff_pl' = Meff_pl_start,
         "Rt_shift" = 0,
         "Rt_shift_scale" = Rt_shift_scale_start),
    list('start_date' = date_start,
         'R0' = R0_start,
         'Meff' = Meff_start,
         'Meff_pl' = Meff_pl_start,
         "Rt_shift" = 0,
         "Rt_shift_scale" = Rt_shift_scale_start),
    list('start_date' = date_start+1,
         'R0' = R0_start*1.02,
         'Meff' = Meff_start,
         'Meff_pl' = Meff_pl_start,
         "Rt_shift" = 0,
         "Rt_shift_scale" = Rt_shift_scale_start))
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

  # here use lower tolerance for countries that have had a really long time without deaths
  if(iso3c %in% c("NZL", "BRN")) {
    pars_obs$atol <- 1e-8
    pars_obs$rtol <- 1e-8
  }

  # add in the spline list
  pars_init <- lapply(pars_init, append, pars_init_rw)
  pars_min <- append(pars_min, pars_min_rw)
  pars_max <- append(pars_max, pars_max_rw)
  pars_discrete <- append(pars_discrete, pars_discrete_rw)

  # Are we doing gibbs sampling
  gibbs_sampling <- as.logical(gibbs_sampling)
  if(gibbs_sampling) {
    gibbs_days <- 1
    # Covriance Matrix
    proposal_kernel <- diag(length(names(pars_init[[1]]))-1) * 0.3
    rownames(proposal_kernel) <- colnames(proposal_kernel) <- names(pars_init[[1]][-1])
  } else {
    gibbs_sampling <- FALSE
    gibbs_days <- NULL
    # Covriance Matrix
    proposal_kernel <- diag(length(names(pars_init[[1]]))) * 0.3
    rownames(proposal_kernel) <- colnames(proposal_kernel) <- names(pars_init[[1]])
    proposal_kernel["start_date", "start_date"] <- 1.5
  }

  # use the old covar matrix and scaling factor if available
  if("covariance_matrix" %in% names(pars_former)) {

    # old proposal kernel
    proposal_kernel_proposed <- pars_former$covariance_matrix[[1]]

    # check if it needs to be expanded
    if(length(grep("Rt_rw", colnames(proposal_kernel_proposed))) == rw_needed) {

      proposal_kernel <- proposal_kernel_proposed

    } else if(length(grep("Rt_rw", colnames(proposal_kernel_proposed))) < rw_needed) {

      add_similar_cr <- function(x) {
        x <- cbind(rbind(x, 0), 0)
        rw_num <- colnames(x)[nrow(x)-1]
        new_rw <- paste0("Rt_rw_", as.numeric(gsub("(.*_)(\\d*)$", "\\2", rw_num)) + 1)
        colnames(x)[ncol(x)] <- rownames(x)[nrow(x)] <- new_rw
        x[nrow(x),] <- x[nrow(x) - 1,]
        x[,ncol(x)] <- x[,ncol(x) - 1]
        return(x)
      }

      # add as needed
      for(i in seq_len(rw_needed - length(grep("Rt_rw", colnames(proposal_kernel_proposed))))) {
        proposal_kernel_proposed <- add_similar_cr(proposal_kernel_proposed)
      }
      proposal_kernel <- proposal_kernel_proposed
    } else {

      # remove as needed
      for(i in seq_len(length(grep("Rt_rw", colnames(proposal_kernel_proposed))) - rw_needed)) {
        proposal_kernel_proposed <- proposal_kernel_proposed[-nrow(proposal_kernel_proposed),-ncol(proposal_kernel_proposed)]
      }
      proposal_kernel <- proposal_kernel_proposed

    }

    if(gibbs_sampling) {

      st_pos <- which(colnames(proposal_kernel) == "start_date")
      proposal_kernel <- proposal_kernel[-st_pos, -st_pos]

    }

  }

  scaling_factor <- 1
  if("scaling_factor" %in% names(pars_former)) {
    scaling_factor <- pars_former$scaling_factor
  }


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
  ## Step 2e: Further country specific changes
  ## -----------------------------------------------------------------------------

  # input country params
  hosp_beds <- squire:::get_hosp_bed_capacity(country)
  icu_beds <- squire:::get_ICU_bed_capacity(country)

  # Increase ICU beds where known to be too low:
  if (iso3c == "BRA") {
    # https://g1.globo.com/bemestar/coronavirus/noticia/2020/06/08/casos-de-coronavirus-e-numero-de-mortes-no-brasil-em-8-de-junho.ghtml - date we predicted ICU to be at capacity and reported to be at 70%
    icu_beds <- icu_beds / 0.7
  }

  # slight mobility data issue here
  if(iso3c == "BGD") {
    R0_change[R0_change > 1.4] <- 1
  }

  # slight hack to enforce transmission through long period with no deaths
  elong_summer_isos <- c("EST", "ISL", "ATG")
  if (iso3c %in% elong_summer_isos) {

    mmr_dates <- seq.Date(as.Date("2020-06-14"), as.Date("2020-07-14"), 21)
    old_deaths <- data$deaths[data$date %in% mmr_dates]
    data$deaths[data$date %in% mmr_dates] <- 1

  }

  elong_deaths_cont_trans <- c("VNM", "TZA", "FJI")
  if (iso3c %in% elong_deaths_cont_trans) {

    mmr_dates <- seq.Date(as.Date("2020-10-01"), as.Date("2021-05-01"), 30)
    old_deaths <- data$deaths[data$date %in% mmr_dates]
    data$deaths[data$date %in% mmr_dates] <- 1

  }

  ## -----------------------------------------------------------------------------
  ## Step 2f: Vacccine Inputs
  ## -----------------------------------------------------------------------------

  # lets read in the vaccine inputs
  vdm <- readRDS("vaccine_doses_by_manufacturer.rds")
  vacc_types <- readRDS("vaccine_agreements.rds")
  who_vacc <- readRDS("who_vacc.rds")
  who_vacc_meta <- readRDS("who_vacc_meta.rds")
  owid <- readRDS("owid.rds")
  owid <- owid %>% filter(countryterritoryCode == iso3c) %>%
    select(date, contains("vacc"))

  vacc_inputs <- get_vaccine_inputs(iso3c, vdm, vacc_types, owid, date_0, who_vacc, who_vacc_meta)

  # Defaults for now.
  strategy <- "HCW, Elderly and High-Risk"
  if(iso3c %in% get_covax_iso3c()) {
    available_doses_proportion <- 0.2
  } else {
    available_doses_proportion <- 0.95
  }
  vaccine_uptake <- 0.8
  vaccine_coverage_mat <- get_coverage_mat(
    iso3c,
    available_doses_proportion = available_doses_proportion,
    strategy = strategy,
    vaccine_uptake = vaccine_uptake
  )

  vaccine_fitting_flag <- TRUE
  squire_model <- nimue:::nimue_deterministic_model()
  init <- init_state_nimue(deaths_removed, iso3c)

  ## -----------------------------------------------------------------------------
  ## Step 2g: PMCMC run
  ## -----------------------------------------------------------------------------

  # sleep so parallel is chill
  Sys.sleep(time = runif(1, 0, sleep))

  out <- squire::pmcmc(data = data,
                       gibbs_sampling = gibbs_sampling,
                       gibbs_days = gibbs_days,
                       n_mcmc = n_mcmc,
                       log_prior = logprior,
                       n_particles = 1,
                       steps_per_day = 1,
                       log_likelihood = NULL,
                       squire_model = squire_model,
                       output_proposals = FALSE,
                       n_chains = n_chains,
                       pars_obs = pars_obs,
                       pars_init = pars_init,
                       pars_min = pars_min,
                       pars_max = pars_max,
                       pars_discrete = pars_discrete,
                       proposal_kernel = proposal_kernel,
                       scaling_factor = scaling_factor,
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
                       baseline_ICU_bed_capacity = icu_beds,
                       init = init,
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

  # Save this dummy one here for debugging purposes
  if (full_scenarios) {

    full_out <- out
    # remove states to keep object memory save down
    for(i in seq_along(full_out$pmcmc_results$chains)) {
      full_out$pmcmc_results$chains[[i]]$states <- NULL
      full_out$pmcmc_results$chains[[i]]$covariance_matrix <- tail(full_out$pmcmc_results$chains$chain1$covariance_matrix,1)
    }

    # Add the prior
    full_out$pmcmc_results$inputs$prior <- as.function(c(formals(logprior),
                                                         body(logprior)),
                                                       envir = new.env(parent = environment(stats::acf)))

    saveRDS(full_out, "pre_grad_out.rds")

  } else {
    file.create("pre_grad_out.rds")
  }


   if (vaccine_fitting_flag) {
    out <- generate_draws_pmcmc_nimue_case_fitted(out = out,
                                                  n_particles = n_particles,
                                                  grad_dur = number_of_last_rw_days)
  }  else {
    out <- generate_draws_pmcmc_case_fitted(out = out,
                                            n_particles = n_particles,
                                            grad_dur = number_of_last_rw_days)
  }


  # Add the prior
  out$pmcmc_results$inputs$prior <- as.function(c(formals(logprior),
                                                  body(logprior)),
                                                envir = new.env(parent = environment(stats::acf)))

  # slight hack to enforce transmission through long period with no deaths
  if (iso3c %in% c(elong_summer_isos, elong_deaths_cont_trans)) {
    data$deaths[data$date %in% mmr_dates] <- old_deaths
    out$pmcmc_results$inputs$data$deaths[out$pmcmc_results$inputs$data$date %in% mmr_dates] <- old_deaths
  }

  # for ones with first waves removed increase their cumulative infection counter accordingly
  if (iso3c %in% reintroduction_iso3cs) {
    index <- squire:::odin_index(out$model)

    if (vaccine_fitting_flag) {

      for(r in seq_len(replicates)) {
        first_non_na <- which(!is.na(out$output[,index$R1[1], r]))[1]
        new_m <- matrix(out$output[first_non_na,index$R1,r], nrow = nrow(out$output), ncol = length(index$R1), byrow = TRUE)
        out$output[,index$infections_cumu,r] <- out$output[,index$infections_cumu,r] + new_m
      }

    } else {

    for(r in seq_len(replicates)) {
      first_non_na <- which(!is.na(out$output[,index$R1[1], r]))[1]
      new_m <- matrix(out$output[first_non_na,index$R1,r], nrow = nrow(out$output), ncol = length(index$R1), byrow = TRUE)
      out$output[,index$cum_infs,r] <- out$output[,index$cum_infs,r] + new_m
    }

  }

  }

  ## -----------------------------------------------------------------------------
  ## Step 2d: Summarise Fits for Interface
  ## -----------------------------------------------------------------------------

  ## and save the info for the interface
  all_chains <- do.call(rbind,lapply(out$pmcmc_results$chains, "[[", "results"))
  if(is.null(all_chains)) {
    all_chains <- out$pmcmc_results$results
  }
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

  if(!is.null(out$pmcmc_results$inputs$interventions$date_R0_change)) {
    tt_beta <- squire:::intervention_dates_for_odin(dates = out$pmcmc_results$inputs$interventions$date_R0_change,
                                                    change = out$pmcmc_results$inputs$interventions$R0_change,
                                                    start_date = start_date,
                                                    steps_per_day = 1)
  } else {
    tt_beta <- 0
  }

  if(!is.null(out$pmcmc_results$inputs$interventions$R0_change)) {
    R0 <- squire:::evaluate_Rt_pmcmc(R0_change = tt_beta$change,
                                     date_R0_change = tt_beta$dates,
                                     R0 = best$R0,
                                     pars = as.list(best[1,-(1:2)]),
                                     Rt_args = out$pmcmc_results$inputs$Rt_args)
  } else {
    R0 <- R0
  }
  beta_set <- squire:::beta_est(squire_model = squire_model,
                                model_params = out$pmcmc_results$inputs$model_params,
                                R0 = R0)


  ## -----------------------------------------------------------------------------


  df <- data.frame(tt_beta = tt_beta$tt, beta_set = beta_set,
                   date = start_date + tt_beta$tt, Rt = R0,
                   grey_bar_start = FALSE)


  ## -----------------------------------------------------------------------------

  # add in uncertainty
  rts <- rt_plot_immunity_vaccine(out)

  # bind these in
  df <- dplyr::left_join(df, rts$rts[, c("date","Rt_min", "Rt_max")], by = "date")
  df <- fill(df, all_of(c("Rt_min", "Rt_max")), .direction = c("down"))
  df <- fill(df, all_of(c("Rt_min", "Rt_max")), .direction = c("up"))


  df$beta_set_min <- squire:::beta_est(squire_model = squire_model,
                                       model_params = out$pmcmc_results$inputs$model_params,
                                       R0 = df$Rt_min)

  df$beta_set_max <- squire:::beta_est(squire_model = squire_model,
                                       model_params = out$pmcmc_results$inputs$model_params,
                                       R0 = df$Rt_max)

  ## -----------------------------------------------------------------------------

  # add in grey bar start for interface
  ox_interventions <- readRDS("oxford_grt.rds")
  ox_interventions_unique <- squire:::interventions_unique(ox_interventions[[iso3c]], "C")
  df$grey_bar_start[which.min(abs(as.numeric(df$date - ox_interventions_unique$dates_change[1])))] <- TRUE

  # mark as unlikely fit because of death issues:
  df$recent_deaths <- TRUE
  if(sum(tail(out$pmcmc_results$inputs$data$deaths, 20)) == 0) {
    df$recent_deaths <- FALSE
  }

  # add in the deaths to the json fits themselves
  df$deaths <- out$pmcmc_results$inputs$data$deaths[match(df$date, out$pmcmc_results$inputs$data$date)]
  df_new_covidsim <- extend_df_for_covidsim(df = df, out = out, ext = 240)
  df_new_covidsim$iso3c <- iso3c

  # and add in the vaccine args
  df_new_covidsim <- ammend_df_covidsim_for_vaccs(df_new_covidsim, out, strategy = strategy, available_doses_proportion = available_doses_proportion)
  writeLines(jsonlite::toJSON(df_new_covidsim, pretty = TRUE), "input_params.json")


  ## -----------------------------------------------------------------------------
  ## Step 3: Summarise Fits
  ## -----------------------------------------------------------------------------

  ## summarise what we have
  g1 <- simple_pmcmc_plot(out) + theme(text = element_text(size = 2))

  title <- cowplot::ggdraw() +
    cowplot::draw_label(
      paste0(country, ",", iso3c),
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

  suppressMessages(suppressWarnings(
    cas_plot <- cases_plot_single(
      df = nimue_format(out, "infections", date_0 = date_0),
      data = out$pmcmc_results$inputs$data,
      date = date,
      date_0 = date_0) +
      theme(legend.position = "none")))


  intervention <- intervention_plot_google(
    interventions[[iso3c]], date, data,
    forecast,
    start_date = min(as.Date(out$replicate_parameters$start_date))) +
    geom_vline(xintercept = as.Date(last_shift_date) + seq(Rt_rw_duration, Rt_rw_duration*rw_needed, by = Rt_rw_duration),
               linetype = "dashed") + xlab("") + theme(axis.text.x = element_text(angle=45, vjust = 0.5))

  rtp <- rt_plot_immunity_vaccine(out, R0_plot = TRUE)
  rtp2 <- rt_plot_immunity_vaccine(out, R0_plot = FALSE)

  date_range <- as.Date(c(min(as.Date(out$replicate_parameters$start_date)),date_0))
  suppressMessages(suppressWarnings(
    bottom <- cowplot::plot_grid(
      intervention + scale_x_date(date_breaks = "1 month", date_labels = "%b" ,limits = date_range),
      d + scale_x_date(date_breaks = "1 month", date_labels = "%b" ,limits = date_range),
      cas_plot  + scale_x_date(date_breaks = "1 month", date_labels = "%b" ,limits = date_range),
      rtp$plot + scale_x_date(date_breaks = "1 month", date_labels = "%b" ,limits = date_range),
      ncol=1,
      rel_heights = c(0.4,0.6,0.6, 0.4))
  ))

  combined <- cowplot::plot_grid(header,
                                 cowplot::plot_grid(g1, line_v, bottom, ncol = 3, rel_widths = c(3,0.1,1)),
                                 ncol=1, rel_heights = c(1,15))



  ggsave("fitting.pdf",width=24, height=12,
         combined)


  ## Save the grid out object

  # remove states to keep object memory save down
  if("chains" %in% names(out$pmcmc_results)) {
  for(i in seq_along(out$pmcmc_results$chains)) {
    out$pmcmc_results$chains[[i]]$states <- NULL
    out$pmcmc_results$chains[[i]]$covariance_matrix <- tail(out$pmcmc_results$chains$chain1$covariance_matrix,1)
  }
  } else {
    out$pmcmc_results$states <- NULL
    out$pmcmc_results$covariance_matrix <- tail(out$pmcmc_results$covariance_matrix, 1)
  }

  ## -----------------------------------------------------------------------------
  ## Step 4: Scenarios
  ## -----------------------------------------------------------------------------

  ## -----------------------------------------------------------------------------
  ## 4.1. Conduct our projections based on no surging
  ## -----------------------------------------------------------------------------

  ## Functions for working out the relative changes in R0 for given scenarios
  time_period <- 365

  ## We need to know work out vaccine doses and efficacy going forwards
  model_user_args <- extend_vaccine_inputs(vacc_inputs, time_period, out)

  # Maintaining the current set of measures for a further 3 months
  maintain_scenario <- squire::projections(out,
                                           R0_change = c(1),
                                           tt_R0 = c(0),
                                           time_period = time_period,
                                           model_user_args = model_user_args)
  r_maintain_scenario <- r_list_format(maintain_scenario, date_0)
  rt_maintain_scenario <- rt_creation_vaccine(maintain_scenario, date_0, date_0+89)

  ## check for capacity being passed here

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

  # has surging been required
  icu <- nim_sq_format(maintain_scenario, "ICU_demand", date_0 = date_0)
  icu <- icu[icu$compartment == "ICU_demand",]
  icu_28 <- group_by(icu, t) %>%
    summarise(i_tot = mean(y, na.rm = TRUE),
              i_min = t_test_safe(y)$conf.int[1],
              i_max = t_test_safe(y)$conf.int[2])

  hosp <- nim_sq_format(maintain_scenario, "hospital_demand", date_0 = date_0)
  hosp <- hosp[hosp$compartment == "hospital_demand",]
  hosp_28 <- group_by(hosp, t) %>%
    summarise(i_tot = mean(y, na.rm = TRUE),
              i_min = t_test_safe(y)$conf.int[1],
              i_max = t_test_safe(y)$conf.int[2])

  if(any(na.omit(icu_28$i_tot > icu_cap)) || any(na.omit(icu_28$i_tot > icu_cap))) {

    surging <- TRUE

  } else {

    surging <- FALSE

  }

  rm(maintain_scenario)

  # Enhancing movement restrictions for 3 months (50% further reduction in contacts) (Mitigation)
  mitigation_scenario <- squire::projections(out,
                                             R0_change = c(0.5),
                                             tt_R0 = c(0),
                                             time_period = time_period,
                                             model_user_args = model_user_args)
  r_mitigation_scenario <- r_list_format(mitigation_scenario, date_0)
  rt_mitigation_scenario <- rt_creation_vaccine(mitigation_scenario, date_0, date_0+89)
  rm(mitigation_scenario)

  # Relax by 50% for 3 months
  reverse_scenario <- squire::projections(out,
                                          R0_change = 1.5,
                                          tt_R0 = c(0),
                                          time_period = time_period,
                                          model_user_args = model_user_args)
  r_reverse_scenario <- r_list_format(reverse_scenario, date_0)
  rt_reverse_scenario <- rt_creation_vaccine(reverse_scenario, date_0, date_0+89)
  rm(reverse_scenario)

  ## now let's trim the out for saving really small
  pmcmc <- out$pmcmc_results
  out_interventions <- out$interventions
  out$output <- NULL
  saveRDS(out, "grid_out.rds")
  rm(out)

  ## -----------------------------------------------------------------------------
  ## 4.2. Investigating a capacity surge
  ## -----------------------------------------------------------------------------

  # update for surged healthcare
  pmcmc$inputs$model_params$hosp_beds <- 1e10
  pmcmc$inputs$model_params$ICU_beds <- 1e10

  out_surged <- generate_draws_pmcmc(pmcmc = pmcmc,
                                     burnin = ceiling(n_mcmc/10),
                                     n_chains = n_chains,
                                     squire_model = pmcmc$inputs$squire_model,
                                     replicates = replicates,
                                     n_particles = n_particles,
                                     forecast = 0,
                                     country = country,
                                     population = squire::get_population(iso3c = iso3c)$n,
                                     interventions = out_interventions,
                                     data = pmcmc$inputs$data)

  out_surged$model$set_user(ICU_beds = 1e10)
  out_surged$model$set_user(hosp_beds = 1e10)

  out_surged$parameters$hosp_bed_capacity <- 1e10
  out_surged$parameters$ICU_bed_capacity <- 1e10


  ## Now simulate the surged scenarios
  out_surged$pmcmc_results$chains <- NULL

  maintain_scenario_surged <- squire::projections(out_surged,
                                                  R0_change = c(1),
                                                  tt_R0 = c(0),
                                                  time_period = time_period,
                                                  model_user_args = model_user_args)

  r_maintain_scenario_surged <- r_list_format(maintain_scenario_surged, date_0)
  rt_maintain_scenario_surged <- rt_creation_vaccine(maintain_scenario_surged, date_0, date_0+89)
  rm(maintain_scenario_surged)

  mitigation_scenario_surged <- squire::projections(out_surged,
                                                    R0_change = c(0.5),
                                                    tt_R0 = c(0),
                                                    time_period = time_period,
                                                    model_user_args = model_user_args)
  r_mitigation_scenario_surged <- r_list_format(mitigation_scenario_surged, date_0)
  rt_mitigation_scenario_surged <- rt_creation_vaccine(mitigation_scenario_surged, date_0, date_0+89)
  rm(mitigation_scenario_surged)

  reverse_scenario_surged <- squire::projections(out_surged,
                                                 R0_change = c(1.5),
                                                 tt_R0 = c(0),
                                                 time_period = time_period,
                                                 model_user_args = model_user_args)
  r_reverse_scenario_surged <- r_list_format(reverse_scenario_surged, date_0)
  rt_reverse_scenario_surged <- rt_creation_vaccine(reverse_scenario_surged, date_0, date_0+89)
  rm(reverse_scenario_surged)
  rm(out_surged)

  o_list <- named_list(
        r_maintain_scenario,
        r_mitigation_scenario,
        r_reverse_scenario,
        r_maintain_scenario_surged,
        r_mitigation_scenario_surged,
        r_reverse_scenario_surged
      )

  rt_list <- named_list(
    rt_maintain_scenario,
    rt_mitigation_scenario,
    rt_reverse_scenario,
    rt_maintain_scenario_surged,
    rt_mitigation_scenario_surged,
    rt_reverse_scenario_surged
  )


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
  data$date <- as.Date(data$date)
  data <- data[order(data$date, decreasing = FALSE),]
  data$deaths <- cumsum(data$deaths)
  data$cases <- cumsum(data$cases)
  data <- data[order(data$date, decreasing = TRUE),]
  data <- data[data$date <= as.Date(date_0), ]

  # prepare reports
  options(tinytex.verbose = TRUE)
  rmarkdown::render("index.Rmd",
                    output_format = c("html_document","pdf_document"),
                    params = list("o_list" = o_list,
                                  "replicates" = replicates,
                                  "data" = data,
                                  "date_0" = date_0,
                                  "country" = country,
                                  "surging" = surging,
                                  "rt" = rtp2),
                    output_options = list(pandoc_args = c(paste0("--metadata=title:",country," COVID-19 report "))))


}

## THIS IS THE ESFT LOOP FOR COUNTRIES WITH NO DEATHS CURRENTLY
if (sum(ecdc_df$deaths) == 0) {

  # What are the Rt values for each income group
  inc_R0s <- income_R0()
  inc_Rts <- income_Rt(date_0 = date)

  # what income group is this country
  wb_metadata <- read.csv("gdp_income_group.csv",
                          fileEncoding="UTF-8-BOM",
                          stringsAsFactors = TRUE)

  income <- wb_metadata$income_group[match(iso3c, wb_metadata$country_code)]

  # And the upper and lower Rs for our scenarios
  R0 <- inc_R0s$R0[inc_R0s$income == income]
  Rt <- max(inc_Rts$Rt[inc_Rts$income == income], 1.2)

  # sim_args
  time_period_esft <- 366
  replicates_esft <- 1
  seeding_cases_esft <- 5

  # lets read in the vaccine inputs
  vdm <- readRDS("vaccine_doses_by_manufacturer.rds")
  vacc_types <- readRDS("vaccine_agreements.rds")
  owid <- readRDS("owid.rds")
  owid <- owid %>% filter(countryterritoryCode == iso3c) %>%
    select(date, contains("vacc"))
  who_vacc <- readRDS("who_vacc.rds")
  who_vacc_meta <- readRDS("who_vacc_meta.rds")

  vacc_inputs <- get_vaccine_inputs(iso3c, vdm, vacc_types, owid, date_0, who_vacc, who_vacc_meta)

  # Defaults for now.
  strategy <- "HCW, Elderly and High-Risk"
  if(iso3c %in% get_covax_iso3c()) {
    available_doses_proportion <- 0.2
  } else {
    available_doses_proportion <- 0.98
  }
  vaccine_uptake <- 0.8
  vaccine_coverage_mat <- get_coverage_mat(
    iso3c,
    available_doses_proportion = available_doses_proportion,
    strategy = strategy,
    vaccine_uptake = vaccine_uptake
  )

  # set up vaccine inits
  vaccine_fitting_flag <- TRUE
  squire_model <- nimue:::nimue_deterministic_model()
  init <- init_state_nimue(deaths_removed = 0, iso3c,
                           seeding_cases = seeding_cases_esft,
                           vaccinated_already = sum(vacc_inputs$max_vaccine))

  # Scenarios with capacity constraints
  # ---------------------------------------------------------------------------
  reverse_scenario <- nimue::run(
    country = country,
    R0 = R0,
    time_period = time_period_esft,
    seeding_cases = seeding_cases_esft,
    init = init,
    max_vaccine = as.integer(mean(tail(vacc_inputs$max_vaccine,7))),
    vaccine_coverage_mat = vaccine_coverage_mat
  )

  mitigation_scenario <- nimue::run(
    country = country,
    R0 = Rt,
    time_period = time_period_esft,
    seeding_cases = seeding_cases_esft,
    init = init,
    max_vaccine = as.integer(mean(tail(vacc_inputs$max_vaccine,7))),
    vaccine_coverage_mat = vaccine_coverage_mat
  )

  maintain_scenario <- nimue::run(
    country = country,
    R0 = Rt + ((R0 - Rt)/2),
    time_period = time_period_esft,
    seeding_cases = seeding_cases_esft,
    init = init,
    max_vaccine = as.integer(mean(tail(vacc_inputs$max_vaccine,7))),
    vaccine_coverage_mat = vaccine_coverage_mat
  )
  saveRDS(maintain_scenario, "grid_out.rds")

  # Scenarios without capacity constraints
  # ---------------------------------------------------------------------------
  reverse_scenario_surged <- nimue::run(
    country = country,
    R0 = R0,
    time_period = time_period_esft,
    seeding_cases = seeding_cases_esft,
    init = init,
    max_vaccine = as.integer(mean(tail(vacc_inputs$max_vaccine,7))),
    vaccine_coverage_mat = vaccine_coverage_mat,
    hosp_bed_capacity = 1e10,
    ICU_bed_capacity = 1e10
  )

  mitigation_scenario_surged <- nimue::run(
    country = country,
    R0 = R0,
    time_period = time_period_esft,
    seeding_cases = seeding_cases_esft,
    init = init,
    max_vaccine = as.integer(mean(tail(vacc_inputs$max_vaccine,7))),
    vaccine_coverage_mat = vaccine_coverage_mat,
    hosp_bed_capacity = 1e10,
    ICU_bed_capacity = 1e10
  )

  maintain_scenario_surged <- nimue::run(
    country = country,
    R0 = R0,
    time_period = time_period_esft,
    seeding_cases = seeding_cases_esft,
    init = init,
    max_vaccine = as.integer(mean(tail(vacc_inputs$max_vaccine,7))),
    vaccine_coverage_mat = vaccine_coverage_mat,
    hosp_bed_capacity = 1e10,
    ICU_bed_capacity = 1e10
  )

  # bundle into a list
  r_list <-
    named_list(
      maintain_scenario,
      mitigation_scenario,
      reverse_scenario,
      maintain_scenario_surged,
      mitigation_scenario_surged,
      reverse_scenario_surged
    )

  # And adjust their time variable so that we have t = 0 as today
  r_list_pass <- lapply(r_list, function(x) {
    for(i in seq_len(dim(x$output)[3])) {
      x$output[,"time",i] <- x$output[,"time",i] - 1
    }
    rownames(x$output) <- as.character(seq.Date(date_0, date_0 + nrow(x$output) -1, 1))
    return(x)
  })

  ## Lastly make up some outputs here to pass orderly
  file.create("index.html", "index.pdf", "index.md",
              "summary_df.rds", "input_params.json",
              "fitting.pdf", "pre_grad_out.rds")

  # major summaries
o_list <- lapply(r_list_pass, r_list_format, date_0)

rt_list <- lapply(r_list_pass, rt_creation_vaccine, date_0, date_0+89)

}

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
    dplyr::filter(.data$date < (as.Date(date_0)+90)) %>%
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

rt_list[[1]]$scenario <- "Maintain Status Quo"
rt_list[[2]]$scenario <- "Additional 50% Reduction"
rt_list[[3]]$scenario <- "Relax Interventions 50%"
rt_list[[4]]$scenario <- "Surged Maintain Status Quo"
rt_list[[5]]$scenario <- "Surged Additional 50% Reduction"
rt_list[[6]]$scenario <- "Surged Relax Interventions 50%"

# combine and annotate
data_sum <- do.call(rbind, data_sum)
data_sum <- rbind(data_sum, do.call(rbind, rt_list)) %>% arrange(date, scenario)
rownames(data_sum) <- NULL

# catch for hong kong and taiwan country name
if(iso3c == "HKG") {
  country <- "Hong Kong"
}
if(iso3c == "TWN") {
  country <- "Taiwan"
}
if(iso3c == "MAC") {
  country <- "Macao"
}

# bring it all together
data_sum$country <- country
data_sum$iso3c <- iso3c
data_sum$report_date <- date
data_sum <- data_sum[data_sum$compartment != "D",]
data_sum$version <- "v8"
data_sum <- dplyr::mutate(data_sum, across(dplyr::starts_with("y_"), ~round(.x,digits = 2)))

# specify if this is calibrated to deaths or just hypothetical forecast for ESFT
if (sum(ecdc_df$deaths) == 0) {
  data_sum$death_calibrated <- FALSE
} else {
  data_sum$death_calibrated <- TRUE
}

write.csv(data_sum, "projections.csv", row.names = FALSE, quote = FALSE)

## -----------------------------------------------------------------------------
## Attack Rates
## -----------------------------------------------------------------------------

# ar_list <- lapply(r_list_pass, function(x) {
#
#   S <- x %>% nim_sq_format(var_select = "infections", date_0 = date) %>%
#     group_by(replicate) %>%
#     mutate(y = cumsum(y)) %>%
#     group_by(t, date) %>%
#     summarise(y = median(y,na.rm = TRUE))
#   S$ar <- S$y/sum(squire::get_population(x$parameters$country)$n)
#   S$iso <- squire::get_population(x$parameters$country)$iso3c[1]
#
#   S <- na.omit(S)
#   return(S)
#
# })
#
# ar_list[[1]]$scenario <- "Maintain Status Quo"
# ar_list[[2]]$scenario <- "Additional 50% Reduction"
# ar_list[[3]]$scenario <- "Relax Interventions 50%"
# ar_list[[4]]$scenario <- "Surged Maintain Status Quo"
# ar_list[[5]]$scenario <- "Surged Additional 50% Reduction"
# ar_list[[6]]$scenario <- "Surged Relax Interventions 50%"
#
# ars <- do.call(rbind, ar_list)
# ars$continent <- countrycode::countrycode(ars$iso, "iso3c", "continent")
# ars$region <- countrycode::countrycode(ars$iso, "iso3c", "region23")
# names(ars)[3] <- "uninfected"
# saveRDS(ars, "attack_rates.rds")
#
