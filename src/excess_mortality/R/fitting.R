# fit spline only model function
#'
#'@export
fit_spline_rt <- function(data,
                          country,
                          pop,
                          delta_characteristics,
                          vaccine_inputs,
                          model = "NIMUE",
                          n_mcmc = 10000,
                          n_chains = 3,
                          replicates = 20,
                          rw_duration = 14,
                          likelihood_version = "Negative Binomial"
) {

  #get iso3c
  iso3c <- squire::population$iso3c[match(country, squire::population$country)]

  ## -----------------------------------------------------------------------------
  ## Step 1 DATA CLEANING AND ORDERING
  ## -----------------------------------------------------------------------------

  # order data
  data <- data[order(data$week_start),]
  data$week_start <- as.Date(data$week_start)
  data$week_end <- as.Date(data$week_end)

  # and remove the rows with no data up to the first date that a death was reported
  first_report <- which(data$deaths>0)[1]
  missing <- which(data$deaths == 0 | is.na(data$deaths))
  to_remove <- missing[missing<first_report]
  if(length(to_remove) > 0) {
    if(length(to_remove) == (nrow(data)-1)) {
      data <- data[-utils::head(to_remove,-1),]
    } else {
      data <- data[-to_remove,]
    }
  }

  ## -----------------------------------------------------------------------------
  ## Step 2a: PMCMC BOUNDS SETUP
  ## -----------------------------------------------------------------------------

  # dat_0 is just the current date now
  date_0 <- max(data$week_end)

  # what is the date of first death
  null_na <- function(x) {if(is.null(x)) {NA} else {x}}
  min_death_date <- data$week_start[which(data$deaths>0)][1]

  # pmcmc args
  n_particles <- 2 # we use the deterministic model now so this does nothing (makes your life quicker and easier too)
  start_adaptation <- max(2, round(n_mcmc/10)) # how long before adapting

  # parallel call
  suppressWarnings(future::plan(future::multiprocess()))

  # Defualt parameter edges for pmcmc
  R0_min <- 1
  R0_max <- 10
  last_start_date <- as.Date(null_na(min_death_date))-10
  first_start_date <- as.Date(null_na(min_death_date))-55
  start_date <- as.Date(null_na(min_death_date))-30

  ## -----------------------------------------------------------------------------
  ## Step 2b: Sourcing suitable starting conditions
  ## -----------------------------------------------------------------------------

  date_start <- data$week_start[which(cumsum(data$deaths)>10)[1]] - 30
  R0_start <- 3

  # These are the the initial conditions now loaded from our previous run.
  R0_start <- min(max(R0_start, R0_min), R0_max)
  date_start <- min(max(as.Date(start_date), as.Date(first_start_date)), as.Date(last_start_date))

  #new variable to say when spline starts
  date_spline_start <- first_start_date + rw_duration

  ## -----------------------------------------------------------------------------
  ## Step 2c: Spline set up
  ## -----------------------------------------------------------------------------

  #calculate how many times our Rt spline changes and at what date
  date_Rt_change <- seq(date_spline_start, as.Date(date_0) - rw_duration, by = rw_duration)

  # how many spline pars do we need
  Rt_rw_duration <- rw_duration # i.e. we fit with a 2 week duration for our random walks.
  rw_needed <- length(date_Rt_change)

  # set up rw pars
  pars_init_rw <- as.list(rep(0, rw_needed))
  pars_min_rw <- as.list(rep(-5, rw_needed))
  pars_max_rw <- as.list(rep(5, rw_needed))
  pars_discrete_rw <- as.list(rep(FALSE, rw_needed))
  names(pars_init_rw) <- names(pars_min_rw) <- names(pars_max_rw) <- names(pars_discrete_rw) <- paste0("Rt_rw_", seq_len(rw_needed))

  ## -----------------------------------------------------------------------------
  ## Step 2d: PMCMC initial parameter set up
  ## -----------------------------------------------------------------------------

  # PMCMC Parameters
  pars_init = list('start_date' = date_start,
                   'R0' = R0_start)
  pars_min = list('start_date' = first_start_date,
                  'R0' = R0_min)
  pars_max = list('start_date' = last_start_date,
                  'R0' = R0_max)
  pars_discrete = list('start_date' = TRUE, 'R0' = FALSE)
  pars_obs = list(phi_cases = 1, k_cases = 2, phi_death = 1, k_death = 7, exp_noise = 1e07)
  #assign this way so they keep NULL if NULL
  pars_obs$dur_R <- delta_characteristics$required_dur_R
  pars_obs$prob_hosp_multiplier <- delta_characteristics$prob_hosp_multiplier
  pars_obs$delta_start_date <- delta_characteristics$start_date
  pars_obs$shift_duration <- delta_characteristics$shift_duration

  #set up likelihood function
  if(likelihood_version == "Negative Binomial"){
    pars_obs$likelihood <- function(model_deaths, data_deaths){
      squire:::ll_nbinom(data_deaths, model_deaths, pars_obs$phi_death,
                         pars_obs$k_death,
                         pars_obs$exp_noise)
    }
  } else if(likelihood_version == "Poisson"){
    pars_obs$likelihood <- function(model_deaths, data_deaths){
      stats::dpois(data_deaths, pars_obs$phi_death*model_deaths,
                   #+ rexp(length(model_deaths), rate = pars_obs$exp_noise),
                   log = TRUE)
    }
  } else {
    stop("likelihood_version, must be one of Poisson & Negative Binomial")
  }

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
    ret <- stats::dunif(x = pars[["start_date"]], min = -55, max = -10, log = TRUE) +
      stats::dunif(x = pars[["R0"]], min = 1, max = 10, log = TRUE)

    # get rw spline parameters
    if(any(grepl("Rt_rw", names(pars)))) {
      Rt_rws <- pars[grepl("Rt_rw", names(pars))]
      for (i in seq_along(Rt_rws)) {
        ret <- ret + stats::dnorm(x = Rt_rws[[i]], mean = 0, sd = 0.2, log = TRUE)
      }
    }
    return(ret)
  }

  if(model == "NIMUE"){
    #get the pre-calculate efficacies and durations/delays
    vacc_inputs <- vaccine_inputs
    dur_V <- vaccine_inputs$dur_V
    dur_vaccine_delay <- vaccine_inputs$dur_vaccine_delay
    vaccine_inputs$dur_V <- NULL
    vaccine_inputs$dur_vaccine_delay <- NULL
  } else {
    vacc_inputs <- NULL
    dur_V <- NULL
    dur_vaccine_delay <- NULL
  }

  # mixing matrix - assume is same as country as whole
  mix_mat <- squire::get_mixing_matrix(country)

  # Now overwrite these with the initial conditions previously found
  pf <- readRDS("pars_init.rds")
  if(!is.null(pf)){
    if(!is.null(pf[[iso3c]])){
      pf <- pf[[iso3c]]
      if(pf$start_date < pars_min$start_date | pf$start_date > pars_max$start_date){
        #if the start date is in compatible then we remove it, these can change as we are using estimates
        pf$start_date <- NULL
      }
    }
  }
  #only use if compatible and present
  pos_mat <- match(names(pars_init), names(pf))
  pars_init[which(!is.na(pos_mat))] <- as.list(pf[stats::na.omit(pos_mat)])
  pars_init$start_date <- as.Date(pars_init$start_date)

  #temp
  #pars_init$R0 <- 3
  #pars_init[grepl("Rt_rw", names(pars_init))] <- 0
  #pars_init[grepl("Rt_rw", names(pars_init))][[1]] <- c(1.609)

  # grab old scaling factor
  scaling_factor <- 1
  if("scaling_factor" %in% names(pf)) {
    scaling_factor <- as.numeric(pf$scaling_factor)
  }

  # grab old covariance matrix
  # use the old covar matrix if available
  if("covariance_matrix" %in% names(pf)) {

    # old proposal kernel
    proposal_kernel_proposed <- pf$covariance_matrix[[1]]

    #reduce to required variables
    proposal_kernel_proposed <- proposal_kernel_proposed[
      rownames(proposal_kernel_proposed) %in% rownames(proposal_kernel),
      colnames(proposal_kernel_proposed) %in% colnames(proposal_kernel)
    ]

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

  }

  ## -----------------------------------------------------------------------------
  ## Step 3: Run PMCMC
  ## -----------------------------------------------------------------------------

  # this is here if the model struggles to fit initially with nimue then we can use squire first
  # to get a better fit firt which can then be used to update the pars_init
  if (model == "SQUIRE") {
    squire_model = squire:::deterministic_model()
  } else if(model == "NIMUE") {
    squire_model = nimue::nimue_deterministic_model(use_dde = TRUE)
  }

  # run the pmcmc
  res <- pmcmc_excess(country = country,
                      data = data,
                      gibbs_days = NULL,
                      gibbs_sampling = FALSE,
                      n_mcmc = n_mcmc,
                      log_prior = logprior,
                      n_particles = 1,
                      steps_per_day = 1,
                      log_likelihood = excess_log_likelihood,
                      squire_model = squire_model,
                      output_proposals = FALSE,
                      n_chains = n_chains,
                      pars_init = pars_init,
                      pars_min = pars_min,
                      pars_max = pars_max,
                      pars_discrete = pars_discrete,
                      pars_obs = pars_obs,
                      proposal_kernel = proposal_kernel,
                      date_Rt_change = date_Rt_change,
                      Rt_args = list(
                        Rt_date_spline_start = date_spline_start,
                        Rt_rw_duration = Rt_rw_duration,
                        date_Rt_change = date_Rt_change),
                      burnin = ceiling(n_mcmc/10),
                      seeding_cases = 5,
                      replicates = replicates,
                      required_acceptance_ratio = 0.20,
                      start_adaptation = start_adaptation,
                      scaling_factor = scaling_factor,
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
                      vaccine_coverage_mat = vacc_inputs$vaccine_coverage_mat,
                      dur_R = 365,
                      dur_V = dur_V,
                      dur_vaccine_delay = dur_vaccine_delay)

  #set the class to be excess_nimue_simulation and inherit from nimue_simulation
  class(res) <- c("excess_nimue_simulation", "nimue_simulation")

  ## remove things so they don't atke up so much memory when you save them :)

  # Add the prior
  res$pmcmc_results$inputs$prior <- as.function(c(formals(logprior),
                                                  body(logprior)),
                                                envir = new.env(parent = environment(stats::acf)))

  res$interventions$delta_adjustments <- as.list(delta_characteristics)
  res$interventions$vaccine_strategy <- vaccine_inputs$strategy

  # remove states to keep object memory save down
  if("chains" %in% names(res$pmcmc_results)) {
    for(i in seq_along(res$pmcmc_results$chains)) {
      res$pmcmc_results$chains[[i]]$states <- NULL
      res$pmcmc_results$chains[[i]]$covariance_matrix <- utils::tail(res$pmcmc_results$chains$chain1$covariance_matrix,1)
    }
  } else {
    res$pmcmc_results$states <- NULL
    res$pmcmc_results$covariance_matrix <- utils::tail(res$pmcmc_results$covariance_matrix, 1)
  }

  return(res)

}
