# fit spline only model function
fit_spline_rt <- function(data,
                          country,
                          pop,
                          min_rf,
                          max_rf, 
                          vacc_inputs,
                          sero_df, 
                          sero_det,
                          pars_obs_dur_R = 365,
                          pars_obs_prob_hosp_multiplier = 1,
                          model = "SQUIRE",
                          n_mcmc = 10000,
                          replicates = 100,
                          rw_duration = 14,
                          hosp_beds = 10000000000,
                          icu_beds = 10000000000
) {
  
  
  
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
  R0_min <- 1.5
  R0_max <- 10
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
                   "Rt_shift_scale" = Rt_shift_scale_start,
                   "rf" = min_rf + (1-min_rf)/2)
  pars_min = list('start_date' = first_start_date, 
                  'R0' = R0_min, 
                  'Meff' = Meff_min, 
                  'Meff_pl' = Meff_pl_min,
                  "Rt_shift" = Rt_shift_min,
                  "Rt_shift_scale" = Rt_shift_scale_min,
                  "rf" = min_rf)
  pars_max = list('start_date' = last_start_date, 
                  'R0' = R0_max, 
                  'Meff' = Meff_max, 
                  'Meff_pl' = Meff_pl_max,
                  "Rt_shift" = Rt_shift_max,
                  "Rt_shift_scale" = Rt_shift_scale_max,
                  "rf" = max_rf)
  pars_discrete = list('start_date' = TRUE, 'R0' = FALSE, 'Meff' = FALSE, 
                       'Meff_pl' = FALSE, "Rt_shift" = FALSE, "Rt_shift_scale" = FALSE,
                       "rf" = FALSE)
  pars_obs = list(phi_cases = 1, k_cases = 2, phi_death = 1, k_death = 2, exp_noise = 1e6,
                  sero_df = sero_df, sero_det = sero_det, 
                  dur_R = pars_obs_dur_R, 
                  prob_hosp_multiplier = pars_obs_prob_hosp_multiplier)
  
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
      dunif(x = pars[["R0"]], min = 1.5, max = 10, log = TRUE) +
      dnorm(x = pars[["Meff"]], mean = 0, sd = 1, log = TRUE) +
      dunif(x = pars[["Meff_pl"]], min = 0, max = 1, log = TRUE) +
      dnorm(x = pars[["Rt_shift"]], mean = 0, sd = 1, log = TRUE) +
      dunif(x = pars[["Rt_shift_scale"]], min = 0.1, max = 10, log = TRUE) + 
      dunif(x = pars[["rf"]], min = 0.01, max = 1, log = TRUE)
    
    
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
  
  # mixing matrix - assume is same as country as whole
  mix_mat <- squire::get_mixing_matrix(country)
  
  ## -----------------------------------------------------------------------------
  ## Step 3: Run PMCMC
  ## -----------------------------------------------------------------------------
  
  pi <- readRDS("pars_init.rds")
  if(pars_obs$dur_R == 365 && pars_obs$prob_hosp_multiplier == 1) {
    pi <- pi$optimistic
  } else if (pars_obs$dur_R == 180 && pars_obs$prob_hosp_multiplier == 1) {
    pi <- pi$central
  } else if (pars_obs$dur_R < 180 && pars_obs$prob_hosp_multiplier > 1) {
    pi <- pi$worst
  } 
  pf <- pi[[state]]
  pf$start_date <- as.Date(pf$start_date)
  pos_mat <- match(names(pars_init), names(pf))
  pars_init[which(!is.na(pos_mat))] <- as.list(pf[na.omit(pos_mat)])
  
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
  
  if (model == "SQUIRE") {
    squire_model = squire:::deterministic_model()
  } else if(model == "NIMUE") {
    squire_model = nimue::nimue_deterministic_model(use_dde = TRUE)
  }
  
  # run the pmcmc
  res <- pmcmc_india(data = data, 
                     gibbs_days = NULL,
                     gibbs_sampling = FALSE,
                     n_mcmc = n_mcmc,
                     log_prior = logprior,
                     n_particles = 1,
                     steps_per_day = 1,
                     log_likelihood = india_log_likelihood,
                     reporting_fraction = pars_init$rf,
                     squire_model = squire_model,
                     output_proposals = FALSE,
                     n_chains = n_chains,
                     pars_init = pars_init,
                     pars_min = pars_min,
                     pars_max = pars_max,
                     pars_discrete = pars_discrete,
                     pars_obs = pars_obs,
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
  
  # set rf to mean of sample
  res$pmcmc_results$inputs$pars_obs$phi_death <- mean(res$replicate_parameters$rf)
  
  return(res)
  
}


india_log_likelihood <- function(pars, data, squire_model, model_params, pars_obs, n_particles, 
                                 forecast_days = 0, return = "ll", Rt_args, interventions, ...) {
  switch(return, full = {
    save_particles <- TRUE
    full_output <- TRUE
    pf_return <- "sample"
  }, ll = {
    save_particles <- FALSE
    forecast_days <- 0
    full_output <- FALSE
    pf_return <- "single"
  }, {
    stop("Unknown return type to calc_loglikelihood")
  })
  squire:::assert_in(c("R0", "start_date"), names(pars), message = "Must specify R0, start date to infer")
  R0 <- pars[["R0"]]
  start_date <- pars[["start_date"]]
  pars_obs$phi_death <- pars[["rf"]]
  squire:::assert_pos(R0)
  squire:::assert_date(start_date)
  R0_change <- interventions$R0_change
  date_R0_change <- interventions$date_R0_change
  date_contact_matrix_set_change <- interventions$date_contact_matrix_set_change
  date_ICU_bed_capacity_change <- interventions$date_ICU_bed_capacity_change
  date_hosp_bed_capacity_change <- interventions$date_hosp_bed_capacity_change
  date_vaccine_change <- interventions$date_vaccine_change
  date_vaccine_efficacy_infection_change <- interventions$date_vaccine_efficacy_infection_change
  date_vaccine_efficacy_disease_change <- interventions$date_vaccine_efficacy_disease_change
  if (is.null(date_R0_change)) {
    tt_beta <- 0
  }
  else {
    tt_list <- squire:::intervention_dates_for_odin(dates = date_R0_change, 
                                                    change = R0_change, start_date = start_date, steps_per_day = round(1/model_params$dt), 
                                                    starting_change = 1)
    model_params$tt_beta <- tt_list$tt
    R0_change <- tt_list$change
    date_R0_change <- tt_list$dates
  }
  if (is.null(date_contact_matrix_set_change)) {
    tt_contact_matrix <- 0
  }
  else {
    tt_list <- squire:::intervention_dates_for_odin(dates = date_contact_matrix_set_change, 
                                                    change = seq_along(interventions$contact_matrix_set)[-1], 
                                                    start_date = start_date, steps_per_day = round(1/model_params$dt), 
                                                    starting_change = 1)
    model_params$tt_matrix <- tt_list$tt
    model_params$mix_mat_set <- model_params$mix_mat_set[tt_list$change, 
                                                         , ]
  }
  if (is.null(date_ICU_bed_capacity_change)) {
    tt_ICU_beds <- 0
  }
  else {
    tt_list <- squire:::intervention_dates_for_odin(dates = date_ICU_bed_capacity_change, 
                                                    change = interventions$ICU_bed_capacity[-1], start_date = start_date, 
                                                    steps_per_day = round(1/model_params$dt), starting_change = interventions$ICU_bed_capacity[1])
    model_params$tt_ICU_beds <- tt_list$tt
    model_params$ICU_beds <- tt_list$change
  }
  if (is.null(date_hosp_bed_capacity_change)) {
    tt_hosp_beds <- 0
  }
  else {
    tt_list <- squire:::intervention_dates_for_odin(dates = date_hosp_bed_capacity_change, 
                                                    change = interventions$hosp_bed_capacity[-1], start_date = start_date, 
                                                    steps_per_day = round(1/model_params$dt), starting_change = interventions$hosp_bed_capacity[1])
    model_params$tt_hosp_beds <- tt_list$tt
    model_params$hosp_beds <- tt_list$change
  }
  if (is.null(date_vaccine_change)) {
    tt_vaccine <- 0
  }
  else {
    tt_list <- squire:::intervention_dates_for_odin(dates = date_vaccine_change, 
                                                    change = interventions$max_vaccine[-1], start_date = start_date, 
                                                    steps_per_day = round(1/model_params$dt), starting_change = interventions$max_vaccine[1])
    model_params$tt_vaccine <- tt_list$tt
    model_params$max_vaccine <- tt_list$change
  }
  if (is.null(date_vaccine_efficacy_infection_change)) {
    tt_vaccine_efficacy_infection <- 0
  }
  else {
    tt_list <- squire:::intervention_dates_for_odin(dates = date_vaccine_efficacy_infection_change, 
                                           change = seq_along(interventions$vaccine_efficacy_infection)[-1], 
                                           start_date = start_date, steps_per_day = round(1/model_params$dt), 
                                           starting_change = 1)
    model_params$tt_vaccine_efficacy_infection <- tt_list$tt
    model_params$vaccine_efficacy_infection <- model_params$vaccine_efficacy_infection[tt_list$change, 
                                                                                       , ]
  }
  if (is.null(date_vaccine_efficacy_disease_change)) {
    tt_vaccine_efficacy_disease <- 0
  }
  else {
    tt_list <- squire:::intervention_dates_for_odin(dates = date_vaccine_efficacy_disease_change, 
                                                    change = seq_along(interventions$vaccine_efficacy_disease)[-1], 
                                                    start_date = start_date, steps_per_day = round(1/model_params$dt), 
                                                    starting_change = 1)
    model_params$tt_vaccine_efficacy_disease <- tt_list$tt
    model_params$prob_hosp <- model_params$prob_hosp[tt_list$change, 
                                                     , ]
  }
  R0 <- squire:::evaluate_Rt_pmcmc(R0_change = R0_change, R0 = R0, date_R0_change = date_R0_change, 
                                   pars = pars, Rt_args = Rt_args)
  beta_set <- squire:::beta_est(squire_model = squire_model, model_params = model_params, 
                                R0 = R0)
  model_params$beta_set <- beta_set
  if (inherits(squire_model, "stochastic")) {
    pf_result <- squire:::run_particle_filter(data = data, squire_model = squire_model, 
                                              model_params = model_params, model_start_date = start_date, 
                                              obs_params = pars_obs, n_particles = n_particles, 
                                              forecast_days = forecast_days, save_particles = save_particles, 
                                              full_output = full_output, return = pf_return)
  }
  else if (inherits(squire_model, "deterministic")) {
    pf_result <- run_deterministic_comparison_india(data = data, 
                                                    squire_model = squire_model, model_params = model_params, 
                                                    model_start_date = start_date, obs_params = pars_obs, 
                                                    forecast_days = forecast_days, save_history = save_particles, 
                                                    return = pf_return)
  }
  pf_result
  
}



run_deterministic_comparison_india <- function(data, squire_model, model_params, model_start_date = "2020-02-02",
                                               obs_params = list(
                                                 phi_cases = 0.1,
                                                 k_cases = 2,
                                                 phi_death = 1,
                                                 k_death = 2,
                                                 exp_noise = 1e+06
                                               ), forecast_days = 0, save_history = FALSE,
                                               return = "ll") {
  
  if (!(return %in% c("full", "ll", "sample", "single"))) {
    stop("return argument must be full, ll, sample", "single")
  }
  if (as.Date(data$date[data$deaths > 0][1], "%Y-%m-%d") < 
      as.Date(model_start_date, "%Y-%m-%d")) {
    stop("Model start date is later than data start date")
  }
  
  # set up as normal
  data <- squire:::particle_filter_data(data = data, start_date = model_start_date, 
                                        steps_per_day = round(1/model_params$dt))
  model_params$tt_beta <- round(model_params$tt_beta * model_params$dt)
  model_params$tt_contact_matrix <- round(model_params$tt_contact_matrix * 
                                            model_params$dt)
  model_params$tt_hosp_beds <- round(model_params$tt_hosp_beds * 
                                       model_params$dt)
  model_params$tt_ICU_beds <- round(model_params$tt_ICU_beds * 
                                      model_params$dt)
  
  # steps as normal
  steps <- c(0, data$day_end)
  fore_steps <- seq(data$day_end[nrow(data)], length.out = forecast_days + 
                      1L)
  steps <- unique(c(steps, fore_steps))
  
  if("dur_R" %in% names(obs_params)) {
    if(obs_params$dur_R != 365) {
      ch_dur_R <- as.integer(as.Date("2021-03-01") - model_start_date)
      model_params$tt_dur_R <- c(0, ch_dur_R, ch_dur_R+60)
      model_params$gamma_R <- c(model_params$gamma_R, 2/obs_params$dur_R, model_params$gamma_R)
    }
  }
  
  if("prob_hosp_multiplier" %in% names(obs_params)) {
    if(obs_params$prob_hosp_multiplier != 1) {
      ch_dur_R <- as.integer(as.Date("2021-03-01") - model_start_date)
      model_params$tt_prob_hosp_multiplier <- c(0, ch_dur_R)
      model_params$prob_hosp_multiplier <- c(model_params$prob_hosp_multiplier, obs_params$prob_hosp_multiplier)
    }
  }
  
  # run model
  model_func <- squire_model$odin_model(user = model_params, 
                                        unused_user_action = "ignore")
  out <- model_func$run(t = seq(0, tail(steps, 1), 1))
  index <- squire:::odin_index(model_func)
  
  # get deaths for comparison
  Ds <- diff(rowSums(out[, index$D]))
  Ds <- Ds[data$day_end[-1]]
  Ds[Ds < 0] <- 0
  deaths <- data$deaths[-1]
  
  # what type of ll for deaths
  if (obs_params$treated_deaths_only) {
    Ds_heathcare <- diff(rowSums(out[, index$D_get]))
    Ds_heathcare <- Ds_heathcare[data$day_end[-1]]
    ll <- squire:::ll_nbinom(deaths, Ds_heathcare, obs_params$phi_death, 
                             obs_params$k_death, obs_params$exp_noise)
  }
  else {
    ll <- squire:::ll_nbinom(deaths, Ds, obs_params$phi_death, obs_params$k_death, 
                             obs_params$exp_noise)
  }
  
  # now the ll for the seroprevalence
  sero_df <- obs_params$sero_df
  if(nrow(sero_df) > 0) {
    
    sero_at_date <- function(date, symptoms, det, dates, N) {
      
      di <- which(dates == date)
      to_sum <- tail(symptoms[seq_len(di)], length(det))
      min(sum(rev(to_sum)*head(det, length(to_sum)), na.rm=TRUE)/N, 0.99)
      
    }
    
    # get symptom incidence
    symptoms <- rowSums(out[,index$E2]) * model_params$gamma_E
    
    # dates of incidence, pop size and dates of sero surveys 
    dates <- data$date[[1]] + seq_len(nrow(out)) - 1L
    N <- sum(model_params$population)
    sero_dates <- list(sero_df$date_end, sero_df$date_start, sero_df$date_start + as.integer((sero_df$date_end - sero_df$date_start)/2))
    unq_sero_dates <- unique(c(sero_df$date_end, sero_df$date_start, sero_df$date_start + as.integer((sero_df$date_end - sero_df$date_start)/2)))
    det <- obs_params$sero_det
    
    # estimate model seroprev
    sero_model <- vapply(unq_sero_dates, sero_at_date, numeric(1), symptoms, det, dates, N)
    sero_model_mat <- do.call(cbind,lapply(sero_dates, function(x) {sero_model[match(x, unq_sero_dates)]}))
    
    # likelihood of model obvs
    lls <- rowMeans(dbinom(sero_df$sero_pos, sero_df$samples, sero_model_mat, log = TRUE))
    
  } else {
    lls <- 0
  }
  
  # and wrap up as normal
  date <- data$date[[1]] + seq_len(nrow(out)) - 1L
  rownames(out) <- as.character(date)
  attr(out, "date") <- date
  pf_results <- list()
  pf_results$log_likelihood <- sum(ll) + sum(lls)
  if (save_history) {
    pf_results$states <- out
  }
  else if (return == "single") {
    pf_results$sample_state <- out[nrow(out), ]
  }
  if (return == "ll") {
    ret <- pf_results$log_likelihood
  }
  else if (return == "sample") {
    ret <- pf_results$states
  }
  else if (return == "single" || return == "full") {
    ret <- pf_results
  }
  ret
}




pmcmc_india <- function(data,
                        n_mcmc,
                        log_likelihood = NULL,
                        log_prior = NULL,
                        n_particles = 1e2,
                        steps_per_day = 4,
                        output_proposals = FALSE,
                        n_chains = 1,
                        squire_model = explicit_model(),
                        pars_obs = list(phi_cases = 1,
                                        k_cases = 2,
                                        phi_death = 1,
                                        k_death = 2,
                                        exp_noise = 1e6),
                        pars_init = list('start_date'     = as.Date("2020-02-07"),
                                         'R0'             = 2.5,
                                         'Meff'           = 2,
                                         'Meff_pl'        = 3,
                                         "R0_pl_shift"    = 0),
                        pars_min = list('start_date'      = as.Date("2020-02-01"),
                                        'R0'              = 0,
                                        'Meff'            = 1,
                                        'Meff_pl'         = 2,
                                        "R0_pl_shift"     = -2),
                        pars_max = list('start_date'      = as.Date("2020-02-20"),
                                        'R0'              = 5,
                                        'Meff'            = 3,
                                        'Meff_pl'         = 4,
                                        "R0_pl_shift"     = 5),
                        pars_discrete = list('start_date' = TRUE,
                                             'R0'         = FALSE,
                                             'Meff'       = FALSE,
                                             'Meff_pl'    = FALSE,
                                             "R0_pl_shift" = FALSE),
                        proposal_kernel = NULL,
                        scaling_factor = 1,
                        reporting_fraction = 1,
                        treated_deaths_only = FALSE,
                        country = NULL,
                        population = NULL,
                        contact_matrix_set = NULL,
                        baseline_contact_matrix = NULL,
                        date_contact_matrix_set_change = NULL,
                        R0_change = NULL,
                        date_R0_change = NULL,
                        hosp_bed_capacity = NULL,
                        baseline_hosp_bed_capacity = NULL,
                        date_hosp_bed_capacity_change = NULL,
                        ICU_bed_capacity = NULL,
                        baseline_ICU_bed_capacity = NULL,
                        date_ICU_bed_capacity_change = NULL,
                        date_vaccine_change = NULL,
                        baseline_max_vaccine = NULL,
                        max_vaccine = NULL,
                        date_vaccine_efficacy_infection_change = NULL,
                        baseline_vaccine_efficacy_infection = NULL,
                        vaccine_efficacy_infection = NULL,
                        date_vaccine_efficacy_disease_change = NULL,
                        baseline_vaccine_efficacy_disease = NULL,
                        vaccine_efficacy_disease = NULL,
                        Rt_args = NULL,
                        burnin = 0,
                        replicates = 100,
                        forecast = 0,
                        required_acceptance_ratio = 0.23,
                        start_adaptation = round(n_mcmc/2),
                        gibbs_sampling = FALSE,
                        gibbs_days = NULL,
                        ...) {
  
  #------------------------------------------------------------
  # Section 1 of pMCMC Wrapper: Checks & Setup
  #------------------------------------------------------------
  
  #--------------------
  # assertions & checks
  #--------------------
  
  # if nimue keep to 1 step per day
  if(inherits(squire_model, "nimue_model")) {
    steps_per_day <- 1
  }
  
  # we work with pars_init being a list of inital conditions for starting
  if(any(c("start_date", "R0") %in% names(pars_init))) {
    pars_init <- list(pars_init)
  }
  
  # make it same length as chains, which allows us to pass in multiple starting points
  if(length(pars_init) != n_chains) {
    pars_init <- rep(pars_init, n_chains)
    pars_init <- pars_init[seq_len(n_chains)]
  }
  
  # data assertions
  squire:::assert_dataframe(data)
  squire:::assert_in("date", names(data))
  squire:::assert_in("deaths", names(data))
  squire:::assert_date(data$date)
  squire:::assert_increasing(as.numeric(as.Date(data$date)),
                             message = "Dates must be in increasing order")
  
  # check input pars df
  squire:::assert_list(pars_init)
  squire:::assert_list(pars_init[[1]])
  squire:::assert_list(pars_min)
  squire:::assert_list(pars_max)
  squire:::assert_list(pars_discrete)
  squire:::assert_eq(names(pars_init[[1]]), names(pars_min))
  squire:::assert_eq(names(pars_min), names(pars_max))
  squire:::assert_eq(names(pars_max), names(pars_discrete))
  squire:::assert_in(c("R0", "start_date"),names(pars_init[[1]]),
                     message = "Params to infer must include R0, start_date")
  squire:::assert_date(pars_init[[1]]$start_date)
  squire:::assert_date(pars_min$start_date)
  squire:::assert_date(pars_max$start_date)
  if (pars_max$start_date >= as.Date(data$date[1])-1) {
    stop("Maximum start date must be at least 2 days before the first date in data")
  }
  
  # check date variables are as Date class
  for(i in seq_along(pars_init)) {
    pars_init[[i]]$start_date <- as.Date(pars_init[[i]]$start_date)
  }
  pars_min$start_date <- as.Date(pars_min$start_date)
  pars_max$start_date <- as.Date(pars_max$start_date)
  
  # check bounds
  for(var in names(pars_init[[1]])) {
    
    squire:::assert_bounded(as.numeric(pars_init[[1]][[var]]),
                            left = as.numeric(pars_min[[var]]),
                            right = as.numeric(pars_max[[var]]),
                            name = paste(var, "init"))
    
    squire:::assert_single_numeric(as.numeric(pars_min[[var]]), name = paste(var, "min"))
    squire:::assert_single_numeric(as.numeric(pars_max[[var]]), name = paste(var, "max"))
    squire:::assert_single_numeric(as.numeric(pars_init[[1]][[var]]), name = paste(var, "init"))
    
  }
  
  # additonal checks that R0 is positive as undefined otherwise
  squire:::assert_pos(pars_min$R0)
  squire:::assert_pos(pars_max$R0)
  squire:::assert_pos(pars_init[[1]]$R0)
  squire:::assert_bounded(pars_init[[1]]$R0, left = pars_min$R0, right = pars_max$R0)
  
  # check proposal kernel
  squire:::assert_matrix(proposal_kernel)
  if (gibbs_sampling) {
    squire:::assert_eq(colnames(proposal_kernel), names(pars_init[[1]][-1]))
    squire:::assert_eq(rownames(proposal_kernel), names(pars_init[[1]][-1]))
  } else {
    squire:::assert_eq(colnames(proposal_kernel), names(pars_init[[1]]))
    squire:::assert_eq(rownames(proposal_kernel), names(pars_init[[1]]))
  }
  
  # check likelihood items
  if ( !(is.null(log_likelihood) | inherits(log_likelihood, "function")) ) {
    stop("Log Likelihood (log_likelihood) must be null or a user specified function")
  }
  if ( !(is.null(log_prior) | inherits(log_prior, "function")) ) {
    stop("Log Likelihood (log_likelihood) must be null or a user specified function")
  }
  squire:::assert_logical(unlist(pars_discrete))
  squire:::assert_list(pars_obs)
  squire:::assert_in(c("phi_cases", "k_cases", "phi_death", "k_death", "exp_noise"), names(pars_obs))
  squire:::assert_numeric(unlist(pars_obs[c("phi_cases", "k_cases", "phi_death", "k_death", "exp_noise")]))
  
  # mcmc items
  squire:::assert_pos_int(n_mcmc)
  squire:::assert_pos_int(n_chains)
  squire:::assert_pos_int(n_particles)
  squire:::assert_logical(output_proposals)
  
  # squire and odin
  squire:::assert_custom_class(squire_model, "squire_model")
  squire:::assert_pos_int(steps_per_day)
  squire:::assert_numeric(reporting_fraction)
  squire:::assert_bounded(reporting_fraction, 0, 1, inclusive_left = FALSE, inclusive_right = TRUE)
  squire:::assert_pos_int(replicates)
  
  # date change items
  squire:::assert_same_length(R0_change, date_R0_change)
  # checks that dates are not in the future compared to our data
  if (!is.null(date_R0_change)) {
    squire:::assert_date(date_R0_change)
    if(as.Date(tail(date_R0_change,1)) > as.Date(tail(data$date, 1))) {
      stop("Last date in date_R0_change is greater than the last date in data")
    }
  }
  
  # ------------------------------------
  # checks on odin interacting variables
  # ------------------------------------
  
  if(!is.null(contact_matrix_set)) {
    squire:::assert_list(contact_matrix_set)
  }
  squire:::assert_same_length(contact_matrix_set, date_contact_matrix_set_change)
  squire:::assert_same_length(ICU_bed_capacity, date_ICU_bed_capacity_change)
  squire:::assert_same_length(hosp_bed_capacity, date_hosp_bed_capacity_change)
  squire:::assert_same_length(max_vaccine, date_vaccine_change)
  squire:::assert_same_length(vaccine_efficacy_infection, date_vaccine_efficacy_infection_change)
  squire:::assert_same_length(vaccine_efficacy_disease, date_vaccine_efficacy_disease_change)
  
  # handle contact matrix changes
  if(!is.null(date_contact_matrix_set_change)) {
    
    squire:::assert_date(date_contact_matrix_set_change)
    squire:::assert_list(contact_matrix_set)
    
    if(is.null(baseline_contact_matrix)) {
      stop("baseline_contact_matrix can't be NULL if date_contact_matrix_set_change is provided")
    }
    if(as.Date(tail(date_contact_matrix_set_change,1)) > as.Date(tail(data$date, 1))) {
      stop("Last date in date_contact_matrix_set_change is greater than the last date in data")
    }
    
    # Get in correct format
    if(is.matrix(baseline_contact_matrix)) {
      baseline_contact_matrix <- list(baseline_contact_matrix)
    }
    
    tt_contact_matrix <- c(0, seq_len(length(date_contact_matrix_set_change)))
    contact_matrix_set <- append(baseline_contact_matrix, contact_matrix_set)
    
  } else {
    tt_contact_matrix <- 0
    contact_matrix_set <- baseline_contact_matrix
  }
  
  # handle ICU changes
  if(!is.null(date_ICU_bed_capacity_change)) {
    
    squire:::assert_date(date_ICU_bed_capacity_change)
    squire:::assert_vector(ICU_bed_capacity)
    squire:::assert_numeric(ICU_bed_capacity)
    
    if(is.null(baseline_ICU_bed_capacity)) {
      stop("baseline_ICU_bed_capacity can't be NULL if date_ICU_bed_capacity_change is provided")
    }
    squire:::assert_numeric(baseline_ICU_bed_capacity)
    if(as.Date(tail(date_ICU_bed_capacity_change,1)) > as.Date(tail(data$date, 1))) {
      stop("Last date in date_ICU_bed_capacity_change is greater than the last date in data")
    }
    
    tt_ICU_beds <- c(0, seq_len(length(date_ICU_bed_capacity_change)))
    ICU_bed_capacity <- c(baseline_ICU_bed_capacity, ICU_bed_capacity)
    
  } else {
    tt_ICU_beds <- 0
    ICU_bed_capacity <- baseline_ICU_bed_capacity
  }
  
  # handle vaccine changes
  if(!is.null(date_vaccine_change)) {
    
    squire:::assert_date(date_vaccine_change)
    squire:::assert_vector(max_vaccine)
    squire:::assert_numeric(max_vaccine)
    squire:::assert_numeric(baseline_max_vaccine)
    
    if(is.null(baseline_max_vaccine)) {
      stop("baseline_max_vaccine can't be NULL if date_vaccine_change is provided")
    }
    if(as.Date(tail(date_vaccine_change,1)) > as.Date(tail(data$date, 1))) {
      stop("Last date in date_vaccine_change is greater than the last date in data")
    }
    
    tt_vaccine <- c(0, seq_len(length(date_vaccine_change)))
    max_vaccine <- c(baseline_max_vaccine, max_vaccine)
    
  } else {
    tt_vaccine <- 0
    if(!is.null(baseline_max_vaccine)) {
      max_vaccine <- baseline_max_vaccine
    } else {
      max_vaccine <- 0
    }
  }
  
  # handle vaccine efficacy disease changes
  if(!is.null(date_vaccine_efficacy_infection_change)) {
    
    squire:::assert_date(date_vaccine_efficacy_infection_change)
    if(!is.list(vaccine_efficacy_infection)) {
      vaccine_efficacy_infection <- list(vaccine_efficacy_infection)
    }
    squire:::assert_vector(vaccine_efficacy_infection[[1]])
    squire:::assert_numeric(vaccine_efficacy_infection[[1]])
    squire:::assert_numeric(baseline_vaccine_efficacy_infection)
    
    if(is.null(baseline_vaccine_efficacy_infection)) {
      stop("baseline_vaccine_efficacy_infection can't be NULL if date_vaccine_efficacy_infection_change is provided")
    }
    if(as.Date(tail(date_vaccine_efficacy_infection_change,1)) > as.Date(tail(data$date, 1))) {
      stop("Last date in date_vaccine_efficacy_infection_change is greater than the last date in data")
    }
    
    tt_vaccine_efficacy_infection <- c(0, seq_len(length(date_vaccine_efficacy_infection_change)))
    vaccine_efficacy_infection <- c(list(baseline_vaccine_efficacy_infection), vaccine_efficacy_infection)
    
  } else {
    tt_vaccine_efficacy_infection <- 0
    if(!is.null(baseline_vaccine_efficacy_infection)) {
      vaccine_efficacy_infection <- baseline_vaccine_efficacy_infection
    } else {
      vaccine_efficacy_infection <- rep(0.8, 17)
    }
  }
  
  # handle vaccine efficacy disease changes
  if(!is.null(date_vaccine_efficacy_disease_change)) {
    
    squire:::assert_date(date_vaccine_efficacy_disease_change)
    if(!is.list(vaccine_efficacy_disease)) {
      vaccine_efficacy_disease <- list(vaccine_efficacy_disease)
    }
    squire:::assert_vector(vaccine_efficacy_disease[[1]])
    squire:::assert_numeric(vaccine_efficacy_disease[[1]])
    squire:::assert_numeric(baseline_vaccine_efficacy_disease)
    
    if(is.null(baseline_vaccine_efficacy_disease)) {
      stop("baseline_vaccine_efficacy_disease can't be NULL if date_vaccine_efficacy_disease_change is provided")
    }
    if(as.Date(tail(date_vaccine_efficacy_disease_change,1)) > as.Date(tail(data$date, 1))) {
      stop("Last date in date_vaccine_efficacy_disease_change is greater than the last date in data")
    }
    
    tt_vaccine_efficacy_disease <- c(0, seq_len(length(date_vaccine_efficacy_disease_change)))
    vaccine_efficacy_disease <- c(list(baseline_vaccine_efficacy_disease), vaccine_efficacy_disease)
    
  } else {
    tt_vaccine_efficacy_disease <- 0
    if(!is.null(baseline_vaccine_efficacy_disease)) {
      vaccine_efficacy_disease <- baseline_vaccine_efficacy_disease
    } else {
      vaccine_efficacy_disease <- rep(0.95, 17)
    }
  }
  
  
  # handle hosp bed changed
  if(!is.null(date_hosp_bed_capacity_change)) {
    
    squire:::assert_date(date_hosp_bed_capacity_change)
    squire:::assert_vector(hosp_bed_capacity)
    squire:::assert_numeric(hosp_bed_capacity)
    
    if(is.null(baseline_hosp_bed_capacity)) {
      stop("baseline_hosp_bed_capacity can't be NULL if date_hosp_bed_capacity_change is provided")
    }
    squire:::assert_numeric(baseline_hosp_bed_capacity)
    if(as.Date(tail(date_hosp_bed_capacity_change,1)) > as.Date(tail(data$date, 1))) {
      stop("Last date in date_hosp_bed_capacity_change is greater than the last date in data")
    }
    
    tt_hosp_beds <- c(0, seq_len(length(date_hosp_bed_capacity_change)))
    hosp_bed_capacity <- c(baseline_hosp_bed_capacity, hosp_bed_capacity)
    
  } else {
    tt_hosp_beds <- 0
    hosp_bed_capacity <- baseline_hosp_bed_capacity
  }
  
  #----------------
  # Generate Odin items
  #----------------
  
  # make the date definitely a date
  data$date <- as.Date(as.character(data$date))
  
  # adjust for reporting fraction
  pars_obs$phi_cases <- reporting_fraction
  pars_obs$phi_death <- reporting_fraction
  pars_obs$treated_deaths_only <- treated_deaths_only
  
  # build model parameters
  model_params <- squire_model$parameter_func(
    country = country,
    population = population,
    dt = 1/steps_per_day,
    contact_matrix_set = contact_matrix_set,
    tt_contact_matrix = tt_contact_matrix,
    hosp_bed_capacity = hosp_bed_capacity,
    tt_hosp_beds = tt_hosp_beds,
    ICU_bed_capacity = ICU_bed_capacity,
    tt_ICU_beds = tt_ICU_beds,
    max_vaccine = max_vaccine,
    tt_vaccine = tt_vaccine,
    vaccine_efficacy_infection = vaccine_efficacy_infection,
    tt_vaccine_efficacy_infection = tt_vaccine_efficacy_infection,
    vaccine_efficacy_disease = vaccine_efficacy_disease,
    tt_vaccine_efficacy_disease = tt_vaccine_efficacy_disease,
    ...)
  
  # collect interventions for odin model likelihood
  interventions <- list(
    R0_change = R0_change,
    date_R0_change = date_R0_change,
    date_contact_matrix_set_change = date_contact_matrix_set_change,
    contact_matrix_set = contact_matrix_set,
    date_ICU_bed_capacity_change = date_ICU_bed_capacity_change,
    ICU_bed_capacity = ICU_bed_capacity,
    date_hosp_bed_capacity_change = date_hosp_bed_capacity_change,
    hosp_bed_capacity = hosp_bed_capacity,
    date_vaccine_change = date_vaccine_change,
    max_vaccine = max_vaccine,
    date_vaccine_efficacy_disease_change = date_vaccine_efficacy_disease_change,
    vaccine_efficacy_disease = vaccine_efficacy_disease,
    date_vaccine_efficacy_infection_change = date_vaccine_efficacy_infection_change,
    vaccine_efficacy_infection = vaccine_efficacy_infection
  )
  
  #----------------..
  # Collect Odin and MCMC Inputs
  #----------------..
  inputs <- list(
    data = data,
    n_mcmc = n_mcmc,
    model_params = model_params,
    interventions = interventions,
    pars_obs = pars_obs,
    Rt_args = Rt_args,
    squire_model = squire_model,
    pars = list(pars_obs = pars_obs,
                pars_init = pars_init,
                pars_min = pars_min,
                pars_max = pars_max,
                proposal_kernel = proposal_kernel,
                scaling_factor = scaling_factor,
                pars_discrete = pars_discrete),
    n_particles = n_particles)
  
  
  #----------------
  # create prior and likelihood functions given the inputs
  #----------------
  
  if(is.null(log_prior)) {
    # set improper, uninformative prior
    log_prior <- function(pars) log(1e-10)
  }
  calc_lprior <- log_prior
  
  if(is.null(log_likelihood)) {
    log_likelihood <- squire:::calc_loglikelihood
  } else if (!('...' %in% names(formals(log_likelihood)))){
    stop('log_likelihood function must be able to take unnamed arguments')
  }
  
  # create shorthand function to calc_ll given main inputs
  calc_ll <- function(pars) {
    X <- log_likelihood(pars = pars,
                        data = data,
                        squire_model = squire_model,
                        model_params = model_params,
                        interventions = interventions,
                        pars_obs = pars_obs,
                        n_particles = n_particles,
                        forecast_days = 0,
                        Rt_args = Rt_args,
                        return = "ll"
    )
    X
  }
  
  #----------------
  # create mcmc run functions depending on whether Gibbs Sampling
  #----------------
  
  if(gibbs_sampling) {
    # checking gibbs days is specified and is an integer
    if (is.null(gibbs_days)) {
      stop("if gibbs_sampling == TRUE, gibbs_days must be specified")
    }
    squire:::assert_int(gibbs_days)
    
    # create our gibbs run func wrapper
    run_mcmc_func <- function(...) {
      force(gibbs_days)
      squire:::run_mcmc_chain_gibbs(..., gibbs_days = gibbs_days)
    }
  } else {
    run_mcmc_func <- squire:::run_mcmc_chain
  }
  
  #----------------
  # proposals
  #----------------
  
  # needs to be a vector to pass to reflecting boundary function
  pars_min <- unlist(pars_min)
  pars_max <- unlist(pars_max)
  pars_discrete <- unlist(pars_discrete)
  
  #--------------------------------------------------------
  # Section 2 of pMCMC Wrapper: Run pMCMC
  #--------------------------------------------------------
  
  # Run the chains in parallel
  message("Running pMCMC...")
  if (Sys.getenv("SQUIRE_PARALLEL_DEBUG") == "TRUE") {
    
    chains <- purrr::pmap(
      .l =  list(n_mcmc = rep(n_mcmc, n_chains),
                 curr_pars = pars_init),
      .f = run_mcmc_func,
      inputs = inputs,
      calc_lprior = calc_lprior,
      calc_ll = calc_ll,
      first_data_date = data$date[1],
      output_proposals = output_proposals,
      required_acceptance_ratio = required_acceptance_ratio,
      start_adaptation = start_adaptation,
      proposal_kernel = proposal_kernel,
      scaling_factor = scaling_factor,
      pars_discrete = pars_discrete,
      pars_min = pars_min,
      pars_max = pars_max)
    
  } else {
    
    chains <- furrr::future_pmap(
      .l =  list(n_mcmc = rep(n_mcmc, n_chains),
                 curr_pars = pars_init),
      .f = run_mcmc_func,
      inputs = inputs,
      calc_lprior = calc_lprior,
      calc_ll = calc_ll,
      first_data_date = data$date[1],
      output_proposals = output_proposals,
      required_acceptance_ratio = required_acceptance_ratio,
      start_adaptation = start_adaptation,
      proposal_kernel = proposal_kernel,
      scaling_factor = scaling_factor,
      pars_discrete = pars_discrete,
      pars_min = pars_min,
      pars_max = pars_max,
      .progress = TRUE,
      .options = furrr::furrr_options(seed = NULL))
    
  }
  
  #----------------
  # MCMC diagnostics and tidy
  #----------------
  if (n_chains > 1) {
    names(chains) <- paste0('chain', seq_len(n_chains))
    
    # calculating rhat
    # convert parallel chains to a coda-friendly format
    chains_coda <- lapply(chains, function(x) {
      
      traces <- x$results
      if('start_date' %in% names(pars_init[[1]])) {
        traces$start_date <- squire:::start_date_to_offset(data$date[1], traces$start_date)
      }
      
      coda::as.mcmc(traces[, names(pars_init[[1]])])
    })
    
    rhat <- tryCatch(expr = {
      x <- coda::gelman.diag(chains_coda)
      x
    }, error = function(e) {
      message('unable to calculate rhat')
    })
    
    
    pmcmc <- list(inputs = chains[[1]]$inputs,
                  rhat = rhat,
                  chains = lapply(chains, '[', -1))
    
    class(pmcmc) <- 'squire_pmcmc_list'
    
  } else {
    
    pmcmc <- chains[[1]]
    class(pmcmc) <- "squire_pmcmc"
    
  }
  #--------------------------------------------------------
  # Section 3 of pMCMC Wrapper: Sample PMCMC Results
  #--------------------------------------------------------
  pmcmc_samples <- squire:::sample_pmcmc(pmcmc_results = pmcmc,
                                         burnin = burnin,
                                         n_chains = n_chains,
                                         n_trajectories = replicates,
                                         log_likelihood = log_likelihood,
                                         n_particles = n_particles,
                                         forecast_days = forecast)
  
  #--------------------------------------------------------
  # Section 4 of pMCMC Wrapper: Tidy Output
  #--------------------------------------------------------
  
  #----------------
  # Pull Sampled results and "recreate" squire models
  #----------------
  # create a fake run object and fill in the required elements
  r <- squire_model$run_func(country = country,
                             contact_matrix_set = contact_matrix_set,
                             tt_contact_matrix = tt_contact_matrix,
                             hosp_bed_capacity = hosp_bed_capacity,
                             tt_hosp_beds = tt_hosp_beds,
                             ICU_bed_capacity = ICU_bed_capacity,
                             tt_ICU_beds = tt_ICU_beds,
                             max_vaccine = max_vaccine,
                             tt_vaccine = tt_vaccine,
                             vaccine_efficacy_infection = vaccine_efficacy_infection,
                             tt_vaccine_efficacy_infection = tt_vaccine_efficacy_infection,
                             vaccine_efficacy_disease = vaccine_efficacy_disease,
                             tt_vaccine_efficacy_disease = tt_vaccine_efficacy_disease,
                             population = population,
                             replicates = 1,
                             day_return = TRUE,
                             time_period = nrow(pmcmc_samples$trajectories),
                             ...)
  
  # and add the parameters that changed between each simulation, i.e. posterior draws
  r$replicate_parameters <- pmcmc_samples$sampled_PMCMC_Results
  
  # as well as adding the pmcmc chains so it's easy to draw from the chains again in the future
  r$pmcmc_results <- pmcmc
  
  # then let's create the output that we are going to use
  names(pmcmc_samples)[names(pmcmc_samples) == "trajectories"] <- "output"
  dimnames(pmcmc_samples$output) <- list(dimnames(pmcmc_samples$output)[[1]], dimnames(r$output)[[2]], NULL)
  r$output <- pmcmc_samples$output
  
  # and adjust the time as before
  full_row <- match(0, apply(r$output[,"time",],2,function(x) { sum(is.na(x)) }))
  saved_full <- r$output[,"time",full_row]
  for(i in seq_len(replicates)) {
    na_pos <- which(is.na(r$output[,"time",i]))
    full_to_place <- saved_full - which(rownames(r$output) == as.Date(max(data$date))) + 1L
    if(length(na_pos) > 0) {
      full_to_place[na_pos] <- NA
    }
    r$output[,"time",i] <- full_to_place
  }
  
  # second let's recreate the output
  r$model <- pmcmc_samples$inputs$squire_model$odin_model(
    user = pmcmc_samples$inputs$model_params, unused_user_action = "ignore"
  )
  # we will add the interventions here so that we know what times are needed for projection
  r$interventions <- interventions
  
  # and fix the replicates
  r$parameters$replicates <- replicates
  r$parameters$time_period <- as.numeric(diff(as.Date(range(rownames(r$output)))))
  r$parameters$dt <- model_params$dt
  
  #--------------------..
  # out
  #--------------------..
  return(r)
  
}


generate_draws <- function(out, draws = 10, parallel = TRUE, burnin = 100, log_likelihood = india_log_likelihood) {
  
  # handle for no death days
  if(!("pmcmc_results" %in% names(out))) {
    message("`out` was not generated by pmcmc as no deaths for this country. \n",
            "Returning the oroginal object, which assumes epidemic seeded on date ",
            "fits were run")
    return(out)
  }
  
  # grab information from the pmcmc run
  pmcmc <- out$pmcmc_results
  squire_model <- out$pmcmc_results$inputs$squire_model
  country <- out$parameters$country
  population <- out$parameters$population
  interventions <- out$interventions
  data <- out$pmcmc_results$inputs$data
  
  # sample parameters
  replicates <- draws
  burnin <- burnin
  if("chains" %in% names(out$pmcmc_results)) {
    n_chains <- length(out$pmcmc_results$chains)
  } else {
    n_chains <- 1
  }
  n_particles <- 2
  forecast <- 0
  
  # are we drawing in parallel
  if (parallel) {
    suppressWarnings(future::plan(future::multisession()))
  }
  
  #--------------------------------------------------------
  # Section 3 of pMCMC Wrapper: Sample PMCMC Results
  #--------------------------------------------------------
  pmcmc_samples <- squire:::sample_pmcmc(pmcmc_results = pmcmc,
                                         burnin = burnin,
                                         n_chains = n_chains,
                                         n_trajectories = replicates,
                                         n_particles = n_particles,
                                         forecast_days = forecast,
                                         log_likelihood = log_likelihood)
  
  #--------------------------------------------------------
  # Section 4 of pMCMC Wrapper: Tidy Output
  #--------------------------------------------------------
  
  # create a fake run object and fill in the required elements
  r <- squire_model$run_func(country = country,
                             contact_matrix_set = pmcmc$inputs$model_params$contact_matrix_set,
                             tt_contact_matrix = pmcmc$inputs$model_params$tt_matrix,
                             hosp_bed_capacity = pmcmc$inputs$model_params$hosp_bed_capacity,
                             tt_hosp_beds = pmcmc$inputs$model_params$tt_hosp_beds,
                             ICU_bed_capacity = pmcmc$inputs$model_params$ICU_bed_capacity,
                             tt_ICU_beds = pmcmc$inputs$model_params$tt_ICU_beds,
                             population = population,
                             day_return = TRUE,
                             replicates = 1,
                             time_period = nrow(pmcmc_samples$trajectories))
  
  # and add the parameters that changed between each simulation, i.e. posterior draws
  r$replicate_parameters <- pmcmc_samples$sampled_PMCMC_Results
  
  # as well as adding the pmcmc chains so it's easy to draw from the chains again in the future
  r$pmcmc_results <- pmcmc
  
  # then let's create the output that we are going to use
  names(pmcmc_samples)[names(pmcmc_samples) == "trajectories"] <- "output"
  dimnames(pmcmc_samples$output) <- list(dimnames(pmcmc_samples$output)[[1]], dimnames(r$output)[[2]], NULL)
  r$output <- pmcmc_samples$output
  
  # and adjust the time as before
  full_row <- match(0, apply(r$output[,"time",],2,function(x) { sum(is.na(x)) }))
  saved_full <- r$output[,"time",full_row]
  for(i in seq_len(replicates)) {
    na_pos <- which(is.na(r$output[,"time",i]))
    full_to_place <- saved_full - which(rownames(r$output) == as.Date(max(data$date))) + 1L
    if(length(na_pos) > 0) {
      full_to_place[na_pos] <- NA
    }
    r$output[,"time",i] <- full_to_place
  }
  
  # second let's recreate the output
  r$model <- pmcmc_samples$inputs$squire_model$odin_model(
    user = pmcmc_samples$inputs$model_params, unused_user_action = "ignore"
  )
  
  # we will add the interventions here so that we know what times are needed for projection
  r$interventions <- interventions
  
  # and fix the replicates
  r$parameters$replicates <- replicates
  r$parameters$time_period <- as.numeric(diff(as.Date(range(rownames(r$output)))))
  r$parameters$dt <- pmcmc$inputs$model_params$dt
  
  return(r)
  
}
