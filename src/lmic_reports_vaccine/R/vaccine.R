get_vaccine_inputs <- function(iso3c, vdm, vacc_types, owid, date_0) {
  
  # filter to just relevant dates
  owid <- owid %>% filter(date <= as.Date(date_0))
  
  # function to interpolate missing vaccine dates
  interp_diffs <- function(date_vacc, tot) {
    
    # if there are breaks then we need to interpolate what has happened for ease
    tot <- approx(as.Date(date_vacc), tot, xout = as.Date(date_vacc))$y
    
    # we should also linearly extend back to 0 vaccinations
    early_vaccs <- predict(
      lm(y~x, data = data.frame(y = head(na.omit(tot), 7), x = 1:7)), 
      newdata = data.frame(x = -30:0), 
      type = "response")
    
    early_vaccs <- round(early_vaccs[early_vaccs > 0])
    tot[tail(which(is.na(tot)), length(early_vaccs))] <- early_vaccs
    
    # create our vacc inputs
    date_vaccine_change <- date_vacc[!is.na(tot)]
    max_vaccines <- tot[!is.na(tot)]
    max_vaccines <- as.integer(c(max_vaccines[1], diff(max_vaccines)))
    
    date_vaccine_change <- date_vaccine_change[max_vaccines > 0]
    max_vaccines <- max_vaccines[max_vaccines > 0]
    
    return(list(date = date_vaccine_change, max = max_vaccines))
    
  }
  
  # get the dates of vaccines being documented
  if(all(is.na(owid$total_vaccinations))) {
    
    
    vaccine_efficacy_infection <- rep(0.6, length(owid$date))
    vaccine_efficacy_disease <- rep(0.8, length(owid$date))
    vaccine_efficacy_infection <- lapply(vaccine_efficacy_infection, rep, 17)
    vaccine_efficacy_disease <- lapply(vaccine_efficacy_disease, rep, 17)
    
    return(list(
      date_vaccine_change = owid$date,
      max_vaccine = rep(0, length(owid$date)),
      vaccine_efficacy_infection = vaccine_efficacy_infection,
      vaccine_efficacy_disease = vaccine_efficacy_disease
    ))
  } else {
    
  tots <- interp_diffs(date_vacc = owid$date, tot = owid$total_vaccinations)
  
  # total vaccinations given out per day
  date_vaccine_change <- tots$date
  max_vaccines <- tots$max
  
  # now for doses
  seconds <- interp_diffs(date_vacc = owid$date, tot = owid$people_fully_vaccinated)
  firsts <- max_vaccines
  firsts[which(date_vaccine_change %in% seconds$date)] <- firsts[which(date_vaccine_change %in% seconds$date)] - seconds$max
  seconds <- max_vaccines - firsts
  firsts <- cumsum(firsts)
  seconds <- cumsum(seconds)
  dose_ratio <- seconds/firsts
  
  # now to work out the efficacy
  vaccine_efficacy_infection <- (1-dose_ratio)*0.6 + dose_ratio*0.8
  vaccine_efficacy_disease <- (1-dose_ratio)*0.8 + dose_ratio*0.98
  vaccine_efficacy_infection <- lapply(vaccine_efficacy_infection, rep, 17)
  vaccine_efficacy_disease <- lapply(vaccine_efficacy_disease, rep, 17)
  
  return(
    list(
      date_vaccine_change = date_vaccine_change,
      max_vaccine = max_vaccines,
      vaccine_efficacy_infection = vaccine_efficacy_infection, 
      vaccine_efficacy_disease = vaccine_efficacy_disease
    )
  )
  }
}

extend_vaccine_inputs <- function(vaccine_inputs, time_period, out) {
  
  # weekly mean vaccine distributions
  max_vaccine <- mean(tail(vacc_inputs$max_vaccine,7))
  tt_vaccine <- 0
  
  # efficacies best to just extend at the same rate
  vei <- vapply(seq_along(vacc_inputs$vaccine_efficacy_infection), 
                function(x) {
                  vacc_inputs$vaccine_efficacy_infection[[x]][1]
                }, numeric(1))
  vei_new <- predict(
    lm(y~x, data.frame("x" = seq_along(vei), "y" = vei)), 
    newdata = data.frame("x" = length(vei)+seq_len(time_period))
  )
  vei_new <- vapply(vei_new, min, numeric(1), 0.8)
  
  vaccine_efficacy_infection <- lapply(vei_new, rep, 17)
  tt_vaccine_efficacy_infection <- seq_along(vaccine_efficacy_infection)-1
  
  vaccine_efficacy_infection_odin_array <- nimue:::format_ve_i_for_odin(
    vaccine_efficacy_infection = vaccine_efficacy_infection,
    tt_vaccine_efficacy_infection = tt_vaccine_efficacy_infection
  )
  
  # efficacies best to just extend at the same rate
  ved <- vapply(seq_along(vacc_inputs$vaccine_efficacy_disease), 
                function(x) {
                  vacc_inputs$vaccine_efficacy_disease[[x]][1]
                }, numeric(1))
  ved_new <- predict(
    lm(y~x, data.frame("x" = seq_along(ved), "y" = ved)), 
    newdata = data.frame("x" = length(ved)+seq_len(time_period))
  )
  ved_new <- vapply(ved_new, min, numeric(1), 0.98)
  vaccine_efficacy_disease <- lapply(ved_new, rep, 17)
  tt_vaccine_efficacy_disease <- seq_along(vaccine_efficacy_disease)-1
  
  vaccine_efficacy_disease_odin_array <- nimue:::format_ve_d_for_odin(
    vaccine_efficacy_disease = vaccine_efficacy_disease,
    tt_vaccine_efficacy_disease = tt_vaccine_efficacy_disease,
    prob_hosp = out$parameters$prob_hosp
  )
  
  # combine into model_user_args for projections
  mua <- list(
    "max_vaccine" = max_vaccine, 
    "tt_vaccine" = 0,
    "vaccine_efficacy_infection" = vaccine_efficacy_infection_odin_array,
    "tt_vaccine_efficacy_infection" = tt_vaccine_efficacy_infection,
    "prob_hosp" = vaccine_efficacy_disease_odin_array,
    "tt_vaccine_efficacy_disease" = tt_vaccine_efficacy_disease
  )
  
  mua <- rep(list(mua), dim(out$output)[3])
  return(mua)  
  
}

get_coverage_mat <- function(iso3c, 
                             available_doses_proportion = 0.98, 
                             strategy = "HCW, Elderly and High-Risk",
                             vaccine_uptake = 0.8) {
  
  strategies <- readRDS("coverage_strategies.rds")
  
  # get cov_mat for strategy
  if(strategy == "HCW and Elderly") {
    cov_mat <- strategies[[iso3c]]$whoPriority * vaccine_uptake
  } else if (strategy == "HCW, Elderly and High-Risk") {
    cov_mat <- strategies[[iso3c]]$etagePriority * vaccine_uptake
  } else if (strategy == "Elderly") {
    cov_mat <- nimue::strategy_matrix("Elderly", max_coverage = vaccine_uptake, 0)
  } else if (strategy == "All") {
    cov_mat <- nimue::strategy_matrix("All", max_coverage = vaccine_uptake, 0)
  } else {
    stop('Incorrect strategy. Must be one of "HCW and Elderly", "HCW, Elderly and High-Risk", "Elderly", "All"')
  }
  
  # scale vaccine coverage for availability function
  scale_cov_mat <- function(cov_mat, vaccine_available, pop) {
    
    # total vaccs available
    tot_vaccines <- sum(pop*vaccine_available)
    
    # step 1, find when max allocation exceeds capacity
    step <- 1
    step_found <- FALSE
    tot_vaccs_steps <- 0
    cov_mat_dup_ex <- rbind(0, cov_mat)
    
    while(!step_found && step <= nrow(cov_mat)) {
      
      if(nrow(cov_mat) == 1) {
        step_found <- TRUE
      }
      
      vaccs_in_step <- sum((cov_mat_dup_ex[step+1, ] - cov_mat_dup_ex[step, ]) * pop)
      tot_vaccs_steps <- tot_vaccs_steps + vaccs_in_step
      if(tot_vaccs_steps > tot_vaccines) {
        step_found <- TRUE
      } else {
        step <- step+1
      }
    }
    
    # if we have enough vaccine return now
    if(step > nrow(cov_mat)) {
      return(cov_mat)
    }
    
    # set steps after max available reached to 0
    if(step < nrow(cov_mat)) {
      cov_mat[(step+1):nrow(cov_mat),] <- 0
    }
    
    # now set this step to be correct for available
    tots_given <- sum(cov_mat[step-1,] %*% pop)
    tots_tried <- sum(cov_mat[step,] %*% pop)
    remaining <- tot_vaccines - tots_given
    
    # next_group
    next_group <- cov_mat[step,]-cov_mat[step-1,]
    new_cov <- remaining/(next_group[which(next_group > 0)] * pop[which(next_group > 0)])
    cov_mat[step, which(next_group > 0)] <- new_cov
    return(cov_mat)
  }
  
  pop <- squire::get_population(iso3c = iso3c)$n
  
  cov_mat <- scale_cov_mat(cov_mat, available_doses_proportion, pop)
  return(cov_mat)
  
}


nimue_format <- function(out, 
                         var_select = NULL, 
                         reduce_age = TRUE,
                         combine_compartments = TRUE, 
                         date_0 = NULL) {
  
  
  # work out what compartments are being plotted
  compartments = c("S", "E",
                   "IMild", "ICase", "IICU", "IHospital",
                   "IRec", "R", "D")
  summaries = c("N",
                "hospitalisations",
                "hospital_demand","hospital_occupancy",
                "ICU_demand", "ICU_occupancy",
                "vaccines", "unvaccinated", "vaccinated", "priorvaccinated",
                "infections", "deaths")
  
  comps <- var_select[var_select %in% compartments]
  summs <- var_select[var_select %in% summaries]
  
  # to match with squire definition
  if("infections" %in% summs) {
    comps <- c(comps, "E2")
    summs <- summs[-which(summs == "infections")]
    inf_fix <- TRUE
  } else {
    inf_fix <- FALSE
  }

  pd <- do.call(rbind, lapply(seq_len(dim(out$output)[3]), function(i) {
    nimue::format(out, compartments = comps, summaries = summs, replicate = i)
  })) %>%
    dplyr::rename(y = .data$value)
  
  pd <- pd[,c("replicate", "compartment", "t", "y")]
  
  # replacing time with date if date_0 is provided
  if(!is.null(date_0)){
    pd$date <- as.Date(pd$t + as.Date(date_0),
                       format = "%Y-%m-%d")
  }
  
  # fix the infection 
  if (inf_fix) {
    pd$y[pd$compartment == "E2"] <- pd$y[pd$compartment == "E2"]*out$odin_parameters$gamma_E
    pd$compartment <- as.character(pd$compartment)
    pd$compartment[pd$compartment == "E2"] <- "infections"
    pd$compartment <- as.factor(pd$compartment)
  }
  
  return(pd)
  
}

nim_sq_format <- function(out, 
                          var_select = NULL, 
                          reduce_age = TRUE,
                          combine_compartments = TRUE, 
                          date_0 = NULL) {
  
  if("tt_vaccine" %in% out$model$.__enclos_env__$private$user) {
    nimue_format(out, var_select, reduce_age, combine_compartments, date_0)
  } else {
    squire::format_output(out, var_select, reduce_age, combine_compartments, date_0)
  }
  
}

init_state_nimue <- function(deaths_removed, iso3c, seeding_cases = 5) {
  
  # get an initial
  pop <- squire::get_population(iso3c = iso3c, simple_SEIR = FALSE)
  init <- nimue:::init(pop$n, seeding_cases = seeding_cases)
  
  if(deaths_removed > 0) {
    # work out how many deaths and where
    probs <- (squire:::probs$prob_hosp * squire:::probs$prob_severe * squire:::probs$prob_severe_death_treatment) +
      (squire:::probs$prob_hosp * (1-squire:::probs$prob_severe * squire:::probs$prob_non_severe_death_treatment))
    probs <- probs*pop$n
    probs <- probs/sum(probs)
    deaths <- as.numeric(t(rmultinom(1, deaths_removed, probs)))
    # approximate IFR for income group
    wb_metadata <- read.csv("gdp_income_group.csv", fileEncoding="UTF-8-BOM", stringsAsFactors = TRUE)
    income <- wb_metadata$income_group[match(iso3c, wb_metadata$country_code)]
    ifrs <- data.frame("income" = c("Low income", "Lower middle income", "Upper middle income", "High income"),
                       "ifr" = c(0.17, 0.31, 0.51, 1.02))
    ifr <- ifrs$ifr[ifrs$income == income]
    R <- rpois(1, deaths_removed*1/ifr/0.01)
    R <- as.numeric(t(rmultinom(1, R, rep(1/length(probs), length(probs)))))
    R <- R - deaths
    # and update the inital to reflect
    init$D[,1] <- deaths
    init$S[,1] <- init$S[,1] - R - deaths
    init$R1_0[,1] <- R
  }
  return(init)
}


generate_draws_pmcmc_nimue_case_fitted <- function(out, n_particles = 10, grad_dur = 21) {
  
  pmcmc <- out$pmcmc_results
  n_chains <- max(length(out$pmcmc_results$chains), 1)
  burnin <- round(out$pmcmc_results$inputs$n_mcmc/10)
  squire_model <- out$pmcmc_results$inputs$squire_model
  replicates <- dim(out$output)[3]
  forecast <- 0
  country <- out$parameters$country
  population <- out$parameters$population
  interventions <- out$interventions
  data <- out$pmcmc_results$inputs$data
  rw_dur <- out$pmcmc_results$inputs$Rt_args$Rt_rw_duration
  
  #--------------------------------------------------------
  # Section 1 # what is our predicted gradient
  #--------------------------------------------------------
  
  # first what is the model predicted infections
  infections <- nimue_format(out, "infections", date_0 = max(data$date))
  infections_end <- infections %>% filter(date > (max(data$date) - grad_dur) & date <= (max(data$date))) %>% 
    group_by(date) %>% summarise(y = median(y))
  
  infections_pre_end <- infections %>% 
    filter(date > (max(data$date) - grad_dur - rw_dur) & date <= (max(data$date) - grad_dur) ) %>% 
    group_by(date) %>% summarise(y = median(y))
  
  # and the observed cases
  cases_end <- tail(data$cases, grad_dur)
  cases_pre_end <- head(tail(data$cases, grad_dur+rw_dur), rw_dur)
  
  # get these gradients
  get_infs <- function(x) {
    sum(x, na.rm = TRUE)
  }
  
  pred_infs_end <- get_infs(infections_end$y)
  pred_infs_pre_end <- get_infs(infections_pre_end$y)
  
  des_infs_end <- get_infs(cases_end)
  des_infs_pre_end <- get_infs(cases_pre_end)
  
  # if there are less than 100 cases in both windowns then don't bother
  if(des_infs_end > 100 && des_infs_pre_end > 100) {
    
    ca_infs_frac <-  pred_infs_pre_end / des_infs_pre_end
    
    # desired model predictd final infs
    wanted_infs <- des_infs_end * ca_infs_frac
    
    # if actual infs available
    if(!is.nan(wanted_infs) || !is.na(wanted_infs) || !is.infinite(wanted_infs)) {

      # do we need to go up or down
      if(wanted_infs < pred_infs_end) {
        alters <- seq(0.025, 0.175, 0.025)
      } else {
        alters <- seq(-0.025, -0.125, -0.025) # more conservative on the way up
      }
      
      # store our grads
      ans <- alters
      last_rw <- ncol(out$pmcmc_results$chains$chain1$results) - 3
      
      #--------------------------------------------------------
      # Section 2 # # find best grad correction
      #--------------------------------------------------------
      
      for(alt in seq_along(alters)) {
        
        message(alt)  
        
        for(ch in seq_along(out$pmcmc_results$chains)) {
          out$pmcmc_results$chains[[ch]]$results[,last_rw] <- out$pmcmc_results$chains[[ch]]$results[,last_rw] + alters[alt]
        }
        
        pmcmc_samples <- squire:::sample_pmcmc(pmcmc_results = out$pmcmc_results,
                                               burnin = burnin,
                                               n_chains = n_chains,
                                               n_trajectories = replicates,
                                               n_particles = n_particles,
                                               forecast_days = forecast)
        
        dimnms <- dimnames(pmcmc_samples$trajectories)
        
        # then let's create the output that we are going to use
        names(pmcmc_samples)[names(pmcmc_samples) == "trajectories"] <- "output"
        dimnames(pmcmc_samples$output) <- list(dimnames(pmcmc_samples$output)[[1]], dimnames(out$output)[[2]], NULL)
        out$output <- pmcmc_samples$output
        
        # and adjust the time as before
        full_row <- match(0, apply(out$output[,"time",],2,function(x) { sum(is.na(x)) }))
        saved_full <- out$output[,"time",full_row]
        for(i in seq_len(replicates)) {
          na_pos <- which(is.na(out$output[,"time",i]))
          full_to_place <- saved_full - which(rownames(out$output) == as.Date(max(data$date))) + 1L
          if(length(na_pos) > 0) {
            full_to_place[na_pos] <- NA
          }
          out$output[,"time",i] <- full_to_place
        }
        
        infections <- nimue_format(out, "infections", date_0 = max(data$date))
        this_infs <- infections %>% filter(date > (max(data$date) - grad_dur) & date <= (max(data$date))) %>% 
          group_by(date) %>% summarise(y = median(y))
        
        ans[alt] <- get_infs(this_infs$y)
        
        
        for(ch in seq_along(out$pmcmc_results$chains)) {
          out$pmcmc_results$chains[[ch]]$results[,last_rw] <- out$pmcmc_results$chains[[ch]]$results[,last_rw] - alters[alt]
        }
        
      }
      
      
      # adapt our whole last chain accordingly
      alts <- which.min(abs(ans-wanted_infs))
      for(ch in seq_along(out$pmcmc_results$chains)) {
        out$pmcmc_results$chains[[ch]]$results[,last_rw] <- out$pmcmc_results$chains[[ch]]$results[,last_rw] + alters[alts]
      }
      
    }
    
  }
  
  #--------------------------------------------------------
  # Section 3 of pMCMC Wrapper: Sample PMCMC Results
  #--------------------------------------------------------
  pmcmc_samples <- squire:::sample_pmcmc(pmcmc_results = out$pmcmc_results,
                                         burnin = burnin,
                                         n_chains = n_chains,
                                         n_trajectories = replicates,
                                         n_particles = n_particles,
                                         forecast_days = forecast)
  
  #--------------------------------------------------------
  # Section 4 of pMCMC Wrapper: Tidy Output
  #--------------------------------------------------------
  dimnms <- dimnames(pmcmc_samples$trajectories)
  
  # then let's create the output that we are going to use
  names(pmcmc_samples)[names(pmcmc_samples) == "trajectories"] <- "output"
  dimnames(pmcmc_samples$output) <- list(dimnames(pmcmc_samples$output)[[1]], dimnames(out$output)[[2]], NULL)
  out$output <- pmcmc_samples$output
  
  # and adjust the time as before
  full_row <- match(0, apply(out$output[,"time",],2,function(x) { sum(is.na(x)) }))
  saved_full <- out$output[,"time",full_row]
  for(i in seq_len(replicates)) {
    na_pos <- which(is.na(out$output[,"time",i]))
    full_to_place <- saved_full - which(rownames(out$output) == as.Date(max(data$date))) + 1L
    if(length(na_pos) > 0) {
      full_to_place[na_pos] <- NA
    }
    out$output[,"time",i] <- full_to_place
  }
  
  return(out)
  
}


ammend_df_covidsim_for_vaccs <- function(df, out, strategy) {

  # add in vaccine pars
  df$max_vaccine <- 0
  df$vaccine_efficacy_infection <- out$interventions$vaccine_efficacy_infection[[1]][1]
  df$vaccine_efficacy_disease <- out$interventions$vaccine_efficacy_disease[[1]][1]
  
  vacc_pos <- which(as.Date(df$date) %in% as.Date(out$interventions$date_vaccine_change))
  df$max_vaccine[vacc_pos] <- as.integer(out$interventions$max_vaccine[-1][seq_along(vacc_pos)])
  df$max_vaccine[(max(vacc_pos)+1):length(df$max_vaccine)] <- as.integer(mean(tail(df$max_vaccine[vacc_pos], 7)))
  
  df$vaccine_efficacy_infection[vacc_pos] <- vapply(
    seq_along(out$interventions$vaccine_efficacy_infection), 
    function(x){ out$interventions$vaccine_efficacy_infection[[x]][1] }, 
    numeric(1)
  )[-1]
  
  df$vaccine_efficacy_disease[vacc_pos] <- vapply(
    seq_along(out$interventions$vaccine_efficacy_disease), 
    function(x){ out$interventions$vaccine_efficacy_disease[[x]][1] }, 
    numeric(1)
  )[-1]
  
  df$vaccine_strategy <- strategy
  df$vaccine_coverage <- max(out$pmcmc_results$inputs$model_params$vaccine_coverage_mat)
  total_vacc <- sum((tail(out$pmcmc_results$inputs$model_params$vaccine_coverage_mat,1) * out$pmcmc_results$inputs$model_params$population))
  df$vaccines_available <- total_vacc / sum(out$pmcmc_results$inputs$model_params$population)
  
  return(df)
  
}

nim_sq_simulation_plot_prep <- function(x,
                                        var_select,
                                        q = c(0.025, 0.975),
                                        summary_f = mean,
                                        x_var = "t",
                                        ...) {
  
  pd <- nim_sq_format(x, var_select = var_select, ...)
  
  pd <- pd %>%
    dplyr::mutate(x = .data[[x_var]])
  
  # t sometimes seems to be being rounded weirdly
  if(x_var == "t") {
    pd$x <- round(pd$x, ceiling(log10(1/x$parameters$dt)))
  }
  
  # remove any NA rows (due to different start dates)
  if(sum(is.na(pd$t) | is.na(pd$y))>0) {
    pd <- pd[-which(is.na(pd$t) | is.na(pd$y)),]
  }
  
  # Format summary data
  pds <- pd %>%
    dplyr::group_by(.data$x, .data$compartment) %>%
    dplyr::summarise(ymin = stats::quantile(.data$y, q[1]),
                     ymax = stats::quantile(.data$y, q[2]),
                     y = summary_f(.data$y))
  
  return(list(pd = pd, pds = pds))
  
}  


get_reff <- function(out, beta) {
  
  # mixing_matrix is already the mixing matrix that we pass to you in the country json files
  mixing_matrix <- squire:::process_contact_matrix_scaled_age(
    out$parameters$contact_matrix_set[[1]],
    out$parameters$population
  )
  
  # these parameters are found in pars_0.json that is imported in index.js
  dur_ICase <- out$parameters$dur_ICase
  dur_IMild <- out$parameters$dur_IMild
  prob_hosp <- out$odin_parameters$prob_hosp
  
  #
  index <- nimue:::odin_index(out$model)
  
  # pop is a 17 length with population sizes in each age category
  pop <- out$parameters$population
  
  # in here we work out each time point the number of individuals in each age category in
  # the S compartment at each time point.
  susceptible <- array(
    out$output[,index$S,],
    dim=c(nrow(out$output), dim(index$S))
  )
  # We divide by the total population
  prop_susc <- sweep(susceptible, 2, pop, FUN='/')
  # We multiply by the effect of vaccines on onward infectiousness
  prop_susc <- sweep(
    prop_susc,
    c(2, 3),
    out$odin_parameters$vaccine_efficacy_infection,
    FUN='*'
  )
  
  # Length 17 with relative R0 in each age category
  relative_R0_by_age <- prob_hosp*dur_ICase + (1-prob_hosp)*dur_IMild
  
  # here we are looping over each time point to calculate the adjusted eigen
  # incorporating the proportion of the susceptible population in each age group
  adjusted_eigens <- vapply(
    seq(nrow(out$output)),
    function(t) {
      Re(eigen(mixing_matrix * rowSums(prop_susc[t,,] * relative_R0_by_age))$values[1])
    },
    numeric(1)
  )
  
  # multiply beta by the adjusted eigen at each time point to get Reff
  beta * adjusted_eigens
}

get_immunity_ratios_vaccine <- function(out, max_date = NULL) {
  
  mixing_matrix <- squire:::process_contact_matrix_scaled_age(
    out$pmcmc_results$inputs$model_params$contact_matrix_set[[1]],
    out$pmcmc_results$inputs$model_params$population
  )
  
  dur_ICase <- out$parameters$dur_ICase
  dur_IMild <- out$parameters$dur_IMild
  prob_hosp <- out$parameters$prob_hosp
  
  # assertions
  squire:::assert_single_pos(dur_ICase, zero_allowed = FALSE)
  squire:::assert_single_pos(dur_IMild, zero_allowed = FALSE)
  squire:::assert_numeric(prob_hosp)
  squire:::assert_numeric(mixing_matrix)
  squire:::assert_square_matrix(mixing_matrix)
  squire:::assert_same_length(mixing_matrix[,1], prob_hosp)
  
  if(sum(is.na(prob_hosp)) > 0) {
    stop("prob_hosp must not contain NAs")
  }
  
  if(sum(is.na(mixing_matrix)) > 0) {
    stop("mixing_matrix must not contain NAs")
  }
  
  index <- squire:::odin_index(out$model)
  pop <- out$parameters$population
  
  if(is.null(max_date)) {
    max_date <- max(out$pmcmc_results$inputs$data$date)
  }
  t_now <- which(as.Date(rownames(out$output)) == max_date)
  
  # make the vaccine time changing args
  nrs <- t_now
  vei <- lapply(seq_len(nrow(out$odin_parameters$vaccine_efficacy_infection)),
                function(x) {out$odin_parameters$vaccine_efficacy_infection[x,,]}
  )
  vei_full <- c(rep(list(vei[[1]]),nrs - length(vei)), vei)
  
  
  # prop susceptible 
  prop_susc <- lapply(seq_len(dim(out$output)[3]), function(x) {
    
    susceptible <- array(
      out$output[,index$S,x],
      dim=c(t_now, dim(index$S))
    )
    
    # We divide by the total population
    prop_susc <- sweep(susceptible, 2, pop, FUN='/')
    
    # We multiply by the effect of vaccines on onward infectiousness
    prop_susc <- vapply(
      seq_len(nrow(prop_susc)),
                        FUN = function(x){prop_susc[x,,]*vei_full[[x]]},
      FUN.VALUE = prop_susc[x,,]
    )
    
    return(prop_susc)
  } )
  
  relative_R0_by_age <- prob_hosp*dur_ICase + (1-prob_hosp)*dur_IMild
  
  adjusted_eigens <- lapply(prop_susc, function(x) {
    
    unlist(lapply(seq_len(nrow(x)), function(y) {
      if(any(is.na(x[y,]))) {
        return(NA)
      } else {
        Re(eigen(mixing_matrix*x[y,]*relative_R0_by_age)$values[1])
      }
    }))
    
  })
  
  betas <- lapply(out$replicate_parameters$R0, function(x) {
    squire:::beta_est(squire_model = out$pmcmc_results$inputs$squire_model, 
                      model_params = out$pmcmc_results$inputs$model_params, 
                      R0 = x)
  })
  
  ratios <- lapply(seq_along(betas), function(x) {
    (betas[[x]] * adjusted_eigens[[x]]) / out$replicate_parameters$R0[[x]]
  })
  
  return(ratios)
}



rt_plot_immunity_vaccine <- function(out) {
  
  iso3c <- squire::get_population(out$parameters$country)$iso3c[1]
  
  if("pmcmc_results" %in% names(out)) {
    wh <- "pmcmc_results"
  } else {
    wh <- "scan_results"
  }
  
  date <- max(as.Date(out$pmcmc_results$inputs$data$date))
  date_0 <- date
  
  # impact of immunity ratios
  ratios <- get_immunity_ratios_vaccine(out)
  
  # create the Rt data frame
  rts <- lapply(seq_len(length(out$replicate_parameters$R0)), function(y) {
    
    tt <- squire:::intervention_dates_for_odin(dates = out$interventions$date_R0_change, 
                                               change = out$interventions$R0_change, 
                                               start_date = out$replicate_parameters$start_date[y],
                                               steps_per_day = 1/out$parameters$dt)
    
    if(wh == "scan_results") {
      Rt <- c(out$replicate_parameters$R0[y], 
              vapply(tt$change, out[[wh]]$inputs$Rt_func, numeric(1), 
                     R0 = out$replicate_parameters$R0[y], Meff = out$replicate_parameters$Meff[y])) 
    } else {
      Rt <- squire:::evaluate_Rt_pmcmc(
        R0_change = tt$change, 
        date_R0_change = tt$dates, 
        R0 = out$replicate_parameters$R0[y], 
        pars = as.list(out$replicate_parameters[y,]),
        Rt_args = out$pmcmc_results$inputs$Rt_args) 
    }
    
    df <- data.frame(
      "Rt" = Rt,
      "Reff" = Rt*tail(na.omit(ratios[[y]]),length(Rt)),
      "R0" = na.omit(Rt)[1]*tail(na.omit(ratios[[y]]),length(Rt)),
      "date" = tt$dates,
      "iso" = iso3c,
      rep = y,
      stringsAsFactors = FALSE)
    df$pos <- seq_len(nrow(df))
    return(df)
  } )
  
  rt <- do.call(rbind, rts)
  rt$date <- as.Date(rt$date)
  
  rt <- rt[,c(5,4,1,2,3,6,7)]
  
  new_rt_all <- rt %>%
    group_by(iso, rep) %>% 
    arrange(date) %>% 
    complete(date = seq.Date(min(rt$date), date_0, by = "days")) 
  
  column_names <- colnames(new_rt_all)[-c(1,2,3)]
  new_rt_all <- fill(new_rt_all, all_of(column_names), .direction = c("down"))
  new_rt_all <- fill(new_rt_all, all_of(column_names), .direction = c("up"))
  
  suppressMessages(sum_rt <- group_by(new_rt_all, iso, date) %>% 
                     summarise(Rt_min = quantile(Rt, 0.025),
                               Rt_q25 = quantile(Rt, 0.25),
                               Rt_q75 = quantile(Rt, 0.75),
                               Rt_max = quantile(Rt, 0.975),
                               Rt_median = median(Rt),
                               Rt = mean(Rt),
                               R0_min = quantile(R0, 0.025),
                               R0_q25 = quantile(R0, 0.25),
                               R0_q75 = quantile(R0, 0.75),
                               R0_max = quantile(R0, 0.975),
                               R0_median = median(R0),
                               R0 = mean(R0),
                               Reff_min = quantile(Reff, 0.025),
                               Reff_q25 = quantile(Reff, 0.25),
                               Reff_q75 = quantile(Reff, 0.75),
                               Reff_max = quantile(Reff, 0.975),
                               Reff_median = median(Reff),
                               Reff = mean(Reff)))
  
  min_date <- min(as.Date(out$replicate_parameters$start_date))
  
  country_plot <- function(vjust = -1.2) {
    ggplot(sum_rt %>% filter(
      date > min_date & date <= as.Date(as.character(date_0+as.numeric(lubridate::wday(date_0)))))) +
      geom_ribbon(mapping = aes(x=date, ymin=R0_min, ymax = R0_max, group = iso), fill = "#8cbbca") +
      geom_ribbon(mapping = aes(x = date, ymin = R0_q25, ymax = R0_q75, group = iso), fill = "#3f8da7") +
      geom_ribbon(mapping = aes(x=date, ymin=Reff_min, ymax = Reff_max, group = iso), fill = "#96c4aa") +
      geom_ribbon(mapping = aes(x = date, ymin = Reff_q25, ymax = Reff_q75, group = iso), fill = "#48996b") +
      geom_line(mapping = aes(x = date, y = Reff_median), color = "#48996b") +
      geom_hline(yintercept = 1, linetype = "dashed") +
      geom_hline(yintercept = sum_rt$R0_median[1], linetype = "dashed") +
      theme_bw() +
      theme(axis.text = element_text(size=12)) +
      xlab("") +
      ylab("Reff") +
      scale_x_date(breaks = "2 weeks",
                   limits = as.Date(c(as.character(min_date),
                                      as.character(date_0+as.numeric(lubridate::wday(date_0))))), 
                   date_labels = "%d %b",
                   expand = c(0,0)) + 
      theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1, colour = "black"),
            panel.grid.major = element_blank(), 
            panel.grid.minor = element_blank(),
            panel.border = element_blank(),
            panel.background = element_blank(), 
            axis.line = element_line(colour = "black")
      )
  }
  
  
  res <- list("plot" = suppressWarnings(country_plot()), "rts" = sum_rt)
  return(res)  
}