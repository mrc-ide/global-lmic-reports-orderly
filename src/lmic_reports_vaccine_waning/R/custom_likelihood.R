#taken from squire, has added ability to update prob hospital multiplier,
# and duration in R. As incorporates adjustment to cases for last few days
calc_loglikelihood_delta <- function(pars, data, squire_model, model_params,
                                     pars_obs, n_particles,
                                     forecast_days = 0, return = "ll",
                                     Rt_args,
                                     interventions,
                                     ...) {
  #----------------..
  # specify particle setup
  #----------------..
  switch(return,
         "full" = {
           save_particles <- TRUE
           full_output <- TRUE
           pf_return <- "sample"
         },
         "ll" = {
           save_particles <- FALSE
           forecast_days <- 0
           full_output <- FALSE
           pf_return <- "single"
         },
         {
           stop("Unknown return type to calc_loglikelihood")
         }
  )

  #----------------..
  # (potentially redundant) assertion
  #----------------..
  squire:::assert_in(c("R0", "start_date"), names(pars),
                     message = "Must specify R0, start date to infer")

  #----------------..
  # unpack current params
  #----------------..
  R0 <- pars[["R0"]]
  start_date <- pars[["start_date"]]

  # reporting fraction par if in pars list
  if("rf" %in% names(pars)) {
    assert_numeric(pars[["rf"]])
    pars_obs$phi_death <- pars[["rf"]]
  }

  #----------------..
  # more assertions
  #----------------..
  squire:::assert_pos(R0)
  squire:::assert_date(start_date)

  #----------------..
  # setup model based on inputs and interventions
  #----------------..
  R0_change <- interventions$R0_change
  date_R0_change <- interventions$date_R0_change
  date_contact_matrix_set_change <- interventions$date_contact_matrix_set_change
  date_ICU_bed_capacity_change <- interventions$date_ICU_bed_capacity_change
  date_hosp_bed_capacity_change <- interventions$date_hosp_bed_capacity_change
  date_vaccine_change <- interventions$date_vaccine_change
  date_vaccine_efficacy_infection_change <- interventions$date_vaccine_efficacy_infection_change
  date_vaccine_efficacy_disease_change <- interventions$date_vaccine_efficacy_disease_change

  # change betas
  if (is.null(date_R0_change)) {
    tt_beta <- 0
  } else {
    tt_list <- squire:::intervention_dates_for_odin(dates = date_R0_change,
                                                    change = R0_change,
                                                    start_date = start_date,
                                                    steps_per_day = round(1/model_params$dt),
                                                    starting_change = 1)
    model_params$tt_beta <- tt_list$tt
    R0_change <- tt_list$change
    date_R0_change <- tt_list$dates
  }

  # and contact matrixes
  if (is.null(date_contact_matrix_set_change)) {
    tt_contact_matrix <- 0
  } else {

    # here just provide positions for change and then use these to index mix_mat_set
    tt_list <- squire:::intervention_dates_for_odin(dates = date_contact_matrix_set_change,
                                                    change = seq_along(interventions$contact_matrix_set)[-1],
                                                    start_date = start_date,
                                                    steps_per_day = round(1/model_params$dt),
                                                    starting_change = 1)


    model_params$tt_matrix <- tt_list$tt
    model_params$mix_mat_set <- model_params$mix_mat_set[tt_list$change,,]
  }

  # and icu beds
  if (is.null(date_ICU_bed_capacity_change)) {
    tt_ICU_beds <- 0
  } else {
    tt_list <- squire:::intervention_dates_for_odin(dates = date_ICU_bed_capacity_change,
                                                    change = interventions$ICU_bed_capacity[-1],
                                                    start_date = start_date,
                                                    steps_per_day = round(1/model_params$dt),
                                                    starting_change = interventions$ICU_bed_capacity[1])
    model_params$tt_ICU_beds <- tt_list$tt
    model_params$ICU_beds <- tt_list$change
  }

  # and hosp beds
  if (is.null(date_hosp_bed_capacity_change)) {
    tt_hosp_beds <- 0
  } else {
    tt_list <- squire:::intervention_dates_for_odin(dates = date_hosp_bed_capacity_change,
                                                    change = interventions$hosp_bed_capacity[-1],
                                                    start_date = start_date,
                                                    steps_per_day = round(1/model_params$dt),
                                                    starting_change = interventions$hosp_bed_capacity[1])
    model_params$tt_hosp_beds <- tt_list$tt
    model_params$hosp_beds <- tt_list$change
  }

  # and vaccine coverage
  if (is.null(date_vaccine_change)) {
    tt_vaccine <- 0
  } else {
    tt_list <- squire:::intervention_dates_for_odin(dates = date_vaccine_change,
                                                    change = interventions$max_vaccine[-1],
                                                    start_date = start_date,
                                                    steps_per_day = round(1/model_params$dt),
                                                    starting_change = interventions$max_vaccine[1])
    model_params$tt_vaccine <- tt_list$tt
    model_params$max_vaccine <- tt_list$change
  }

  # and vaccine efficacy infection
  if (is.null(date_vaccine_efficacy_infection_change)) {
    tt_vaccine_efficacy_infection <- 0
  } else {

    # here we just pass the change as a position vector as we need to then
    # index the array of vaccine efficacies
    tt_list <- squire:::intervention_dates_for_odin(dates = date_vaccine_efficacy_infection_change,
                                                    change = seq_along(interventions$vaccine_efficacy_infection)[-1],
                                                    start_date = start_date,
                                                    steps_per_day = round(1/model_params$dt),
                                                    starting_change = 1)

    model_params$tt_vaccine_efficacy_infection <- tt_list$tt

    # here we have to not index the array by the postion vectors that are reutrned by intervention_dates_for_odin
    model_params$vaccine_efficacy_infection <- model_params$vaccine_efficacy_infection[tt_list$change,,]
  }

  # and vaccine efficacy disease
  if (is.null(date_vaccine_efficacy_disease_change)) {
    tt_vaccine_efficacy_disease <- 0
  } else {

    # here we just pass the change as a position vector as we need to then
    # index the array of vaccine efficacies
    tt_list <- squire:::intervention_dates_for_odin(dates = date_vaccine_efficacy_disease_change,
                                                    change = seq_along(interventions$vaccine_efficacy_disease)[-1],
                                                    start_date = start_date,
                                                    steps_per_day = round(1/model_params$dt),
                                                    starting_change = 1)

    model_params$tt_vaccine_efficacy_disease <- tt_list$tt

    # here we have to not index the array by the position vectors that are returned by intervention_dates_for_odin
    model_params$prob_hosp <- model_params$prob_hosp[tt_list$change,,]
  }

  #nimue specific functions (currently stored in pars_obs, may one day be interventions with the rest)
  if(!is.null(pars_obs$delta_adjust)){
    #check if we are using NIMUE
    if(!("nimue_model" %in% class(nimue:::nimue_deterministic_model()))){
      stop("Can only specify delta adjustments (pars_obs and dur_R change) for Nimue")
    }
    #update dur_R if needed
    if(!is.null(pars_obs$delta_adjust$dur_R)){
      gamma_R <- 2/pars_obs$delta_adjust$dur_R
      if(any(!gamma_R %in% model_params$gamma_R)){
        tt_list <- squire:::intervention_dates_for_odin(dates = pars_obs$delta_adjust$date_dur_R_change,
                                                        change = seq_along(gamma_R)[-1],
                                                        start_date = start_date,
                                                        steps_per_day = round(1/model_params$dt),
                                                        starting_change = 1)
        model_params$tt_dur_R <- tt_list$tt
        model_params$gamma_R <- gamma_R[tt_list$change]
      }
    }
    #update prob_hosp if needed
    if(!is.null(pars_obs$delta_adjust$prob_hosp_multiplier)){
      if(any(!pars_obs$delta_adjust$prob_hosp_multiplier %in% model_params$prob_hosp_multiplier)){
        tt_list <- squire:::intervention_dates_for_odin(dates = pars_obs$delta_adjust$date_prob_hosp_multiplier_change,
                                                        change = seq_along(pars_obs$delta_adjust$prob_hosp_multiplier)[-1],
                                                        start_date = start_date,
                                                        steps_per_day = round(1/model_params$dt),
                                                        starting_change = 1)
        model_params$tt_prob_hosp_multiplier <- tt_list$tt
        model_params$prob_hosp_multiplier <- pars_obs$delta_adjust$prob_hosp_multiplier[tt_list$change]
      }
    }
  }

  #--------------------..
  # update new R0s based on R0_change and R0_date_change, and Meff_date_change
  #--------------------..
  # and now get new R0s for the R0
  R0 <- squire:::evaluate_Rt_pmcmc(R0_change = R0_change,
                                   R0 = R0,
                                   date_R0_change = date_R0_change,
                                   pars = pars,
                                   Rt_args = Rt_args)

  # which allow us to work out our beta
  beta_set <- squire:::beta_est(squire_model = squire_model,
                                model_params = model_params,
                                R0 = R0)
  if(pars_obs$cases_fitting){
    #we'll save the date that fitting should start and the date from which we'll
    #estimate the reporting fraction
    cases_fitting_start_date <- min(
      data$date[
        data$date > (max(data$date) - pars_obs$cases_days)
      ]
    )
    cases_reporting_start_date <- min(
      data$date[
        data$date > cases_fitting_start_date - Rt_args$Rt_rw_duration
      ]
    )
  }

  #----------------..
  # update the model params accordingly from new inputs
  #----------------..
  model_params$beta_set <- beta_set

  #----------------..
  # run the particle filter
  #----------------..
  if(inherits(squire_model, "deterministic") & pars_obs$cases_fitting){
    pf_result <- run_deterministic_comparison_cases(
      data = data,
      squire_model = squire_model,
      model_params = model_params,
      model_start_date = start_date,
      cases_fitting_start_date,
      cases_reporting_start_date,
      obs_params = pars_obs,
      forecast_days = forecast_days,
      save_history = save_particles,
      return = pf_return
    )
  } else if (inherits(squire_model, "stochastic")) {

    pf_result <- squire:::run_particle_filter(data = data,
                                              squire_model = squire_model,
                                              model_params = model_params,
                                              model_start_date = start_date,
                                              obs_params = pars_obs,
                                              n_particles = n_particles,
                                              forecast_days = forecast_days,
                                              save_particles = save_particles,
                                              full_output = full_output,
                                              return = pf_return)

  } else if (inherits(squire_model, "deterministic")) {

    pf_result <- squire:::run_deterministic_comparison(data = data,
                                                       squire_model = squire_model,
                                                       model_params = model_params,
                                                       model_start_date = start_date,
                                                       obs_params = pars_obs,
                                                       forecast_days = forecast_days,
                                                       save_history = save_particles,
                                                       return = pf_return)

  }

  # out
  pf_result
}

run_deterministic_comparison_cases <- function(data,
                                               squire_model,
                                               model_params,
                                               model_start_date = "2020-02-02",
                                               cases_fitting_start_date,
                                               cases_reporting_start_date,
                                               obs_params = list(phi_cases = 0.1,
                                                                 k_cases = 2,
                                                                 phi_death = 1,
                                                                 k_death = 2,
                                                                 exp_noise = 1e6),
                                               forecast_days = 0,
                                               save_history = FALSE,
                                               return = "ll") {

  # parameter checks
  if (!(return %in% c("full", "ll", "sample", "single"))) {
    stop("return argument must be full, ll, sample", "single")
  }
  if (as.Date(data$date[data$deaths > 0][1], "%Y-%m-%d") < as.Date(model_start_date, "%Y-%m-%d")) {
    stop("Model start date is later than data start date")
  }

  # convert data into particle-filter form
  data <- squire:::particle_filter_data(data = data,
                                        start_date = model_start_date,
                                        steps_per_day = round(1 / model_params$dt))

  model_params$tt_beta <- round(model_params$tt_beta*model_params$dt)
  model_params$tt_contact_matrix <- round(model_params$tt_contact_matrix*model_params$dt)
  model_params$tt_hosp_beds <- round(model_params$tt_hosp_beds*model_params$dt)
  model_params$tt_ICU_beds <- round(model_params$tt_ICU_beds*model_params$dt)

  #set up model
  model_func <- squire_model$odin_model(user = model_params,
                                        unused_user_action = "ignore")

  # steps for the deterministic
  steps <- c(0, data$day_end)
  fore_steps <- seq(data$day_end[nrow(data)], length.out = forecast_days + 1L)
  steps <- unique(c(steps,fore_steps))

  # model run
  if("atol" %in% names(obs_params) && "rtol" %in% names(obs_params)) {
    assert_numeric(obs_params$atol)
    atol <- obs_params$atol
    assert_numeric(obs_params$rtol)
    rtol <- obs_params$rtol
  } else {
    atol <- 1e-6
    rtol <- 1e-6
  }

  out <- model_func$run(t = seq(0, tail(steps,1), 1), atol = atol, rtol = rtol)
  index <- squire:::odin_index(model_func)

  # get deaths for comparison
  Ds <- diff(rowSums(out[,index$D]))
  Ds <- Ds[data$day_end[-1]]
  Ds[Ds < 0] <- 0
  deaths <- data$deaths[-1]

  # calculate ll for deaths
  if (obs_params$treated_deaths_only) {

    Ds_heathcare <- diff(rowSums(out[,index$D_get]))
    Ds_heathcare <- Ds_heathcare[data$day_end[-1]]
    ll <- squire:::ll_nbinom(deaths, Ds_heathcare, obs_params$phi_death, obs_params$k_death, obs_params$exp_noise)

  } else {

    ll <- squire:::ll_nbinom(deaths, Ds, obs_params$phi_death, obs_params$k_death, obs_params$exp_noise)

  }

  # calculate ll for the seroprevalence
  lls <- 0
  if("sero_df" %in% names(obs_params) && "sero_det" %in% names(obs_params)) {

    sero_df <- obs_params$sero_df
    sero_det <- obs_params$sero_det

    # put some checks here that sero_df is correctly formatted
    check_sero_df(sero_df)

    # were there actually seroprevalence data points to compare against
    if(nrow(sero_df) > 0) {

      sero_at_date <- function(date, symptoms, det, dates, N) {

        di <- which(dates == date)
        if(length(di) > 0) {
          to_sum <- tail(symptoms[seq_len(di)], length(det))
          min(sum(rev(to_sum)*head(det, length(to_sum)), na.rm=TRUE)/N, 0.99)
        } else {
          0
        }

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

    }

  }

  # calculate ll for the PCR prevalence
  llp <- 0
  if("pcr_df" %in% names(obs_params) && "pcr_det" %in% names(obs_params)) {

    pcr_df <- obs_params$pcr_df
    pcr_det <- obs_params$pcr_det

    # put some checks here that pcr_df is correctly formatted
    check_pcr_df(pcr_df)

    # were there actually pcr prevalence data points to compare against
    if(nrow(pcr_df) > 0) {

      pcr_at_date <- function(date, infections, det, dates, N) {

        di <- which(dates == date)
        if(length(di) > 0) {
          to_sum <- tail(infections[seq_len(di)], length(det))
          min(sum(rev(to_sum)*head(det, length(to_sum)), na.rm=TRUE)/N, 0.99)
        } else {
          0
        }

      }

      # get infection incidence
      infections <- c(0,rowSums(out[-nrow(out),index$S]-out[-1,index$S]))

      # dates of incidence, pop size and dates of pcr surveys
      dates <- data$date[[1]] + seq_len(nrow(out)) - 1L
      N <- sum(model_params$population)
      pcr_dates <- list(pcr_df$date_end, pcr_df$date_start, pcr_df$date_start + as.integer((pcr_df$date_end - pcr_df$date_start)/2))
      unq_pcr_dates <- unique(c(pcr_df$date_end, pcr_df$date_start, pcr_df$date_start + as.integer((pcr_df$date_end - pcr_df$date_start)/2)))
      det <- obs_params$pcr_det

      # estimate model pcrprev
      pcr_model <- vapply(unq_pcr_dates, pcr_at_date, numeric(1), infections, det, dates, N)
      pcr_model_mat <- do.call(cbind,lapply(pcr_dates, function(x) {pcr_model[match(x, unq_pcr_dates)]}))

      # likelihood of model obvs
      llp <- rowMeans(dbinom(pcr_df$pcr_pos, pcr_df$samples, pcr_model_mat, log = TRUE))

    }

  }

  # calculate ll for the cases for last few days
  llc <- 0
  if(obs_params$cases_fitting){
    #get the relevant data
    data_reporting <- data %>%
      dplyr::filter(
        date >= cases_reporting_start_date & date < cases_fitting_start_date
      )
    data_fitting <- data %>%
      dplyr::filter(
        date >= cases_fitting_start_date
      )
    #get infections
    model_infections <- rowSums(out[,index$E2]) * model_params$gamma_E
    #check that cases are formatted correctly and exists
    if(all(is.na(data_reporting$cases)) |
       all(is.na(data_fitting$cases)) |
       any(c(data_reporting$cases, data_fitting$cases) < 0)){
      stop("Data for cases not formatted correctly or negative, please fix or disable fitting to cases")
    }
    #estimate reporting fraction
    est_reporting_fraction <- mean(
      data_reporting$cases/model_infections[data_reporting$day_end],
      na.rm = TRUE
    )
    #note that though this can be larger than 1, then the model output is small,
    #this should correct as the model fits, this should still follow the shape of
    #the infections
    #calculate log likelihood with estimated reporting fraction and given variance
    cases_data <- na.omit(cbind(
      data_fitting$cases,
      model_infections[data_fitting$day_end]
    ))
    llc <- squire:::ll_nbinom(
      cases_data[,1], cases_data[,2], est_reporting_fraction,
      obs_params$k_cases, obs_params$exp_noise
    )
    #note that we still fit to deaths in this time period
  }

  # format the out object
  date <- data$date[[1]] + seq_len(nrow(out)) - 1L
  rownames(out) <- as.character(date)
  attr(out, "date") <- date

  # format similar to particle_filter nomenclature
  pf_results <- list()
  pf_results$log_likelihood <- sum(ll) + sum(lls) + sum(llp) + sum(llc)

  # single returns final state
  if (save_history) {
    pf_results$states <- out
  } else if (return == "single") {
    pf_results$sample_state <- out[nrow(out), ]
  }

  # create returned object
  if (return == "ll") {
    ret <- pf_results$log_likelihood
  } else if (return == "sample") {
    ret <- pf_results$states
  } else if (return == "single" || return == "full") {
    ret <- pf_results
  }

  ret
}


#a function to make a plot looking at adjusted cases fitting
compare_adjustment_plot <- function(out, grad_dur, extend_past = 10){
  #get the model infections
  plotting_data <- nimue_format(out, var_select = "infections") %>%
    filter(t >= (max(t, na.rm = TRUE) - grad_dur - extend_past)) %>%
    rename(model_infections = y) %>%
    select(replicate, t, model_infections) %>%
    mutate(
      reporting = t < max(t, na.rm = TRUE) - grad_dur
    ) %>%
    #add the real data
    left_join(
      out$pmcmc_results$inputs$data %>%
        mutate(t = as.numeric(date - max(date))) %>%
        select(!deaths),
      by = "t"
    ) %>%
    #calculate reporting fraction
    group_by(replicate) %>%
    mutate(
      reporting_fraction = mean(
        if_else(reporting,
                cases/model_infections,
                as.numeric(NA)),
        na.rm = TRUE
      )
    ) %>%
    #convert cases to the ones we fit too
    mutate(
      adjusted_cases = cases/reporting_fraction
    )
  #plot
  ggplot(data = plotting_data, aes(x = date, y = adjusted_cases)) +
    geom_vline(aes(xintercept = max(date)-grad_dur),
               linetype = "dashed") +
    geom_point(alpha = 0.25) +
    geom_line(
      inherit.aes = FALSE,
      data = plotting_data %>%
        group_by(date) %>%
        summarise(
          infections_med = median(model_infections , na.rm = TRUE)
        ),
      aes(x = date, y = infections_med),
      colour = "red"
    ) +
    geom_ribbon(
      inherit.aes = FALSE,
      data = plotting_data %>%
        group_by(date) %>%
        summarise(
          infections_025 = quantile(model_infections , 0.025, na.rm = TRUE),
          infections_975 = quantile(model_infections , 0.975, na.rm = TRUE)
        ),
      aes(x = date, ymin = infections_025,
          ymax = infections_975),
      alpha =0.25,
      fill = "red"
    ) +
    ggplot2::theme_bw() +
    labs(x = "Date", y = "Infections, cases adjusted\nby reporting fraction") +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1, colour = "black"),
                   axis.title.x = ggplot2::element_blank(),
                   panel.grid.major.x = ggplot2::element_blank(),
                   panel.grid.minor.x = ggplot2::element_blank(),
                   panel.border = ggplot2::element_blank(),
                   panel.background = ggplot2::element_blank(),
                   axis.line = ggplot2::element_line(colour = "black")
    ) + ggplot2::theme(legend.position = "top",
                       legend.justification = c(0,1),
                       legend.direction = "horizontal") +
    lims(y = c(0, max(c(plotting_data$adjusted_cases,
                        plotting_data$model_infections))))
}
