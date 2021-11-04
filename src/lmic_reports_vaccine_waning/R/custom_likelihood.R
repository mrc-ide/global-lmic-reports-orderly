#taken from squire, has added ability to update prob hospital multiplier,
# and duration in R
calc_loglikelihood_delta <- function(pars, data, squire_model, model_params,
                               pars_obs, n_particles,
                               forecast_days = 0, return = "ll",
                               Rt_args,
                               interventions, ...) {
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

  #----------------..
  # update the model params accordingly from new inputs
  #----------------..
  model_params$beta_set <- beta_set

  #----------------..
  # run the particle filter
  #----------------..
  if (inherits(squire_model, "stochastic")) {

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
