orderly_id <- tryCatch(orderly::orderly_run_info()$id,
                       error = function(e) "<id>")

print(sessionInfo())
RhpcBLASctl::blas_set_num_threads(1L)
RhpcBLASctl::omp_set_num_threads(1L)

document <- as.logical(document)

# pandoc linking
if(file.exists("N:\\lmic_fitting\\pandoc") & document) {
  rmarkdown:::set_pandoc_info("N:\\lmic_fitting\\pandoc")
  Sys.setenv(RSTUDIO_PANDOC="N:\\lmic_fitting\\pandoc")
  tinytex::use_tinytex("N:\\lmic_fitting\\TinyTex")
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
system(paste0("echo Vaccine Reports for  ",iso3c))
if(!identical(seed, FALSE)){
  set.seed(seed)
}

# format user provided arguments correctly
date <- as.Date(date)

## Get the excess mortality estimates from the economist
excess_deaths <- readRDS("excess_deaths.Rds") %>%
  rename(iso = iso3c) %>%
  filter(iso == iso3c)

country <- squire::population$country[match(iso3c, squire::population$iso3c)[1]]

## MAIN LOOP IS ONLY FOR THOSE WITH DEATHS
if(sum(excess_deaths$deaths) > 10) {

  # get the raw data correct
  data <- excess_deaths %>%
    select(!iso) %>%
    arrange(date_start)

  if(iso3c == "TKM"){
    #remove starting deaths before major waves
    if(sum(cumsum(data$deaths) < 40) > 15){
      data <- data %>%
        filter(cumsum(deaths) > 40)
    }
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

  # conduct unmitigated
  pop <- squire::get_population(country)

  ## -----------------------------------------------------------------------------
  ## Step 2: Fit model
  ## -----------------------------------------------------------------------------

  squire_model <- squire.page:::nimue_booster_model()
  start_date <-  data$date_start[1] - 30
  end_date <- date

  ## -----------------------------------------------------------------------------
  ## Step 2a: Invariant Parameters
  ## -----------------------------------------------------------------------------

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

  tt_vaccine <- c(0, as.numeric(vacc_inputs$date_vaccine_change - start_date))

  parameters$first_doses <- c(0, vacc_inputs$first_doses)
  parameters$tt_first_doses <- tt_vaccine
  parameters$second_doses <- c(0, vacc_inputs$second_doses)
  parameters$tt_second_doses <- tt_vaccine
  parameters$booster_doses <- c(0, vacc_inputs$booster_doses)
  parameters$tt_booster_doses <- tt_vaccine
  rm(tt_vaccine)
  parameters$vaccine_coverage_mat <- vacc_inputs$vaccine_coverage_mat
  parameters$rel_infectiousness_vaccinated <- c(0.5)

  ## -----------------------------------------------------------------------------
  ## Step 2b: Calculate sampled parameters
  ## -----------------------------------------------------------------------------

  variants_to_model <- c("Delta", "Omicron", "Omicron Sub-Variant")

  #load inputs
  sample_vaccine_efficacies <- readRDS("vaccine_params.Rds")$sample_vaccine_efficacies
  variant_timings <- readRDS("variant_timings.Rds")[[iso3c]] %>%
    filter(variant %in% variants_to_model & !is.na(start_date)) %>%
    select(!mid_date) %>%
    arrange(start_date) %>%
    mutate(end_date = if_else(
      .data$end_date > lead(.data$start_date, n = 1, default = as_date(max(.data$end_date) + 1)),
      lead(.data$start_date, 1, default = as_date(max(.data$end_date) + 1)),
      .data$end_date
    ))
  variants_to_model <- variant_timings$variant

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
                                          dur_R, start_date)
    pars$dur_R <- dur_R_change$var
    pars$tt_dur_R <- dur_R_change$tt

    prob_hosp_multiplier_change <- multiplier_changes_over_time(variant_timings,
                                                                prob_hosp_multiplier, x, start_date)
    pars$prob_hosp_multiplier <- prob_hosp_multiplier_change$var
    pars$tt_prob_hosp_multiplier <- prob_hosp_multiplier_change$tt

    prob_severe_multiplier_change <- multiplier_changes_over_time(variant_timings,
                                                                  prob_severe_multiplier, x, start_date)
    pars$prob_severe_multiplier <- prob_severe_multiplier_change$var
    pars$tt_prob_severe_multiplier <- prob_severe_multiplier_change$tt

    #Vaccine Efficacies

    vacc_inputs <- vaccine_eff_over_time(variant_timings, variant_ve, x, start_date)
    pars$dur_V <- vacc_inputs$dur_V
    pars$vaccine_efficacy_infection <- vacc_inputs$vaccine_efficacy_infection
    pars$vaccine_efficacy_disease <- vacc_inputs$vaccine_efficacy_disease
    pars$tt_dur_V <- pars$tt_vaccine_efficacy_infection <-
      pars$tt_vaccine_efficacy_disease <- vacc_inputs$tt

    pars
  })

  ## -----------------------------------------------------------------------------
  ## Step 2c: Fitting Parameters
  ## -----------------------------------------------------------------------------

  #load in country specific default parameters
  fitting_params <- readRDS("fitting_params.Rds")[[iso3c]]
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

  ## -----------------------------------------------------------------------------
  ## Step 2d: Actual run
  ## -----------------------------------------------------------------------------

  #use parallel if asked for
  if(parallel){
    future::plan(future::multisession())
  }

  out <- rt_optimise(
    data = data,
    distribution = distribution,
    squire_model = squire_model,
    parameters = parameters,
    start_date = start_date,
    parallel = parallel,
    rt_spacing = 14,
    initial_infections_interval = initial_infections_interval,
    n_particles = n_particles,
    k = k,
    rt_interval = rt_interval
  )

  out_trim <- squire.page::trim_rt_optimise(out, 0.5)

  if(length(out_trim$samples) == 0){
    warning("No suitable trajectories calculated")
  } else {
    out <- out_trim
  }
  rm(out_trim)

  #Stop using parallel, furrr doesn't like something (maybe model object)
  Sys.setenv(SQUIRE_PARALLEL_DEBUG = "TRUE")
  if(parallel){
    future::plan(future::sequential())
  }

  ## -----------------------------------------------------------------------------
  ## Step 3: Summarise Fits
  ## -----------------------------------------------------------------------------

  ## summarise what we have
  #originally had a series of plots to summarise parameter distributions,
  #at the moment there are far to many so its something to work on

  title <- cowplot::ggdraw() +
    cowplot::draw_label(
      paste0(country, ", ", iso3c),
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
      date_0 = start_date, forecast = forecast,
      single = TRUE, full = TRUE) +
      theme(legend.position = "none")))

  suppressMessages(suppressWarnings(
    cas_plot <- cases_plot_single(
      df = nimue_format(out, "infections", date_0 = start_date),
      data = out$inputs$data,
      date = date,
      date_0 = start_date) +
      theme(legend.position = "none")))

  rtp <- rt_plot_immunity(out, R0_plot = TRUE)
  rtp2 <- rt_plot_immunity(out, R0_plot = FALSE)

  date_range <- as.Date(c(start_date, end_date))

  vt <- variant_timings_plot(variant_timings, date_range)

  suppressMessages(suppressWarnings(
    combined <- cowplot::plot_grid(
      header,
      vt + scale_x_date(date_breaks = "3 month", date_labels = "%b"),
      rtp2$plot + scale_x_date(date_breaks = "3 month", date_labels = "%b" ,limits = date_range),
      d + scale_x_date(date_breaks = "3 month", date_labels = "%b" ,limits = date_range),
      ncol=1,
      rel_heights = c(1.25, 5, 10, 10))
  ))

  ggsave("fitting.pdf",width=24, height=12,
         combined)


  ## Save the grid out object

  #remove all dr jacoby data, no reason we should need it

  ## now let's trim the out for saving really small
  output_temp <- out$output
  out$output <- NULL
  if(inherits(out, "rt_optimised_trimmed")){
    output_excluded_temp <- out$excluded$output
    out$excluded$output <- NULL
  }
  saveRDS(out, "grid_out.rds")
  #reattach output
  out$output <- output_temp
  rm(output_temp)
  if(inherits(out, "rt_optimised_trimmed")){
    out$excluded$output <- output_excluded_temp
    rm(output_excluded_temp)
  }

  if(document){
    #the following is only relevant if producing documentation

    ## -----------------------------------------------------------------------------
    ## Step 4: Scenarios
    ## -----------------------------------------------------------------------------

    ## -----------------------------------------------------------------------------
    ## 4.1. Conduct our projections based on no surging
    ## -----------------------------------------------------------------------------

    ## Functions for working out the relative changes in R0 for given scenarios
    time_period <- 365

    ## We need to know work out vaccine doses and efficacy going forwards
    model_user_args <- extend_vaccine_inputs(vacc_inputs, time_period, out, end_date = end_date)
    #For now these are all identicall we so'll add them to $parameters.
    #I've left this framework in case we want to add randomness to the doses
    out_extended <- out
    for(x in names(model_user_args[[1]])){
      if(str_detect(x, "tt")){
        out_extended$parameters[[x]] <- c(out$parameters[[x]], as.numeric(end_date - start_date) + 1)
      } else {
        out_extended$parameters[[x]] <- c(out$parameters[[x]], model_user_args[[1]][[x]])
      }
    }
    model_user_args <- map(model_user_args, ~list())
    rm(x, out)

    # Maintaining the current set of measures for a further 3 months
    maintain_scenario_leg <- squire.page::projections(out_extended,
                                     R0_change = c(1),
                                     tt_R0 = c(0),
                                     time_period = time_period,
                                     model_user_args = model_user_args)
    r_maintain_scenario_leg <- r_list_format(maintain_scenario_leg, date_0 = start_date)
    rt_maintain_scenario_leg <- rt_creation_vaccine(maintain_scenario_leg, end_date + time_period)

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
    icu <- nimue_format(maintain_scenario_leg, "ICU_demand", date_0 = start_date)
    icu <- icu[icu$compartment == "ICU_demand",]
    icu_28 <- group_by(icu, t) %>%
      summarise(i_tot = mean(y, na.rm = TRUE),
                i_min = t_test_safe(y)$conf.int[1],
                i_max = t_test_safe(y)$conf.int[2])

    hosp <- nimue_format(maintain_scenario_leg, "hospital_demand", date_0 = start_date)
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

    rm(maintain_scenario_leg)

    # Get the optimistic/pessimistic scenarios
    Rt_futures <- get_future_Rt_optimised(out_extended, forcast_days = 28)

    optimistic_scenario <- squire.page::projections(out_extended,
                                       R0 = map(Rt_futures$optimistic, ~.x$R0),
                                       tt_R0 = map(Rt_futures$optimistic, ~.x$tt_R0),
                                       time_period = time_period,
                                       model_user_args = model_user_args)

    r_optimistic_scenario <- r_list_format(optimistic_scenario, start_date)
    rt_optimistic_scenario <- rt_creation_vaccine(optimistic_scenario, end_date + time_period)

    rm(optimistic_scenario)

    pessimistic_scenario <- squire.page::projections(out_extended,
                                    R0 = map(Rt_futures$pessimistic, ~.x$R0),
                                    tt_R0 = map(Rt_futures$pessimistic, ~.x$tt_R0),
                                    time_period = time_period,
                                    model_user_args = model_user_args)
    r_pessimistic_scenario <- r_list_format(pessimistic_scenario, start_date)
    rt_pessimistic_scenario <- rt_creation_vaccine(pessimistic_scenario, end_date + time_period)
    rm(pessimistic_scenario)

    central_scenario <- squire.page::projections(out_extended,
                                    R0 = map(Rt_futures$central, ~.x$R0),
                                    tt_R0 = map(Rt_futures$central, ~.x$tt_R0),
                                    time_period = time_period,
                                    model_user_args = model_user_args)
    r_central_scenario <- r_list_format(central_scenario, start_date)
    rt_central_scenario <- rt_creation_vaccine(central_scenario, end_date + time_period)
    rm(central_scenario)

    #legacy scenarios
    mitigation_scenario_leg <- squire.page::projections(out_extended,
                                           R0_change = c(0.5),
                                           tt_R0 = 0,
                                           time_period = time_period,
                                           model_user_args = model_user_args)
    r_mitigation_scenario_leg <- r_list_format(mitigation_scenario_leg, start_date)
    rt_mitigation_scenario_leg <- rt_creation_vaccine(mitigation_scenario_leg, end_date + time_period)
    rm(mitigation_scenario_leg)

    reverse_scenario_leg <- squire.page::projections(out_extended,
                                        R0_change = c(1.5),
                                        tt_R0 = c(0),
                                        time_period = time_period,
                                        model_user_args = model_user_args)
    r_reverse_scenario_leg <- r_list_format(reverse_scenario_leg, start_date)
    rt_reverse_scenario_leg <- rt_creation_vaccine(reverse_scenario_leg, end_date + time_period)
    rm(reverse_scenario_leg)

    ## -----------------------------------------------------------------------------
    ## 4.2. Investigating a capacity surge
    ## -----------------------------------------------------------------------------

    # update for surged healthcare
    out_surged <- out_extended
    rm(out_extended)

    out_surged$parameters$hosp_bed_capacity <- 1e10
    out_surged$parameters$ICU_bed_capacity <- 1e10

    out_surged <- generate_draws(out_surged)

    ## Now simulate the surged scenarios

    optimistic_scenario_surged <- projections(out_surged,
                                       R0 = map(Rt_futures$optimistic, ~.x$R0),
                                       tt_R0 = map(Rt_futures$optimistic, ~.x$tt_R0),
                                       time_period = time_period,
                                       model_user_args = model_user_args)

    r_optimistic_scenario_surged <- r_list_format(optimistic_scenario_surged, start_date)
    rt_optimistic_scenario_surged <- rt_creation_vaccine(optimistic_scenario_surged, end_date + time_period)

    rm(optimistic_scenario_surged)

    pessimistic_scenario_surged <- projections(out_surged,
                                        R0 = map(Rt_futures$pessimistic, ~.x$R0),
                                        tt_R0 = map(Rt_futures$pessimistic, ~.x$tt_R0),
                                        time_period = time_period,
                                        model_user_args = model_user_args)
    r_pessimistic_scenario_surged <- r_list_format(pessimistic_scenario_surged, start_date)
    rt_pessimistic_scenario_surged <- rt_creation_vaccine(pessimistic_scenario_surged, end_date + time_period)
    rm(pessimistic_scenario_surged)

    central_scenario_surged <- projections(out_surged,
                                    R0 = map(Rt_futures$central, ~.x$R0),
                                    tt_R0 = map(Rt_futures$central, ~.x$tt_R0),
                                    time_period = time_period,
                                    model_user_args = model_user_args)
    r_central_scenario_surged <- r_list_format(central_scenario_surged, start_date)
    rt_central_scenario_surged <- rt_creation_vaccine(central_scenario_surged, end_date + time_period)
    rm(central_scenario_surged)

    #legacy scenarios
    maintain_scenario_surged_leg <- projections(out_surged,
                                            R0_change = c(1),
                                            tt_R0 = c(0),
                                            time_period = time_period,
                                            model_user_args = model_user_args)

    r_maintain_scenario_surged_leg <- r_list_format(maintain_scenario_surged_leg, start_date)
    rt_maintain_scenario_surged_leg <- rt_creation_vaccine(maintain_scenario_surged_leg, end_date + time_period)
    rm(maintain_scenario_surged_leg)

    mitigation_scenario_surged_leg <- projections(out_surged,
                                                  R0_change = c(0.5),
                                                  tt_R0 = c(0),
                                                  time_period = time_period,
                                                  model_user_args = model_user_args)
    r_mitigation_scenario_surged_leg <- r_list_format(mitigation_scenario_surged_leg, start_date)
    rt_mitigation_scenario_surged_leg <- rt_creation_vaccine(mitigation_scenario_surged_leg, end_date + time_period)
    rm(mitigation_scenario_surged_leg)


    reverse_scenario_surged_leg <- projections(out_surged,
                                               R0_change = c(1.5),
                                               tt_R0 = c(0),
                                               time_period = time_period,
                                               model_user_args = model_user_args)
    r_reverse_scenario_surged_leg <- r_list_format(reverse_scenario_surged_leg, start_date)
    rt_reverse_scenario_surged_leg <- rt_creation_vaccine(reverse_scenario_surged_leg, end_date + time_period)
    rm(reverse_scenario_surged_leg)

    rm(out_surged)

    o_list <- named_list(
      r_central_scenario,
      r_optimistic_scenario,
      r_pessimistic_scenario,
      r_central_scenario_surged,
      r_optimistic_scenario_surged,
      r_pessimistic_scenario_surged,
      r_maintain_scenario_leg,
      r_mitigation_scenario_leg,
      r_reverse_scenario_leg,
      r_maintain_scenario_surged_leg,
      r_mitigation_scenario_surged_leg,
      r_reverse_scenario_surged_leg
    )

    rt_list <- named_list(
      rt_central_scenario,
      rt_optimistic_scenario,
      rt_pessimistic_scenario,
      rt_central_scenario_surged,
      rt_optimistic_scenario_surged,
      rt_pessimistic_scenario_surged,
      rt_maintain_scenario_leg,
      rt_mitigation_scenario_leg,
      rt_reverse_scenario_leg,
      rt_maintain_scenario_surged_leg,
      rt_mitigation_scenario_surged_leg,
      rt_reverse_scenario_surged_leg
    )

    rt_futures_df <-
      summarise_rt_futures(Rt_futures)

    ## -----------------------------------------------------------------------------
    ## Step 5: Report
    ## -----------------------------------------------------------------------------

    # get data in correct format for plotting
    df_excess_deaths <- excess_deaths %>%
      filter(iso == iso3c) %>%
      # mutate(
      #   deaths = deaths/as.numeric(date_end - date_start)
      # ) %>%
      rename(iso3c = iso) %>%
      arrange(
        date_start
      )

    #also get data on cases
    df_cases <- readRDS("reported_covid.Rds") %>%
      rename(iso = iso3c) %>%
      filter(iso == iso3c) %>%
      rename(iso3c = iso) %>%
      arrange(date)

    # prepare reports
    options(tinytex.verbose = TRUE)
    rmarkdown::render("index.Rmd",
                      output_format = c("html_document","pdf_document"),
                      params = list("o_list" = o_list,
                                    "excess" = TRUE,
                                    "df_excess" = df_excess_deaths,
                                    "df_cases" = df_cases,
                                    "date_0" = start_date,
                                    "country" = country,
                                    "surging" = surging,
                                    "rt" = rtp2,
                                    "variants" = variant_timings,
                                    "date_range" = date_range,
                                    "rt_futures" = rt_futures_df),
                      output_options = list(pandoc_args = c(paste0("--metadata=title:",country," COVID-19 report "))))

    } else {
    # create empty files
    file.create("index.html", "index.pdf", "index.md",
                "summary_df.rds")
  }
} else {
## THIS IS THE ESFT LOOP FOR COUNTRIES WITH NO DEATHS CURRENTLY
## THIS COULD DO WITH UPDATING

  # What are the Rt values for each income group
  inc_R0s <- income_R0()
  inc_Rts <- income_Rt(date_0 = date)

  # what income group is this country
  income <- as.character(squire.page::get_income_group(iso3c))

  # And the upper and lower Rs for our scenarios
  R0 <- inc_R0s$R0[inc_R0s$income == income]
  Rt <- max(inc_Rts$Rt[inc_Rts$income == income], 1.2)

  # sim_args
  time_period_esft <- 366
  replicates_esft <- 1
  seeding_cases_esft <- 5

  vacc_inputs <- get_vaccine_inputs(iso3c)

  # set up vaccine inits
  vaccine_fitting_flag <- TRUE
  squire_model <- squire.page::nimue_booster_model()
  init <- init_state_nimue(deaths_removed = 0, iso3c,
                           seeding_cases = seeding_cases_esft,
                           vaccinated_already = sum(vacc_inputs$first_doses))
  init <- map(init, ~cbind(.x, rep(0, 17)))

  # Scenarios with capacity constraints
  # ---------------------------------------------------------------------------
  r_pessimistic_scenario <- squire.page:::run_booster(
    country = country,
    R0 = R0,
    time_period = time_period_esft,
    seeding_cases = seeding_cases_esft,
    init = init,
    first_doses =  as.integer(mean(tail(vacc_inputs$first_doses,7))),
    second_doses =  as.integer(mean(tail(vacc_inputs$second_doses,7))),
    booster_doses =  as.integer(mean(tail(vacc_inputs$booster_doses,7))),
    vaccine_coverage_mat = vacc_inputs$vaccine_coverage_mat
  )

  r_optimistic_scenario <- squire.page:::run_booster(
    country = country,
    R0 = Rt,
    time_period = time_period_esft,
    seeding_cases = seeding_cases_esft,
    init = init,
    first_doses =  as.integer(mean(tail(vacc_inputs$first_doses,7))),
    second_doses =  as.integer(mean(tail(vacc_inputs$second_doses,7))),
    booster_doses =  as.integer(mean(tail(vacc_inputs$booster_doses,7))),
    vaccine_coverage_mat = vacc_inputs$vaccine_coverage_mat
  )

  r_central_scenario <- squire.page:::run_booster(
    country = country,
    R0 = Rt + ((R0 - Rt)/2),
    time_period = time_period_esft,
    seeding_cases = seeding_cases_esft,
    init = init,
    first_doses =  as.integer(mean(tail(vacc_inputs$first_doses,7))),
    second_doses =  as.integer(mean(tail(vacc_inputs$second_doses,7))),
    booster_doses =  as.integer(mean(tail(vacc_inputs$booster_doses,7))),
    vaccine_coverage_mat = vacc_inputs$vaccine_coverage_mat
  )
  output_temp <- r_central_scenario$output
  saveRDS(r_central_scenario, "grid_out.rds")
  r_central_scenario$output <- output_temp
  rm(output_temp)

  # Scenarios without capacity constraints
  # ---------------------------------------------------------------------------
  r_pessimistic_scenario_surged <- squire.page:::run_booster(
    country = country,
    R0 = R0,
    time_period = time_period_esft,
    seeding_cases = seeding_cases_esft,
    init = init,
    first_doses =  as.integer(mean(tail(vacc_inputs$first_doses,7))),
    second_doses =  as.integer(mean(tail(vacc_inputs$second_doses,7))),
    booster_doses =  as.integer(mean(tail(vacc_inputs$booster_doses,7))),
    vaccine_coverage_mat = vacc_inputs$vaccine_coverage_mat,
    hosp_bed_capacity = 1e10,
    ICU_bed_capacity = 1e10
  )

  r_optimistic_scenario_surged <- squire.page:::run_booster(
    country = country,
    R0 = R0,
    time_period = time_period_esft,
    seeding_cases = seeding_cases_esft,
    init = init,
    first_doses =  as.integer(mean(tail(vacc_inputs$first_doses,7))),
    second_doses =  as.integer(mean(tail(vacc_inputs$second_doses,7))),
    booster_doses =  as.integer(mean(tail(vacc_inputs$booster_doses,7))),
    vaccine_coverage_mat = vacc_inputs$vaccine_coverage_mat,
    hosp_bed_capacity = 1e10,
    ICU_bed_capacity = 1e10
  )

  r_central_scenario_surged <- squire.page:::run_booster(
    country = country,
    R0 = R0,
    time_period = time_period_esft,
    seeding_cases = seeding_cases_esft,
    init = init,
    first_doses =  as.integer(mean(tail(vacc_inputs$first_doses,7))),
    second_doses =  as.integer(mean(tail(vacc_inputs$second_doses,7))),
    booster_doses =  as.integer(mean(tail(vacc_inputs$booster_doses,7))),
    vaccine_coverage_mat = vacc_inputs$vaccine_coverage_mat,
    hosp_bed_capacity = 1e10,
    ICU_bed_capacity = 1e10
  )

  # bundle into a list
  r_list <-
    named_list(
      r_central_scenario,
      r_optimistic_scenario,
      r_pessimistic_scenario,
      r_central_scenario_surged,
      r_optimistic_scenario_surged,
      r_pessimistic_scenario_surged
    )

  # And adjust their time variable so that we have t = 0 as today
  r_list_pass <- lapply(r_list, function(x) {
    for(i in seq_len(dim(x$output)[3])) {
      x$output[,"time",i] <- x$output[,"time",i] - 1
    }
    rownames(x$output) <- as.character(seq.Date(date, date + nrow(x$output) -1, 1))
    return(x)
  })

  ## Lastly make up some outputs here to pass orderly
  file.create("index.html", "index.pdf", "index.md",
              "summary_df.rds",
              "fitting.pdf")

  # major summaries
  o_list <- lapply(r_list_pass, r_list_format, date)
  rt_list <- lapply(
    r_list_pass,
    rt_creation_simple_out_vaccine, date, date+89
  )
  names(rt_list) <- str_replace(names(rt_list), "r_", "rt_")
}

if(document){
  sim_end_date <- as.Date(date)+90

  # summarise the projections
  data_sum <- map(o_list, function(pd){

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
      dplyr::filter(.data$date < sim_end_date) %>%
      dplyr::group_by(.data$date, .data$compartment) %>%
      dplyr::summarise(y_025 = stats::quantile(.data$y, 0.025),
                       y_25 = stats::quantile(.data$y, 0.25),
                       y_median = median(.data$y),
                       y_mean = mean(.data$y),
                       y_75 = stats::quantile(.data$y, 0.75),
                       y_975 = stats::quantile(.data$y, 0.975))

    return(as.data.frame(pds, stringsAsFactors = FALSE))
  })

  data_sum$r_central_scenario$scenario <- "Central"
  data_sum$r_optimistic_scenario$scenario <- "Optimistic"
  data_sum$r_pessimistic_scenario$scenario <- "Pessimistic"
  data_sum$r_central_scenario_surged$scenario <- "Surged Central"
  data_sum$r_optimistic_scenario_surged$scenario <- "Surged Optimistic"
  data_sum$r_pessimistic_scenario_surged$scenario <- "Surged Pessimistic"

  if(sum(excess_deaths$deaths) > 10){
    #these won't be present for the non-fitting countries
    data_sum$r_maintain_scenario_leg$scenario <- "Maintain Status Quo"
    data_sum$r_mitigation_scenario_leg$scenario <- "Additional 50% Reduction"
    data_sum$r_reverse_scenario_leg$scenario <- "Relax Interventions 50%"
    data_sum$r_maintain_scenario_surged_leg$scenario <- "Surged Maintain Status Quo"
    data_sum$r_mitigation_scenario_surged_leg$scenario <- "Surged Additional 50% Reduction"
    data_sum$r_reverse_scenario_surged_leg$scenario <- "Surged Relax Interventions 50%"
  } else {
    data_sum$r_maintain_scenario_leg <- data_sum$r_central_scenario
    data_sum$r_maintain_scenario_leg$scenario <- "Maintain Status Quo"
    data_sum$r_mitigation_scenario_leg <- data_sum$r_optimistic_scenario
    data_sum$r_mitigation_scenario_leg$scenario <- "Additional 50% Reduction"
    data_sum$r_reverse_scenario_leg <- data_sum$r_pessimistic_scenario
    data_sum$r_reverse_scenario_leg$scenario <- "Relax Interventions 50%"
    data_sum$r_maintain_scenario_surged_leg <- data_sum$r_central_scenario_surged
    data_sum$r_maintain_scenario_surged_leg$scenario <- "Surged Maintain Status Quo"
    data_sum$r_mitigation_scenario_surged_leg <- data_sum$r_optimistic_scenario_surged
    data_sum$r_mitigation_scenario_surged_leg$scenario <- "Surged Additional 50% Reduction"
    data_sum$r_reverse_scenario_surged_leg <- data_sum$r_pessimistic_scenario_surged
    data_sum$r_reverse_scenario_surged_leg$scenario <- "Surged Relax Interventions 50%"
  }

  rt_list$rt_central_scenario$scenario <- "Central"
  rt_list$rt_optimistic_scenario$scenario <- "Optimistic"
  rt_list$rt_pessimistic_scenario$scenario <- "Pessimistic"
  rt_list$rt_central_scenario_surged$scenario <- "Surged Central"
  rt_list$rt_optimistic_scenario_surged$scenario <- "Surged Optimistic"
  rt_list$rt_pessimistic_scenario_surged$scenario <- "Surged Pessimistic"
  if(sum(excess_deaths$deaths) > 10){
    #these won't be present for the non-fitting countries
    rt_list$rt_maintain_scenario_leg$scenario <- "Maintain Status Quo"
    rt_list$rt_mitigation_scenario_leg$scenario <- "Additional 50% Reduction"
    rt_list$rt_reverse_scenario_leg$scenario <- "Relax Interventions 50%"
    rt_list$rt_maintain_scenario_surged_leg$scenario <- "Surged Maintain Status Quo"
    rt_list$rt_mitigation_scenario_surged_leg$scenario <- "Surged Additional 50% Reduction"
    rt_list$rt_reverse_scenario_surged_leg$scenario <- "Surged Relax Interventions 50%"
  } else {
    rt_list$rt_maintain_scenario_leg <- rt_list$rt_central_scenario
    rt_list$rt_maintain_scenario_leg$scenario <- "Maintain Status Quo"
    rt_list$rt_mitigation_scenario_leg <- rt_list$rt_optimistic_scenario
    rt_list$rt_mitigation_scenario_leg$scenario <- "Additional 50% Reduction"
    rt_list$rt_reverse_scenario_leg <- rt_list$rt_pessimistic_scenario
    rt_list$rt_reverse_scenario_leg$scenario <- "Relax Interventions 50%"
    rt_list$rt_maintain_scenario_surged_leg <- rt_list$rt_central_scenario_surged
    rt_list$rt_maintain_scenario_surged_leg$scenario <- "Surged Maintain Status Quo"
    rt_list$rt_mitigation_scenario_surged_leg <- rt_list$rt_optimistic_scenario_surged
    rt_list$rt_mitigation_scenario_surged_leg$scenario <- "Surged Additional 50% Reduction"
    rt_list$rt_reverse_scenario_surged_leg <- rt_list$rt_pessimistic_scenario_surged
    rt_list$rt_reverse_scenario_surged_leg$scenario <- "Surged Relax Interventions 50%"
  }

  # combine and annotate
  data_sum <- do.call(rbind, data_sum)
  data_sum <- rbind(
    data_sum,
    do.call(rbind, rt_list) %>% filter(date <= max(data_sum$date))
  ) %>% arrange(date, scenario)
  rownames(data_sum) <- NULL

  # catch for hong kong and taiwan country name
  if (iso3c == "HKG") {
    country <- "Hong Kong"
  } else if (iso3c == "TWN") {
    country <- "Taiwan"
  } else if (iso3c == "MAC") {
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
  data_sum$death_calibrated <- sum(excess_deaths$deaths) == 0

  write.csv(data_sum, "projections.csv", row.names = FALSE, quote = FALSE)
} else {
  file.create("projections.csv")
}
