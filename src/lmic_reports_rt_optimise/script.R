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

country <- squire::population$country[match(iso3c, squire::population$iso3c)[1]]

## Get the excess mortality estimates from the economist
excess_deaths <- readRDS("excess_deaths.Rds") %>%
  rename(iso = iso3c) %>%
  filter(iso == iso3c) %>%
  arrange(date_start)

##Get the reported deaths
reported_deaths <- readRDS("reported_covid.Rds") %>%
  rename(iso = iso3c) %>%
  filter(iso == iso3c) %>%
  filter(cumsum(deaths) > 0) %>%
  transmute(
    iso = iso,
    date_start = date - 1,
    date_end = date,
    deaths = deaths
  ) %>%
  arrange(date_start)

death_limit <- 10

fit_excess <- sum(excess_deaths$deaths) > death_limit
fit_reported <- sum(reported_deaths$deaths) > death_limit

cases <- readRDS("reported_covid.Rds") %>%
  rename(iso = iso3c, d_date = date) %>%
  filter(iso == iso3c & d_date > date - 30*3 & d_date <= date) %>%
  transmute(date = d_date, detected_infections = cases) %>%
  arrange(date)
estimate_reported_cases <- sum(cases$detected_infections) > 0

#setup model and parameters
if(fit_excess | fit_reported){

  #get model start dates
  if(fit_excess){
    if(iso3c == "TKM"){
      #remove starting deaths before major waves
      if(sum(cumsum(excess_deaths$deaths) < 40) > 15){
        excess_deaths <- excess_deaths %>%
          filter(cumsum(excess_deaths) > 40)
      }
    }

    first_report <- which(excess_deaths$deaths>0)[1]
    missing <- which(excess_deaths$deaths == 0 | is.na(excess_deaths$deaths))
    to_remove <- missing[missing<first_report]
    if(length(to_remove) > 0) {
      if(length(to_remove) == (nrow(excess_deaths)-1)) {
        excess_deaths <- excess_deaths[-head(to_remove,-1),]
      } else {
        excess_deaths <- excess_deaths[-to_remove,]
      }
    }

    excess_start_date <-  excess_deaths$date_start[1] - 30
    first_start_date <- excess_start_date
  }
  if(fit_reported){
    first_report <- which(reported_deaths$deaths>0)[1]
    missing <- which(reported_deaths$deaths == 0 | is.na(reported_deaths$deaths))
    to_remove <- missing[missing<first_report]
    if(length(to_remove) > 0) {
      if(length(to_remove) == (nrow(reported_deaths)-1)) {
        reported_deaths <- reported_deaths[-head(to_remove,-1),]
      } else {
        reported_deaths <- reported_deaths[-to_remove,]
      }
    }

    reported_start_date <-  reported_deaths$date_start[1] - 30
    first_start_date <- reported_start_date
  }
  if(fit_excess & fit_reported){
    first_start_date <- min(c(reported_start_date, excess_start_date))
  }

  pop <- squire::get_population(country)
  squire_model <- squire.page:::nimue_booster_model()
  end_date <- date
  projection_period <- 365

  #Invariant parameters
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
  tt_vaccine <- as.numeric(c(first_start_date, vacc_inputs$date_vaccine_change) - first_start_date)

  parameters$primary_doses <- c(0, vacc_inputs$primary_doses)
  parameters$tt_primary_doses <- tt_vaccine
  parameters$second_dose_delay <- vacc_inputs$second_dose_delay
  parameters$booster_doses <- c(0, vacc_inputs$booster_doses)
  parameters$tt_booster_doses <- tt_vaccine
  rm(tt_vaccine)
  parameters$vaccine_coverage_mat <- vacc_inputs$vaccine_coverage_mat
  parameters$rel_infectiousness_vaccinated <- c(0.5)
  parameters$vaccine_booster_follow_up_coverage <- c(rep(0, 12), rep(1, 5))
  parameters$protection_delay_rate <- 1/7
  parameters$protection_delay_shape <- 2
  parameters$protection_delay_time <- as.numeric(date - first_start_date)

  #Variant dependant parameters
  variants_to_model <- c("Delta", "Omicron", "Omicron Sub-Variant")

  #load inputs
  sample_vaccine_efficacies <- readRDS("vaccine_params.Rds")$sample_vaccine_efficacies
  variant_timings <- readRDS("variant_timings.Rds")[[iso3c]] %>%
    filter(variant %in% variants_to_model & !is.na(start_date)) %>%
    select(!mid_date) %>%
    arrange(.data$start_date) %>%
    mutate(end_date = if_else(
      .data$end_date > lead(.data$start_date, n = 1, default = as_date(max(.data$end_date) + 1)),
      lead(.data$start_date, 1, default = as_date(max(.data$end_date) + 1)),
      .data$end_date
    ))
  variants_to_model <- variant_timings$variant

  # #add the changes to generation time to the fixed parameters
  # parameters <- append(parameters, estimate_generation_time(variant_timings, start_date)) %>%
  #                 append(estimate_healthcare_durations(variant_timings, start_date))
  parameters <-
    append(parameters, estimate_healthcare_durations(variant_timings, first_start_date))

  #Sample from random parameter
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
                                          dur_R, first_start_date)
    pars$dur_R <- dur_R_change$var
    pars$tt_dur_R <- dur_R_change$tt

    prob_hosp_multiplier_change <- multiplier_changes_over_time(variant_timings,
                                                                prob_hosp_multiplier, x, first_start_date)
    pars$prob_hosp_multiplier <- prob_hosp_multiplier_change$var
    pars$tt_prob_hosp_multiplier <- prob_hosp_multiplier_change$tt

    prob_severe_multiplier_change <- multiplier_changes_over_time(variant_timings,
                                                                  prob_severe_multiplier, x, first_start_date)
    pars$prob_severe_multiplier <- prob_severe_multiplier_change$var
    pars$tt_prob_severe_multiplier <- prob_severe_multiplier_change$tt

    #Vaccine Efficacies

    vacc_inputs <- vaccine_eff_over_time(variant_timings, variant_ve, x, first_start_date)
    pars$dur_V <- vacc_inputs$dur_V
    pars$vaccine_efficacy_infection <- vacc_inputs$vaccine_efficacy_infection
    pars$vaccine_efficacy_disease <- vacc_inputs$vaccine_efficacy_disease
    pars$tt_dur_V <- pars$tt_vaccine_efficacy_infection <-
      pars$tt_vaccine_efficacy_disease <- vacc_inputs$tt

    pars
  })
}

## Excess deaths loop
if(fit_excess){
  #ensure parameters are in correct format
  if(first_start_date != excess_start_date){
    difference_in_t <- as.numeric(excess_start_date - first_start_date)
    excess_parameters <- update_parameters(parameters, difference_in_t)
    excess_distribution <- update_distribution(distribution, difference_in_t)
  } else {
    excess_parameters <- parameters
    excess_distribution <- distribution
  }

  #fitting parameters
  #load in country specific default parameters
  fitting_params <- readRDS("fitting_params.Rds")[[iso3c]]$excess
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

  #run model
  #use parallel if asked for
  if(parallel){
    future::plan(future::multisession())
  }

  excess_out <- rt_optimise(
    data = excess_deaths,
    distribution = excess_distribution,
    squire_model = squire_model,
    parameters = excess_parameters,
    start_date = excess_start_date,
    parallel = parallel,
    rt_spacing = 14,
    initial_infections_interval = initial_infections_interval,
    n_particles = n_particles,
    k = k,
    rt_interval = rt_interval
  )

  excess_out <- trim_output(excess_out)

  #Stop using parallel, furrr doesn't like something (maybe model object)
  Sys.setenv(SQUIRE_PARALLEL_DEBUG = "TRUE")
  if(parallel){
    future::plan(future::sequential())
  }

  #save fitting plot
  summarise_fit("excess_fitting.pdf", excess_out, country, iso3c, end_date, excess_start_date)

  save_output(excess_out, "excess_out.Rds")

  if(document){
    #projections
    project_excess <- as.numeric(end_date + projection_period - max(excess_deaths$date_end))

    excess_out <- extend_output(excess_out, vacc_inputs, project_excess, excess_start_date, end_date)
    model_user_args <- map(seq_len(samples), ~list())

    # Maintaining the current set of measures for a further 3 months
    forwards_projection_excess <- squire.page::projections(excess_out,
                                                    R0 = map(excess_out$samples, ~tail(.x$R0, 1)),
                                                    tt_R0 = c(0),
                                                    time_period = project_excess,
                                                    model_user_args = model_user_args)
    r_forwards_projection_excess <- r_list_format(forwards_projection_excess, date_0 = excess_start_date)
    rt_forwards_projection_excess <- rt_creation_vaccine(forwards_projection_excess, end_date + projection_period)

    surging <- check_breached_capacity(forwards_projection_excess, country, excess_start_date)

    cases_milds_excess <- summarise_infection_types(forwards_projection_excess, end_date)

    age_cases_milds <- cases_milds_excess$age
    cases_milds_excess <- cases_milds_excess$total

    if(estimate_reported_cases){
      cases_milds_excess <- detected_cases_estimation(cases_milds_excess, cases)
      age_cases_milds <- estimate_cases_age(age_cases_milds, cases_milds_excess, end_date)
    }

    # age_split_excess <- get_age_output(forwards_projection_excess, excess_start_date, end_date) %>%
    #   left_join(
    #     age_cases_milds,
    #     by = c("replicate", "date", "age_group")
    #   ) %>%
    #   summarise_age_dependent() %>%
    #   mutate(
    #     fit_type = "Excess Mortality"
    #   )
    age_split_excess <- age_cases_milds %>%
      rename(non_hospitalised_infections = new_Mild,
             hospitalisations = new_hosp) %>%
      select(any_of(c("replicate", "date", "age_group", "non_hospitalised_infections", "non_hospitalised_cases", "hospitalisations")))%>%
      summarise_age_dependent() %>%
      mutate(
        fit_type = "Excess Mortality"
      )
    rm(age_cases_milds)
  }
} else {
  file.create(c("excess_fitting.pdf", "excess_out.Rds"))
}

## Report deaths loop
if(fit_reported){
  #ensure parameters are in correct format
  if(first_start_date != reported_start_date){
    difference_in_t <- as.numeric(reported_start_date - first_start_date)
    reported_parameters <- update_parameters(parameters, difference_in_t)
    reported_distribution <- update_distribution(distribution, difference_in_t)
  } else {
    reported_parameters <- parameters
    reported_distribution <- distribution
  }

  #fitting parameters
  #load in country specific default parameters
  fitting_params <- readRDS("fitting_params.Rds")[[iso3c]]$reported
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

  #run model
  #use parallel if asked for
  if(parallel){
    future::plan(future::multisession())
  }

  reported_out <- rt_optimise(
    data = reported_deaths,
    distribution = reported_distribution,
    squire_model = squire_model,
    parameters = reported_parameters,
    start_date = reported_start_date,
    parallel = parallel,
    rt_spacing = 14,
    initial_infections_interval = initial_infections_interval,
    n_particles = n_particles,
    k = k,
    rt_interval = rt_interval
  )

  reported_out <- trim_output(reported_out)

  #Stop using parallel, furrr doesn't like something (maybe model object)
  Sys.setenv(SQUIRE_PARALLEL_DEBUG = "TRUE")
  if(parallel){
    future::plan(future::sequential())
  }

  #save fitting plot
  summarise_fit("reported_fitting.pdf", reported_out, country, iso3c, end_date, reported_start_date)

  save_output(reported_out, "reported_out.Rds")

  if(document){
    #projections
    project_reported <- as.numeric(end_date + projection_period - max(reported_deaths$date_end))

    reported_out <- extend_output(reported_out, vacc_inputs, project_reported, reported_start_date, end_date)
    model_user_args <- map(seq_len(samples), ~list())

    # Maintaining the current set of measures for a further 3 months
    forwards_projection_reported <- squire.page::projections(reported_out,
                                                           R0 = map(reported_out$samples, ~tail(.x$R0, 1)),
                                                           tt_R0 = c(0),
                                                           time_period = project_reported,
                                                           model_user_args = model_user_args)
    r_forwards_projection_reported <- r_list_format(forwards_projection_reported, date_0 = reported_start_date)
    rt_forwards_projection_reported <- rt_creation_vaccine(forwards_projection_reported, end_date + projection_period)

    surging <- check_breached_capacity(forwards_projection_reported, country, reported_start_date)

    cases_milds_reported <- summarise_infection_types(forwards_projection_reported, end_date)

    age_cases_milds <- cases_milds_reported$age
    cases_milds_reported <- cases_milds_reported$total

    if(estimate_reported_cases){
      cases_milds_reported <- detected_cases_estimation(cases_milds_reported, cases)
      age_cases_milds <- estimate_cases_age(age_cases_milds, cases_milds_reported, end_date)
    }

    # age_split_reported <- get_age_output(forwards_projection_reported, reported_start_date, end_date) %>%
    #   left_join(
    #     age_cases_milds,
    #     by = c("replicate", "date", "age_group")
    #   ) %>%
    #   summarise_age_dependent() %>%
    #   mutate(
    #     fit_type = "Reported Deaths"
    #   )

    age_split_reported <- age_cases_milds %>%
      rename(non_hospitalised_infections = new_Mild,
             hospitalisations = new_hosp) %>%
      select(any_of(c("replicate", "date", "age_group", "non_hospitalised_infections", "non_hospitalised_cases", "hospitalisations")))%>%
      summarise_age_dependent() %>%
      mutate(
        fit_type = "Reported Deaths"
      )
    rm(age_cases_milds)
  }
} else {
  file.create(c("reported_fitting.pdf", "reported_out.Rds"))
}

#create rmd report
if(document & (fit_excess | fit_reported)){
  # get data in correct format for plotting
  df_excess_deaths <- excess_deaths %>%
    arrange(
      date_start
    )

  #also get data on cases
  df_cases <- readRDS("reported_covid.Rds") %>%
    rename(iso = iso3c) %>%
    filter(iso == iso3c) %>%
    rename(iso3c = iso) %>%
    arrange(date)

  list_projections <- list()
  list_date_range <- list()
  list_cases_milds <- list()
  list_rtp2 <- list()
  if(fit_excess){
    list_projections$excess <- r_forwards_projection_excess
    list_date_range$excess <- c(excess_start_date, end_date)
    list_cases_milds$excess <- cases_milds_excess
    list_rtp2$excess <- rt_plot_immunity(excess_out, R0_plot = FALSE)
  }
  if(fit_reported){
    list_projections$reported <- r_forwards_projection_reported
    list_date_range$reported <- c(reported_start_date, end_date)
    list_cases_milds$reported <- cases_milds_reported
    list_rtp2$reported <- rt_plot_immunity(reported_out, R0_plot = FALSE)
  }
  date_range <- list_date_range[[which.min(unlist(map(list_date_range, ~.x[1])))]]

  # prepare reports
  options(tinytex.verbose = TRUE)
  rmarkdown::render("index.Rmd",
                    output_format = c("html_document","pdf_document"),
                    params = list("list_projections" = list_projections,
                                  "date_range" = date_range,
                                  "list_cases_milds" = list_cases_milds,
                                  "list_rtp2" = list_rtp2,
                                  "fit_excess" = fit_excess,
                                  "fit_reported" = fit_reported,
                                  "df_excess" = df_excess_deaths,
                                  "df_cases" = df_cases,
                                  "date_0" = first_start_date,
                                  "date" = date,
                                  "country" = country,
                                  "surging" = surging,
                                  "variants" = variant_timings),
                    output_options = list(pandoc_args = c(paste0("--metadata=title:",country," COVID-19 report "))))

} else {
  file.create("index.html", "index.pdf", "index.md", "summary_df.rds")
}

if(!fit_excess | !fit_reported){
  #ESFT loop for countries with no deaths,

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

  #catch if no vaccines
  if(
    length(vacc_inputs$primary_doses) == 0
  ){
    vacc_inputs$primary_doses <- 0
  }
  if(
    length(vacc_inputs$booster_doses) == 0
  ){
    vacc_inputs$booster_doses <- 0
  }
  second_dose_delay <- 60
  # Scenarios with capacity constraints
  # ---------------------------------------------------------------------------

  forwards_projection_esft <- squire.page:::run_booster(
    country = country,
    R0 = Rt + ((R0 - Rt)/2),
    time_period = time_period_esft,
    seeding_cases = seeding_cases_esft,
    init = init,
    primary_doses =  as.integer(mean(tail(vacc_inputs$primary_doses,7))),
    second_dose_delay = second_dose_delay,
    booster_doses =  as.integer(mean(tail(vacc_inputs$booster_doses,7))),
    vaccine_coverage_mat = vacc_inputs$vaccine_coverage_mat
  )

  # And adjust their time variable so that we have t = 0 as today
  forwards_projection_esft <- lapply(list(forwards_projection_esft), function(x) {
    for(i in seq_len(dim(x$output)[3])) {
      x$output[,"time",i] <- x$output[,"time",i] - 1
    }
    rownames(x$output) <- as.character(seq.Date(date, date + nrow(x$output) -1, 1))
    return(x)
  })[[1]]

  # major summaries
  r_forwards_projection_esft <- lapply(list(forwards_projection_esft), r_list_format, date)[[1]]
  rt_forwards_projection_esft <- lapply(
    list(forwards_projection_esft),
    rt_creation_simple_out_vaccine, date, date+89
  )[[1]]


  cases_milds_esft <- map(list(forwards_projection_esft), function(x){
    x$samples <- list(list())
    x$squire_model <- squire.page::nimue_booster_model()
    x$parameters$replicates <- x$parameters$seed <- x$parameters$use_dde <- NULL
    x$inputs$start_date <- end_date
    summarise_infection_types(x, end_date)
  })[[1]]

  age_cases_milds <- cases_milds_esft$age
  cases_milds_esft <- cases_milds_esft$total

  age_split_esft <- get_age_output(forwards_projection_esft, end_date, end_date) %>%
    left_join(
      age_cases_milds,
      by = c("replicate", "date", "age_group")
    ) %>%
    summarise_age_dependent()
  rm(age_cases_milds)
}

if(document){
  sim_end_date <- as.Date(date)+90

  # summarise the projections
  scenario_list <- list()
  if(fit_excess){
    scenario_list$r_excess_scenario <- r_forwards_projection_excess
  } else {
    scenario_list$r_excess_scenario <- r_forwards_projection_esft
  }
  if(fit_reported){
    scenario_list$r_reported_scenario <- r_forwards_projection_reported
  } else {
    scenario_list$r_reported_scenario <- r_forwards_projection_esft
  }
  data_sum <- map(scenario_list, function(pd){

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

  data_sum$r_excess_scenario$fit_type <- "Excess Mortality"
  data_sum$r_excess_scenario$death_calibrated <- sum(excess_deaths$deaths) > death_limit
  data_sum$r_reported_scenario$fit_type <- "Reported Deaths"
  data_sum$r_reported_scenario$death_calibrated <- sum(reported_deaths$deaths) > death_limit

  rt_list <- list()
  if(fit_excess){
    rt_list$rt_excess_scenario <- rt_forwards_projection_excess
  } else {
    rt_list$rt_excess_scenario <- rt_forwards_projection_esft
  }
  if(fit_reported){
    rt_list$rt_reported_scenario <- rt_forwards_projection_reported
  } else {
    rt_list$rt_reported_scenario <- rt_forwards_projection_esft
  }

  rt_list$rt_excess_scenario$fit_type <- "Excess Mortality"
  rt_list$rt_excess_scenario$death_calibrated <- sum(excess_deaths$deaths) > death_limit
  rt_list$rt_reported_scenario$fit_type <- "Reported Deaths"
  rt_list$rt_reported_scenario$death_calibrated <- sum(reported_deaths$deaths) > death_limit

  # combine and annotate
  data_sum <- do.call(rbind, data_sum)
  data_sum <- rbind(
    data_sum,
    do.call(rbind, rt_list) %>% filter(date <= max(data_sum$date))
  ) %>% arrange(date, fit_type)
  rownames(data_sum) <- NULL

  if(fit_excess){
    age_stratified <- age_split_excess
  } else {
    age_stratified <- age_split_esft %>%
      mutate(fit_type = "Excess Mortality")
  }
  if(fit_reported){
    age_stratified <- rbind(
      age_stratified,
      age_split_reported
    )
  } else {
    age_stratified <- rbind(
      age_stratified,
      age_split_esft %>%
        mutate(fit_type = "Reported Deaths")
    )
  }
  age_stratified$iso3c <- iso3c
  age_stratified$age_group <- paste0("'", as.character(age_stratified$age_group), "'")
  write.csv(age_stratified, "age_stratified.csv", row.names = FALSE, quote = FALSE)

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
  data_sum$version <- "v10"
  data_sum <- dplyr::mutate(data_sum, across(dplyr::starts_with("y_"), ~round(.x,digits = 2)))

  write.csv(data_sum, "projections.csv", row.names = FALSE, quote = FALSE)
} else {
  file.create(c("projections.csv", "age_stratified.csv"))
}

