orderly_id <- tryCatch(orderly::orderly_run_info()$id,
                       error = function(e) "<id>")

print(sessionInfo())
RhpcBLASctl::blas_set_num_threads(1L)
RhpcBLASctl::omp_set_num_threads(1L)

document <- as.logical(document)

# pandoc linking
if(file.exists("Q:\\COVID-Fitting\\pandoc") & document) {
  rmarkdown:::set_pandoc_info("Q:\\COVID-Fitting\\pandoc")
  Sys.setenv(RSTUDIO_PANDOC="Q:\\COVID-Fitting\\pandoc")
  tinytex::use_tinytex("Q:\\COVID-Fitting\\TinyTex")
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
date_0 <- date

## Get the excess mortality estimates from the economist
excess_deaths <- readRDS("excess_deaths.Rds") %>%
  rename(iso = iso3c) %>%
  filter(iso == iso3c)

country <- squire::population$country[match(iso3c, squire::population$iso3c)[1]]

## MAIN LOOP IS ONLY FOR THOSE WITH DEATHS
if(sum(excess_deaths$deaths) > 0) {

  # get the raw data correct
  data <- excess_deaths[, c("week_start", "week_end", "deaths")] %>%
    arrange(week_start) %>%
    rename(date_start = week_start, date_end = week_end)

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

  # conduct unmitigated
  pop <- squire::get_population(country)

  ## -----------------------------------------------------------------------------
  ## Step 2: Fit model
  ## -----------------------------------------------------------------------------

  squire_model <- squire.page:::nimue_booster_model()
  start_date <-  data$date_start[1] - 30

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
      parameters$icu_beds <- parameters$icu_beds / 0.7
    }
  ## -----------------------------------------------------------------------------
  ## Step 2b: Calculate sampled parameters
  ## -----------------------------------------------------------------------------

  #load vaccine inputs
  vacc_inputs <- get_vaccine_inputs(iso3c)
  variant_characteristics <- readRDS("variant_characteristics.rds")[[iso3c]]
  #check if omicron starts during delta
  if(variant_characteristics$Omicron$start_date <
    variant_characteristics$Delta$start_date +
    variant_characteristics$Delta$shift_duration) {
    variant_characteristics$Omicron$shift_duration <-
      as.numeric(
        variant_characteristics$Omicron$start_date + variant_characteristics$Omicron$shift_duration
      ) -
      as.numeric(
        variant_characteristics$Delta$start_date + variant_characteristics$Delta$shift_duration
      ) - 1
    variant_characteristics$Omicron$start_date <- variant_characteristics$Delta$start_date +
      variant_characteristics$Delta$shift_duration + 1
  }

  samples <- 100
  distribution <- map(seq_len(samples), function(x){
    pars <- list()
    ## disease parameters:
    #duration of recovery
    pars$dur_R <- rpois(1, 365)
    #add the others

    ## vaccine parameters, should probably refit this each time:
    pars$dur_V <- map_dbl(vacc_inputs$dur_V, ~rpois(1, .x))
    pars$rel_infectiousness_vaccinated <- rbeta(1, 1, 1)
    #could simulate some level of under/over reporting on vaccine usage
    tt_vaccine_change <- as.numeric(vacc_inputs$date_vaccine_change - start_date)
    if(tt_vaccine_change[1] > 0){
      tt_vaccine_change <- c(0, tt_vaccine_change)
    } else {
      tt_vaccine_change <- c(tt_vaccine_change[1] - 1, tt_vaccine_change)
    }
    pars$first_doses <- vacc_inputs$first_doses
    pars$tt_first_doses <- tt_vaccine_change
    pars$second_doses <- vacc_inputs$second_doses
    pars$tt_second_doses <- tt_vaccine_change
    pars$booster_doses <- vacc_inputs$booster_doses
    pars$tt_booster_doses <- tt_vaccine_change
    pars$vaccine_coverage_mat <- vacc_inputs$vaccine_coverage_mat

    ## variant parameters:
    #add some kind of sampling
    #calculate changes
    dur_R_change <- variant_immune_escape(variant_characteristics,
                                          pars$dur_R, start_date)
    prob_hosp_multiplier <- variant_changes_over_time(variant_characteristics,
                                                    "prob_hosp_multiplier", start_date)
    prob_severe_multiplier <- variant_changes_over_time(variant_characteristics,
                                                      "prob_severe_multiplier", start_date)
    dur_ICU <- variant_changes_over_time(variant_characteristics,
                                       "dur_ICU", start_date)
    dur_ICU_death <- variant_changes_over_time(variant_characteristics,
                                             "dur_ICU_death", start_date)
    dur_hosp <- variant_changes_over_time(variant_characteristics,
                                        "dur_hosp", start_date)
    dur_hosp_death <- variant_changes_over_time(variant_characteristics,
                                              "dur_hosp_death", start_date)
    pars$dur_R <- dur_R_change$var
    pars$tt_dur_R <- dur_R_change$tt
    pars$prob_hosp_multiplier <- prob_hosp_multiplier$var
    pars$tt_prob_hosp_multiplier <- prob_hosp_multiplier$tt
    pars$prob_severe_multiplier <- prob_severe_multiplier$var
    pars$tt_prob_severe_multiplier <- prob_severe_multiplier$tt
    pars$dur_get_mv_survive <-  dur_ICU$var
    pars$tt_dur_get_mv_survive <- dur_ICU$tt
    pars$dur_get_mv_die <-  dur_ICU_death$var
    pars$tt_dur_get_mv_die <- dur_ICU_death$tt
    pars$dur_get_ox_survive <-  dur_hosp$var
    pars$tt_dur_get_ox_survive <- dur_hosp$tt
    pars$dur_get_ox_die <-  dur_hosp_death$var
    pars$tt_dur_get_ox_die <- dur_hosp_death$tt

    vacc_inputs <- vaccine_eff_over_time(vacc_inputs, variant_characteristics, start_date)
    pars$dur_V <- vacc_inputs$dur_V
    pars$tt_dur_V <- vacc_inputs$tt_vaccine_efficacy_change
    pars$vaccine_efficacy_infection <- vacc_inputs$vaccine_efficacy_infection
    pars$tt_vaccine_efficacy_infection <- vacc_inputs$tt_vaccine_efficacy_change
    pars$vaccine_efficacy_disease <- vacc_inputs$vaccine_efficacy_disease
    pars$tt_vaccine_efficacy_disease <- vacc_inputs$tt_vaccine_efficacy_change

    pars
  })

  ## -----------------------------------------------------------------------------
  ## Step 2c: Actual run
  ## -----------------------------------------------------------------------------

  #use parallel if asked for
  if(parallel){
    future::plan(future::multisession())
  }

  progressr::handlers(global = TRUE)

  out <- rt_optimise(
    data = data,
    distribution = distribution[1:3],
    squire_model = squire_model,
    parameters = parameters,
    start_date = start_date,
    parallel = parallel,
    rt_spacing = 14,
    initial_r = 3,
    initial_infections = 100,
    proposal_width = 0.9,
    n_particles = n_particles,
    k = 14
  ) #2 hours!!

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

  rtp <- rt_plot_immunity(out, R0_plot = TRUE)
  rtp2 <- rt_plot_immunity(out, R0_plot = FALSE)

  date_range <- as.Date(c(min(as.Date(out$replicate_parameters$start_date)),date_0))
  if(pars_obs$cases_fitting){
    suppressMessages(suppressWarnings(
      bottom <- cowplot::plot_grid(
        intervention + scale_x_date(date_breaks = "1 month", date_labels = "%b" ,limits = date_range),
        d + scale_x_date(date_breaks = "1 month", date_labels = "%b" ,limits = date_range),
        drjacoby::plot_mc_acceptance(out$pmcmc_results$drjacoby_out),
        infection_adjustment_plot,
        ncol=1,
        rel_heights = c(0.4,0.6, 0.6, 0.4))
    ))
  } else {
    suppressMessages(suppressWarnings(
      bottom <- cowplot::plot_grid(
        intervention + scale_x_date(date_breaks = "1 month", date_labels = "%b" ,limits = date_range),
        d + scale_x_date(date_breaks = "1 month", date_labels = "%b" ,limits = date_range),
        drjacoby::plot_mc_acceptance(out$pmcmc_results$drjacoby_out),
        ncol=1,
        rel_heights = c(0.4,0.6,0.6))
    ))

  }

  combined <- cowplot::plot_grid(header,
                                 cowplot::plot_grid(g1, line_v, bottom, ncol = 3, rel_widths = c(3,0.1,1)),
                                 ncol=1, rel_heights = c(1,15))



  ggsave("fitting.pdf",width=24, height=12,
         combined)


  ## Save the grid out object

  #remove all dr jacoby data, no reason we should need it

  out$pmcmc_results$drjacoby_out <- NULL

  ## now let's trim the out for saving really small
  output_temp <- out$output
  out$output <- NULL
  saveRDS(out, "grid_out.rds")
  #reattach output
  out$output <- output_temp
  rm(output_temp)


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
    model_user_args <- extend_vaccine_inputs(vacc_inputs, time_period, out, date_0 = date_0)

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
    icu <- nimue_format(maintain_scenario, "ICU_demand", date_0 = date_0)
    icu <- icu[icu$compartment == "ICU_demand",]
    icu_28 <- group_by(icu, t) %>%
      summarise(i_tot = mean(y, na.rm = TRUE),
                i_min = t_test_safe(y)$conf.int[1],
                i_max = t_test_safe(y)$conf.int[2])

    hosp <- nimue_format(maintain_scenario, "hospital_demand", date_0 = date_0)
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

    # Get the optimistic/pessimistic scenarios
    Rt_futures <- get_future_Rt(out, forcast_days = 28)

    #update the model_user_args to set next Rt to our values
    model_user_args_optim <- update_Rt(model_user_args, Rt_futures$optimistic, out)

    mitigation_scenario <- squire::projections(out,
                                               R0_change = c(1),
                                               tt_R0 = c(0),
                                               time_period = time_period,
                                               model_user_args = model_user_args_optim)

    r_mitigation_scenario <- r_list_format(mitigation_scenario, date_0)
    rt_mitigation_scenario <- rt_creation_vaccine(mitigation_scenario, date_0, date_0+89)

    rm(mitigation_scenario)
    rm(model_user_args_optim)

    # Relax by 50% for 3 months
    model_user_args_pess <- update_Rt(model_user_args, Rt_futures$pessimistic, out)
    reverse_scenario <- squire::projections(out,
                                            R0_change = c(1),
                                            tt_R0 = c(0),
                                            time_period = time_period,
                                            model_user_args = model_user_args_pess)
    r_reverse_scenario <- r_list_format(reverse_scenario, date_0)
    rt_reverse_scenario <- rt_creation_vaccine(reverse_scenario, date_0, date_0+89)
    rm(reverse_scenario)
    rm(model_user_args_pess)

    #legacy scenarios
    mitigation_scenario_leg <- squire::projections(out,
                                                   R0_change = c(0.5),
                                                   tt_R0 = c(0),
                                                   time_period = time_period,
                                                   model_user_args = model_user_args)
    r_mitigation_scenario_leg <- r_list_format(mitigation_scenario_leg, date_0)
    rt_mitigation_scenario_leg <- rt_creation_vaccine(mitigation_scenario_leg, date_0, date_0+89)
    rm(mitigation_scenario_leg)


    reverse_scenario_leg <- squire::projections(out,
                                                R0_change = c(1.5),
                                                tt_R0 = c(0),
                                                time_period = time_period,
                                                model_user_args = model_user_args)
    r_reverse_scenario_leg <- r_list_format(reverse_scenario_leg, date_0)
    rt_reverse_scenario_leg <- rt_creation_vaccine(reverse_scenario_leg, date_0, date_0+89)
    rm(reverse_scenario_leg)

    ## -----------------------------------------------------------------------------
    ## 4.2. Investigating a capacity surge
    ## -----------------------------------------------------------------------------

    # update for surged healthcare
    out_surged <- out
    rm(out)

    out_surged$pmcmc_results$inputs$model_params$hosp_beds <- 1e10
    out_surged$pmcmc_results$inputs$model_params$ICU_beds <- 1e10

    # out_surged <- generate_draws_pmcmc(pmcmc = pmcmc,
    #                                    burnin = ceiling(n_mcmc/10),
    #                                    n_chains = n_chains,
    #                                    squire_model = pmcmc$inputs$squire_model,
    #                                    replicates = replicates,
    #                                    n_particles = n_particles,
    #                                    forecast = 0,
    #                                    country = country,
    #                                    population = squire::get_population(iso3c = iso3c)$n,
    #                                    interventions = out_interventions,
    #                                    data = pmcmc$inputs$data)
    out_surged <- generate_draws(out_surged,
                                 draws = replicates,
                                 pars.list = NULL)

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

    Rt_futures_surged <- get_future_Rt(out_surged, forcast_days = 28)

    model_user_args_optim <- update_Rt(model_user_args, Rt_futures_surged$optimistic, out_surged)

    mitigation_scenario_surged <- squire::projections(out_surged,
                                                      R0_change = c(1),
                                                      tt_R0 = c(0),
                                                      time_period = time_period,
                                                      model_user_args = model_user_args_optim)
    r_mitigation_scenario_surged <- r_list_format(mitigation_scenario_surged, date_0)
    rt_mitigation_scenario_surged <- rt_creation_vaccine(mitigation_scenario_surged, date_0, date_0+89)
    rm(mitigation_scenario_surged)
    rm(model_user_args_optim)

    model_user_args_pess <- update_Rt(model_user_args, Rt_futures_surged$pessimistic, out_surged)

    reverse_scenario_surged <- squire::projections(out_surged,
                                                   R0_change = c(1),
                                                   tt_R0 = c(0),
                                                   time_period = time_period,
                                                   model_user_args = model_user_args_pess)
    r_reverse_scenario_surged <- r_list_format(reverse_scenario_surged, date_0)
    rt_reverse_scenario_surged <- rt_creation_vaccine(reverse_scenario_surged, date_0, date_0+89)
    rm(reverse_scenario_surged)
    rm(model_user_args_pess)

    #legacy scenarios
    mitigation_scenario_surged_leg <- squire::projections(out_surged,
                                                          R0_change = c(0.5),
                                                          tt_R0 = c(0),
                                                          time_period = time_period,
                                                          model_user_args = model_user_args)
    r_mitigation_scenario_surged_leg <- r_list_format(mitigation_scenario_surged_leg, date_0)
    rt_mitigation_scenario_surged_leg <- rt_creation_vaccine(mitigation_scenario_surged_leg, date_0, date_0+89)
    rm(mitigation_scenario_surged_leg)


    reverse_scenario_surged_leg <- squire::projections(out_surged,
                                                       R0_change = c(1.5),
                                                       tt_R0 = c(0),
                                                       time_period = time_period,
                                                       model_user_args = model_user_args)
    r_reverse_scenario_surged_leg <- r_list_format(reverse_scenario_surged_leg, date_0)
    rt_reverse_scenario_surged_leg <- rt_creation_vaccine(reverse_scenario_surged_leg, date_0, date_0+89)
    rm(reverse_scenario_surged_leg)

    rm(out_surged)

    o_list <- named_list(
      r_maintain_scenario,
      r_mitigation_scenario,
      r_reverse_scenario,
      r_maintain_scenario_surged,
      r_mitigation_scenario_surged,
      r_reverse_scenario_surged,
      r_mitigation_scenario_leg,
      r_reverse_scenario_leg,
      r_mitigation_scenario_surged_leg,
      r_reverse_scenario_surged_leg
    )

    rt_list <- named_list(
      rt_maintain_scenario,
      rt_mitigation_scenario,
      rt_reverse_scenario,
      rt_maintain_scenario_surged,
      rt_mitigation_scenario_surged,
      rt_reverse_scenario_surged,
      rt_mitigation_scenario_leg,
      rt_reverse_scenario_leg,
      rt_mitigation_scenario_surged_leg,
      rt_reverse_scenario_surged_leg
    )

    rt_futures_df <- rbind(
      summarise_rt_futures(Rt_futures) %>%
        mutate(surged = FALSE),
      summarise_rt_futures(Rt_futures_surged) %>%
        mutate(surged = TRUE)
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

    #delta stuff
    adjust_delta_index <- variant_characteristics

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
                                    "rt" = rtp2,
                                    "adjust_delta" = adjust_delta_index,
                                    "rt_futures" = rt_futures_df),
                      output_options = list(pandoc_args = c(paste0("--metadata=title:",country," COVID-19 report "))))

    } else {
    # create empty files
    file.create("index.html", "index.pdf", "index.md",
                "summary_df.rds", "input_params.json")
  }
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

  vacc_inputs <- get_vaccine_inputs(iso3c)

  # set up vaccine inits
  vaccine_fitting_flag <- TRUE
  squire_model <- squire.page::nimue_booster_model()
  init <- init_state_nimue(deaths_removed = 0, iso3c,
                           seeding_cases = seeding_cases_esft,
                           vaccinated_already = sum(vacc_inputs$max_vaccine))

  # Scenarios with capacity constraints
  # ---------------------------------------------------------------------------
  reverse_scenario <- squire.page:::run_booster(
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

  mitigation_scenario <- squire.page:::run_booster(
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

  maintain_scenario <- squire.page:::run_booster(
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
  saveRDS(maintain_scenario, "grid_out.rds")

  # Scenarios without capacity constraints
  # ---------------------------------------------------------------------------
  reverse_scenario_surged <- squire.page:::run_booster(
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

  mitigation_scenario_surged <- squire.page:::run_booster(
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

  maintain_scenario_surged <- squire.page:::run_booster(
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
              "fitting.pdf")

  # major summaries
o_list <- lapply(r_list_pass, r_list_format, date_0)

rt_list <- lapply(r_list_pass, rt_creation_vaccine, date_0, date_0+89)

}

if(document){
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
  data_sum[[2]]$scenario <- "Optimistic"
  data_sum[[3]]$scenario <- "Pessimistic"
  data_sum[[4]]$scenario <- "Surged Maintain Status Quo"
  data_sum[[5]]$scenario <- "Surged Optimistic"
  data_sum[[6]]$scenario <- "Surged Pessimistic"
  if(sum(ecdc_df$deaths) > 0){
    #these won't be present for the non-fitting countries
    data_sum[[7]]$scenario <- "Additional 50% Reduction"
    data_sum[[8]]$scenario <- "Relax Interventions 50%"
    data_sum[[9]]$scenario <- "Surged Additional 50% Reduction"
    data_sum[[10]]$scenario <- "Surged Relax Interventions 50%"
  } else {
    data_sum[[7]] <- data_sum[[2]]
    data_sum[[7]]$scenario <- "Additional 50% Reduction"
    data_sum[[8]] <- data_sum[[3]]
    data_sum[[8]]$scenario <- "Relax Interventions 50%"
    data_sum[[9]] <- data_sum[[5]]
    data_sum[[9]]$scenario <- "Surged Additional 50% Reduction"
    data_sum[[10]] <- data_sum[[6]]
    data_sum[[10]]$scenario <- "Surged Relax Interventions 50%"
  }

  rt_list[[1]]$scenario <- "Maintain Status Quo"
  rt_list[[2]]$scenario <- "Optimistic"
  rt_list[[3]]$scenario <- "Pessimistic"
  rt_list[[4]]$scenario <- "Surged Maintain Status Quo"
  rt_list[[5]]$scenario <- "Surged Optimistic"
  rt_list[[6]]$scenario <- "Surged Pessimistic"
  if(sum(ecdc_df$deaths) > 0){
    #these won't be present for the non-fitting countries
    rt_list[[7]]$scenario <- "Additional 50% Reduction"
    rt_list[[8]]$scenario <- "Relax Interventions 50%"
    rt_list[[9]]$scenario <- "Surged Additional 50% Reduction"
    rt_list[[10]]$scenario <- "Surged Relax Interventions 50%"
  } else {
    rt_list[[7]] <- rt_list[[2]]
    rt_list[[7]]$scenario <- "Additional 50% Reduction"
    rt_list[[8]] <- rt_list[[3]]
    rt_list[[8]]$scenario <- "Relax Interventions 50%"
    rt_list[[9]] <- rt_list[[5]]
    rt_list[[9]]$scenario <- "Surged Additional 50% Reduction"
    rt_list[[10]] <- rt_list[[6]]
    rt_list[[10]]$scenario <- "Surged Relax Interventions 50%"
  }

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
} else {
  file.create("projections.csv")
}
