r_list_format <- function(out, date_0) {

  df <- nimue_format(out,
                      var_select = c("infections","deaths","hospital_demand",
                                     "ICU_demand", "D", "hospital_incidence","ICU_incidence"),
                      date_0 = date_0)

  pr <- nimue_format(out, var_select = c("S","R","D"), date_0 = date_0) %>%
    na.omit %>%
    pivot_wider(names_from = compartment, values_from = y) %>%
    mutate(y = sum(squire.page:::get_parameters(out)$population)-D-R-S,
           compartment = "prevalence") %>%
    select(replicate, compartment, t, y, date)

  rbind(df, pr)

}
named_list <- function(...) {
  get <- as.character(match.call())
  l <- list(...)
  names(l) <- utils::tail(get, -1)
  return(l)
}

summarise_rt_futures <- function(Rt_futures){
  #for each trend type
  map_dfr(Rt_futures, function(trend){
    #for each replicate
    values <- map_dbl(trend, ~tail(.x$R0, 1))
    tibble(
      q_025 = quantile(values, 0.025),
      med = quantile(values, 0.5),
      q_975 = quantile(values, 0.975)
    )
  }, .id = "scenario")
}


multiplier_changes_over_time <- function(variant_timings, variable, x, start_date) {
  variant_timings <- arrange(variant_timings, start_date)

  #adjust for historic changes
  var <- map_dbl(variable, ~.x[[x]])[variant_timings$variant]
  var <- var/lag(var, 1, 1)
  var <- c(Wild = 1, var)

  variant_par_changes_over_time(variant_timings, as.list(var), start_date)
}

variant_par_changes_over_time <- function(variant_timings, variable, start_date) {
  change_tt <- map(transpose(variant_timings), function(l){
    period <- as.numeric(as_date(c(l$start_date, l$end_date)) - start_date)
    tt <- seq(period[1], period[2], by = 1)
    list(
      change = seq(0, 1, length.out = length(tt) + 1)[-1],
      tt = tt
    )
  })
  tt <- map(change_tt, ~.x$tt) %>% unlist()
  var <- map(seq_along(change_tt), function(i){
    if(i == 1){
      start_value <- variable[["Wild"]]
    } else {
      start_value <- variable[[variant_timings$variant[i - 1]]]
    }
    end_value <- variable[[variant_timings$variant[i]]]
    start_value * (1 - change_tt[[i]]$change) + end_value * change_tt[[i]]$change
  }) %>% unlist()

  #add the zero entry if needed
  if(all(tt > 0)){
    tt <- c(0, tt)
    var <- c(variable[["Wild"]], var)
  }

  return(
    list(
      var = var,
      tt = tt
    )
  )
}

variant_immune_escape <- function(variant_timings, immune_escape, x, dur_R, start_date) {
  variant_timings <- arrange(variant_timings, start_date)
  #adjust immune escape for previous escape levels
  immune_escapes <- map_dbl(immune_escape, ~.x[x])
  immune_escapes <- (immune_escapes - lag(immune_escapes, 1, default = 0))/(
    1 - lag(immune_escapes, 1, default = 0)
  )
  immune_escapes <- immune_escapes[variant_timings$variant]

  dur_Rs <- map(transpose(variant_timings), function(l){
    shift_duration <- as.numeric(l$end_date - l$start_date)
    c(1 / (
      (shift_duration / dur_R[x] - log(1 - immune_escape[[l$variant]][[x]])) /
        shift_duration
    ), dur_R[x])
  }) %>% unlist()

  tt <- map(transpose(variant_timings), function(l){
    as.numeric(as_date(c(l$start_date, l$end_date)) - start_date)
  }) %>% unlist()

  if(all(tt > 0)){
    tt <- c(0, tt)
    dur_Rs <- c(dur_R[x], dur_Rs)
  }

  return(
    list(
      var = dur_Rs,
      tt = tt
    )
  )
}

get_future_Rt_optimised <- function(model_out, forcast_days){
  if("excess_nimue_simulation" %in% class(model_out)){
    stop("This function cannot be used with excess_nimue_simulation at the moment")
  }
  Rt_futures <- purrr::map(model_out$samples, function(sample){
    #look at the changes over length of time of the forcast period
    Rt_changes <- map_dbl(seq(sample$tt_R0[2], tail(sample$tt_R0, 1) - forcast_days), function(t){
      squire.page:::block_interpolate(t + forcast_days, sample$R0, sample$tt_R0)/
        squire.page:::block_interpolate(t, sample$R0, sample$tt_R0)
    })
    #calculate the 20, 50 and 80% quantiles to get the optimisitc, central, and
    #pessimistic changes
    changes <- as.numeric(quantile(Rt_changes, c(0.2, 0.5, 0.8)))
    #expand these out to cover the forcast period and reach the final value by
    #this time
    change_time <- mean(diff(sample$tt_R0[-1]))
    map(changes, function(change){
      tt_R0 <- seq(0, forcast_days, by = 1)
      R0 <- seq(tail(sample$R0, 1), tail(sample$R0, 1)*change, length.out = length(tt_R0) + 1)[-1]
      list(
        R0 = R0,
        tt_R0 = tt_R0
      )
    })
  })
  #reformat
  new_names <- c(1, 2, 3)
  names(new_names) <- c("optimistic", "central", "pessimistic")
  map(new_names, function(x) map(Rt_futures, ~.x[[x]]))
}

update_Rt_optimised <- function(model_user_args, Rt_future){
  map(seq_along(Rt_future), function(sample_index){
    c(model_user_args[[sample_index]], Rt_future[[sample_index]])
  })
}

summarise_fit <- function(file, out, country, iso3c, end_date, start_date){
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
      out, out$inputs$data, date = end_date,
      date_0 = out$inputs$start_date, forecast = 0,
      single = TRUE, full = TRUE) +
      theme(legend.position = "none")))

  suppressMessages(suppressWarnings(
    cas_plot <- cases_plot_single(
      df = nimue_format(out, "infections", date_0 = out$inputs$start_date),
      data = out$inputs$data,
      date = end_date,
      date_0 = out$inputs$start_date) +
      theme(legend.position = "none")))

  rtp <- rt_plot_immunity(out, R0_plot = TRUE)
  rtp2 <- rt_plot_immunity(out, R0_plot = FALSE)

  date_range <- as.Date(c(out$inputs$start_date, end_date))

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

  ggsave(file,width=24, height=12,
         combined)
}

trim_output <- function(out, trimming){
  out_trim <- squire.page::trim_rt_optimise(out, trimming)

  if(length(out_trim$samples) == 0){
    warning("No suitable trajectories calculated")
    out
  } else {
    out_trim
  }
}

save_output <- function(out, file){
  ## now let's trim the output for saving
  out$output <- NULL
  if(inherits(out, "rt_optimised_trimmed")){
    out$excluded$output <- NULL
  }
  saveRDS(out, file)
}

extend_output <- function(out, vacc_inputs, time_period, start_date, end_date){
  ## We need to know work out vaccine doses and efficacy going forwards
  #catch for no vaccines
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
  model_user_args <- extend_vaccine_inputs(vacc_inputs, time_period, out, end_date = end_date)
  #For now these are all identicall we so'll add them to $parameters.
  #I've left this framework in case we want to add randomness to the doses
  out_extended <- out
  out_extended$parameters$protection_delay_time <- out_extended$parameters$protection_delay_time + time_period
  for(x in names(model_user_args[[1]])){
    if(str_detect(x, "tt")){
      out_extended$parameters[[x]] <- c(out$parameters[[x]], as.numeric(end_date - start_date) + 1)
    } else {
      out_extended$parameters[[x]] <- c(out$parameters[[x]], model_user_args[[1]][[x]])
    }
  }
  out_extended
}

check_breached_capacity <- function(forwards_projection, country, start_date){
  ## Is capacity passed prior to today
  icu_cap <- squire:::get_ICU_bed_capacity(country)
  hosp_cap <- squire:::get_hosp_bed_capacity(country)
  t_test_safe <- function(x, ...) {
    out <- try(t.test(x, ...), silent = TRUE)
    if (inherits(out, "try-error")) {
      out <- list("conf.int"=c(mean(x),mean(x)))
    }
    return(out)
  }
  # has surging been required
  icu <- nimue_format(forwards_projection, "ICU_demand", date_0 = start_date)
  icu <- icu[icu$compartment == "ICU_demand",]
  icu_28 <- group_by(icu, t) %>%
    summarise(i_tot = mean(y, na.rm = TRUE),
              i_min = t_test_safe(y)$conf.int[1],
              i_max = t_test_safe(y)$conf.int[2])
  hosp <- nimue_format(forwards_projection, "hospital_demand", date_0 = start_date)
  hosp <- hosp[hosp$compartment == "hospital_demand",]
  hosp_28 <- group_by(hosp, t) %>%
    summarise(i_tot = mean(y, na.rm = TRUE),
              i_min = t_test_safe(y)$conf.int[1],
              i_max = t_test_safe(y)$conf.int[2])

  any(na.omit(icu_28$i_tot > icu_cap)) || any(na.omit(icu_28$i_tot > icu_cap))
}

get_compartment_array <- function(output, compartment, end_date, replicate){
  comp <- output[
    lubridate::as_date(rownames(output)) > end_date - 30*3,
    str_detect(colnames(output), compartment),
    replicate
  ]
  #convert to 3d (date, age, vaccine status)
  array_ <- array(NA, c(nrow(comp), 17, 7))
  for(row in seq_len(17)){
    array_[,row,] <- comp[, stringr::str_detect(colnames(comp), paste0(compartment,"\\[", row, ","))]
  }
  rownames(array_) <- rownames(comp)
  array_
}

summarise_infection_types <- function(forwards_projection, end_date){
  df <- purrr::map_dfr(seq_along(forwards_projection$samples), function(replicate){
    params <- squire.page:::setup_parameters(forwards_projection$squire_model,
                                             c(
                                               forwards_projection$parameters,
                                               forwards_projection$samples[[replicate]]
                                             )
    )
    E2_array <- get_compartment_array(forwards_projection$output, "E2", end_date, replicate)
    dates <- lubridate::as_date(rownames(E2_array))
    #get just the transitions
    E2_array <- E2_array * params$gamma_E
    #get the time/vaccine varying proability of hopsitalisation
    tt <- as.numeric(dates - forwards_projection$inputs$start_date)
    prob_hosp_multiplier <- squire.page:::block_interpolate(tt, params$prob_hosp_multiplier, params$tt_prob_hosp_multiplier)
    prob_hosp <- params$prob_hosp[
      squire.page:::block_interpolate(tt, seq_len(dim(params$prob_hosp)[1]), params$tt_vaccine_efficacy_disease), ,
    ] * prob_hosp_multiplier
    new_Case <- E2_array * prob_hosp
    new_Mild <- E2_array * (1 - prob_hosp)
    #get hospitalisations
    ICase2_array <- get_compartment_array(forwards_projection$output, "ICase2", end_date, replicate)
    new_hosp <- params$gamma_ICase * ICase2_array
    #sum over ages and vaccine status
    tibble(
      replicate = replicate,
      date = unlist(map(dates, ~rep(.x, 17))),
      age_group = rep(ordered(c(
        "0-5",
        "5-10",
        "10-15",
        "15-20",
        "20-25",
        "25-30",
        "30-35",
        "35-40",
        "40-45",
        "45-50",
        "50-55",
        "55-60",
        "60-65",
        "65-70",
        "70-75",
        "75-80",
        "80+"
      )), length(dates)),
      new_Case = unlist(array_tree(rowSums(new_Case, dims = 2))),
      new_Mild = unlist(array_tree(rowSums(new_Mild, dims = 2))),
      new_hosp = unlist(array_tree(rowSums(new_hosp, dims = 2)))
    )
  })
  list(
    total = df %>%
      group_by(replicate, date) %>%
      summarise(
        new_Case = sum(new_Case),
        new_Mild = sum(new_Mild),
        .groups = "drop"
      ) %>%
      mutate(
        date = lubridate::as_date(date),
        projection = date > end_date
      ),
    age = df %>%
      mutate(
        date = lubridate::as_date(date),
        projection = date > end_date
      )
  )
}

detected_cases_estimation <- function(cases_milds, cases){
    #could split this further into ICU, HOSP, Mild
    #should probably add tail cases fitting back into to get this closer
    l_post <- function(detected_infections, p_case_eff, p_mild, Cases, Mild){
      p_case <- p_mild + (1 - p_mild) * p_case_eff
      #prior that pcase is about 10x higher than pmild
      lp <- dgamma(p_case/p_mild, shape = 100, rate = 10, log = TRUE)
      detected_infections[detected_infections > Mild + Cases] <- (Mild + Cases)[detected_infections > Mild + Cases]
      prop_mild <- Mild * p_mild /(Mild * p_mild + Cases * p_case) #expectation of this, since we're doing MLE doesn't matter
      detected_mild <- (detected_infections*prop_mild) %>% round()
      ll <- (dbinom(detected_mild, Mild, p_mild, log = TRUE) +
               dbinom(detected_infections - detected_mild, Cases, p_case, log = TRUE)) %>%
        sum()
      if(ll < -10^9){
        ll <- -10^9
      }
      lp + ll
    }
    p_detection <- cases_milds %>%
      left_join(cases, by = c("date")) %>%
      filter(!is.na(detected_infections)) %>%
      mutate(
        across(
          c(new_Case, new_Mild, detected_infections),
          ~as.integer(round(.x))
        )
      ) %>%
      group_by(replicate) %>%
      mutate(
        initial_value = min(c(1, sum(detected_infections)/sum(new_Case + new_Mild))),
      ) %>%
      summarise(
        #non optimial, often barely moves from initial values when detection is very low
        out = list(optim(c(initial_value, initial_value), function(x, detected_infections, new_Case, new_Mild, l_post){
          l_post(detected_infections, x[1], x[2], new_Case, new_Mild)
        }, method = "L-BFGS-B", lower = rep(10^-10, 2), upper = rep(0.9999999, 2),
        control = list(fnscale = -1),#, lmm = 10, factr = 1e10, pgtol = 1e2),
        detected_infections = detected_infections,
        new_Case = new_Case, new_Mild = new_Mild, l_post = l_post)$par),
        .groups = "drop",
        initial_value = initial_value[1]
      ) %>%
      transmute(
        replicate = replicate,
        p_mild = purrr::map_dbl(out, ~.x[2]),
        p_case = p_mild + (1 - p_mild) * purrr::map_dbl(out, ~.x[1])
      )

    cases_milds <- cases_milds %>%
      left_join(cases, by = c("date")) %>%
      left_join(p_detection, by = "replicate") %>%
      mutate(
        new_detected_Case = new_Case * p_case,
        new_detected_Mild = new_Mild * p_mild
      )

    # temp <- cases_milds %>%
    #   filter(!projection) %>%
    #   group_by(date) %>%
    #   summarise(
    #     est_detected_025 = quantile(new_detected_Case + new_detected_Mild, c(0.025)),
    #     est_detected_50 = quantile(new_detected_Case + new_detected_Mild, c(0.5)),
    #     est_detected_975 = quantile(new_detected_Case + new_detected_Mild, c(0.975))
    #   )
    # ggplot() +
    #   geom_line(data = cases_milds %>% filter(!projection & replicate == 1),
    #             aes(x = date, y = detected_infections), linetype = "dashed") +
    #   geom_ribbon(data = temp, aes(x = date, ymin = est_detected_025 , ymax = est_detected_975),
    #               alpha = 0.2) +
    #   geom_line(data = temp, aes(x = date, y = est_detected_50))


  cases_milds %>%
      select(!detected_infections)
}

update_parameters <- function(parameters, difference_in_t, squire_model){
  tt_pars <- stringr::str_subset(names(parameters), "tt_")
  for(var in tt_pars){
    parameters[[var]] <- parameters[[var]] - difference_in_t
    parameters[[var]][1] <- 0
    if(sum(parameters[[var]] < 0) > 0){
      if(stringr::str_detect(var, "(doses|max_vaccine)")){
        parameters$prefit_vaccines <- TRUE
      }
      parameters[[var]][1] <- parameters[[var]][2] - 1
    }
  }
  if(!is.null(parameters$protection_delay_time)){
    parameters$protection_delay_time <- parameters$protection_delay_time - difference_in_t
  }
  parameters
}

update_distribution <- function(distribution, difference_in_t){
  purrr::map(distribution, function(parameters){
    pars <- update_parameters(parameters, difference_in_t)
  })
}

prefit_vaccines <- function(parameters, distribution, squire_model){
  warning("When pre-fitting vaccines, the compartments are hardcoded and so may need changing in future, if the underlying model changes structure!")
  map(seq_along(distribution), function(i){
    i_parameters <- c(parameters, distribution[[i]])
    i_parameters$prefit_vaccines <- NULL
    tt_pars <- stringr::str_subset(names(i_parameters), "tt_")
    tt_pars_vaccine <- stringr::str_subset(tt_pars, "(doses|max_vaccine)")
    #fix leading zeros
    first_t <- min(unlist(i_parameters[tt_pars_vaccine]), na.rm = TRUE) - 1
    i_parameters <- imap(i_parameters, function(x, par){
      if(par %in% tt_pars){
        x[1] <- first_t
        x <- x - first_t
      }
      x
    })

    odin_pars <- squire.page:::setup_parameters(squire_model, i_parameters)
    odin_pars$beta_set <- 0
    odin_pars$tt_beta <- first_t
    odin_pars$S_0 <- odin_pars$S_0 + odin_pars$E1_0
    odin_pars$E1_0[,1] <- rep(0, length(odin_pars$E1_0[,1]))

    model_instance <- squire_model$odin_model(odin_pars)
    output <- model_instance$run(c(0, -first_t))[2,]

    #extract the compartment values
    S_array <- array(NA, c(17, 7))
    for(age in seq_len(17)){
      S_array[age, ] <- output[paste0("S[", age, ",", seq_len(7), "]")]
    }

    new_dist <- distribution[[i]]
    new_dist$S_0 <- S_array
    new_dist
  })
}

get_age_output <- function(out, date_0, end_date){
  nimue_format(out,
                     var_select = c("infections", "hospital_demand", "ICU_demand", "hospital_incidence", "ICU_incidence"),
                     date_0 = date_0, reduce_age = FALSE) %>%
    select(replicate, date, compartment, age_group, y) %>%
    pivot_wider(names_from = compartment, values_from = y) %>%
    mutate(age_group = ordered(as.character(age_group))) %>%
    filter(date >= end_date)
}

estimate_cases_age <- function(age_cases_milds, cases_milds, end_date){
  age_cases_milds %>%
    left_join(
      cases_milds %>%
        group_by(replicate) %>%
        summarise(
          p_mild = p_mild[1],
          p_case = p_case[1],
        ),
      by = "replicate"
    ) %>%
    mutate(
      new_detected_Case = new_Case * p_case,
      new_detected_Mild = new_Mild * p_mild,
    ) %>%
    select(!c(p_mild, p_case, projection))%>%
    filter(date >= end_date) %>%
    rename(non_hospitalised_cases = new_detected_Mild)
}

summarise_age_dependent <- function(df){
  df %>%
    group_by(date, age_group) %>%
    summarise(
      across(
        !replicate,
        mean,
      )
    )
}
