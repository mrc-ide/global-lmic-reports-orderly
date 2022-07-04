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
      start_value <- 1
    } else {
      start_value <- var[[variant_timings$variant[i - 1]]]
    }
    end_value <- var[[variant_timings$variant[i]]]
    start_value * (1 - change_tt[[i]]$change) + end_value * change_tt[[i]]$change
  }) %>% unlist()

  #add the zero entry if needed
  if(all(tt > 0)){
    tt <- c(0, tt)
    var <- c(1, var)
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
