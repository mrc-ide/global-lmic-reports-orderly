r_list_format <- function(out, date_0) {

  df <- nimue_format(out,
                      var_select = c("infections","deaths","hospital_demand",
                                     "ICU_demand", "D", "hospital_incidence","ICU_incidence"),
                      date_0 = date_0)

  pr <- nimue_format(out, var_select = c("S","R","D"), date_0 = date_0) %>%
    na.omit %>%
    pivot_wider(names_from = compartment, values_from = y) %>%
    mutate(y = sum(out$parameters$population)-D-R-S,
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


variant_changes_over_time <- function(variant_characteristics,
                                      var, start_date) {

  #delta
  dates <- seq(
    variant_characteristics$Delta$start_date,
    variant_characteristics$Delta$start_date +
      variant_characteristics$Delta$shift_duration,
    by = 1
  )
  variable <- c(seq(variant_characteristics$Wild[[var]],
                    variant_characteristics$Delta[[var]],
                    length.out = length(dates)+1))
  #omicron
  dates_omicron <- seq(
    variant_characteristics$Omicron$start_date,
    variant_characteristics$Omicron$start_date +
      variant_characteristics$Omicron$shift_duration,
    by = 1
  )
  variable <- c(
    variable,
    seq(variant_characteristics$Delta[[var]],
        variant_characteristics$Omicron[[var]],
        length.out = length(dates_omicron)+1)[-1]
  )
  dates <- c(dates, dates_omicron)
  tt <- as.numeric(dates - start_date)
  if(tt[1] > 0){
    tt <- c(0, tt)
  } else {
    tt <- c(tt[1] - 1, tt)
  }
  return(
    list(
      var = variable,
      tt = tt
    )
  )
}

variant_immune_escape <- function(variant_characteristics, dur_R, start_date) {

  #delta
  dur_R_d <- c(dur_R, 1 / (
    (variant_characteristics$Delta$shift_duration / dur_R - log(1 - variant_characteristics$Delta$immune_escape)) /
      variant_characteristics$Delta$shift_duration
  ), dur_R)
  date_dur_R_change_d <- c(variant_characteristics$Delta$start_date,
                           variant_characteristics$Delta$start_date +
                             variant_characteristics$Delta$shift_duration)
  #omicron
  dur_R_o <- c( 1 / (
    (variant_characteristics$Omicron$shift_duration / dur_R - log(1 - variant_characteristics$Omicron$immune_escape)) /
      variant_characteristics$Omicron$shift_duration
  ), dur_R)
  date_dur_R_change_o <- c(variant_characteristics$Omicron$start_date,
                           variant_characteristics$Omicron$start_date +
                             variant_characteristics$Omicron$shift_duration)

  tt <- as.numeric(c(date_dur_R_change_d, date_dur_R_change_o) - start_date)
  if(tt[1] > 0){
    tt <- c(0, tt)
  } else {
    tt <- c(tt[1] - 1, tt)
  }
  return(
    list(
      variable = c(dur_R_d, dur_R_o),
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
