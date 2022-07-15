get_vaccine_inputs <- function(iso3c) {
  vacc_inputs <- readRDS("vacc_inputs.Rds")[[iso3c]]
  vacc_inputs$date_vaccine_change <- vacc_inputs$date_vaccine_change
  vacc_inputs$primary_doses <- as.numeric(vacc_inputs$primary_doses)
  if(length(vacc_inputs$second_dose_delay) == 0){
    vacc_inputs$second_dose_delay <- 60
  } else {
    vacc_inputs$second_dose_delay <- as.numeric(vacc_inputs$second_dose_delay)
  }
  vacc_inputs$booster_doses <- as.numeric(vacc_inputs$booster_doses)
  return(vacc_inputs)
}


vaccine_eff_over_time <- function(variant_timings, variant_ve, x, start_date) {
  variant_timings <- arrange(variant_timings, start_date)

  change_tt <- map(transpose(variant_timings), function(l){
    period <- as.numeric(as_date(c(l$start_date, l$end_date)) - start_date)
    tt <- seq(period[1], period[2], by = 1)
    list(
      change = seq(0, 1, length.out = length(tt) + 1)[-1],
      tt = tt
    )
  })
  tt <- map(change_tt, ~.x$tt) %>% unlist()
  vars <- c("dur_V", "vaccine_efficacy_infection", "vaccine_efficacy_disease")
  names(vars) <- vars
  vars <- map(vars, function(var){
    per_variant <- map(seq_along(change_tt), function(i){
      if(i == 1){
        start_values <- variant_ve$Wild[[var]][[x]]
      } else {
        start_values <- variant_ve[[variant_timings$variant[i-1]]][[var]][[x]]
      }
      end_values <- variant_ve[[variant_timings$variant[i]]][[var]][[x]]
      map(
        seq_along(start_values),
        ~start_values[.x] * (1 - change_tt[[i]]$change) + end_values[.x] * change_tt[[i]]$change
      ) %>% transpose() %>%
        map(unlist)
    })
    #now append these togther
    unlist(per_variant, recursive = FALSE)
  })

  #add the zero entry if needed
  if(all(tt > 0)){
    tt <- c(0, tt)
    new_vars <- map(names(vars), ~c(list(variant_ve$Wild[[.x]][[x]]), vars[[.x]]))
    names(new_vars) <- names(vars)
    vars <- new_vars
    rm(new_vars)
  }

  return(
    c(vars, list(tt = tt))
  )
}
