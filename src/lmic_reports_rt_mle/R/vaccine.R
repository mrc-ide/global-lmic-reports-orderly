get_vaccine_inputs <- function(iso3c) {
  vacc_inputs <- readRDS("vacc_inputs.Rds")[[iso3c]]
  vacc_inputs$date_vaccine_change <- vacc_inputs$date_vaccine_change
  vacc_inputs$first_doses <- as.numeric(vacc_inputs$first_doses)
  vacc_inputs$second_doses <- as.numeric(vacc_inputs$second_doses)
  vacc_inputs$booster_doses <- as.numeric(vacc_inputs$booster_doses)
  vacc_inputs$dur_V <- as.numeric(vacc_inputs$dur_V)
  return(vacc_inputs)
}


vaccine_eff_over_time <- function(vacc_inputs, variant_characteristics, start_date) {
  #time varying variant factors
  delta_shift <- case_when(
    vacc_inputs$date_vaccine_change < variant_characteristics$Delta$start_date ~
      0,
    vacc_inputs$date_vaccine_change > variant_characteristics$Delta$start_date +
      variant_characteristics$Delta$shift_duration~
      1,
    TRUE ~ as.numeric((vacc_inputs$date_vaccine_change -
                         variant_characteristics$Delta$start_date)/
                        variant_characteristics$Delta$shift_duration)
  )
  omicron_shift <- case_when(
    vacc_inputs$date_vaccine_change < variant_characteristics$Omicron$start_date ~
      0,
    vacc_inputs$date_vaccine_change > variant_characteristics$Omicron$start_date +
      variant_characteristics$Omicron$shift_duration~
      1,
    TRUE ~ as.numeric((vacc_inputs$date_vaccine_change -
                         variant_characteristics$Omicron$start_date)/
                        variant_characteristics$Omicron$shift_duration)
  )
  #reduce to values where things change
  indexes <- diff(c(-100, delta_shift + omicron_shift)) != 0
  date_vaccine_efficacy_change <- vacc_inputs$date_vaccine_change[indexes][-1]
  if(length(date_vaccine_efficacy_change) == 0){
    #if no changes then we don't modify it
    date_vaccine_efficacy_change <- NULL
  }
  delta_shift <- delta_shift[indexes]
  omicron_shift <- omicron_shift[indexes]

  #scale for break through protection
  variant_characteristics$Wild$ve_disease <- (variant_characteristics$Wild$ve_disease - variant_characteristics$Wild$ve_infection)/
    (1- variant_characteristics$Wild$ve_infection)
  variant_characteristics$Delta$ve_disease <- (variant_characteristics$Delta$ve_disease - variant_characteristics$Delta$ve_infection)/
    (1- variant_characteristics$Delta$ve_infection)
  variant_characteristics$Omicron$ve_disease <- (variant_characteristics$Omicron$ve_disease - variant_characteristics$Omicron$ve_infection)/
    (1- variant_characteristics$Omicron$ve_infection)

  inputed_values <- purrr::map(c(ve_disease = "ve_disease", ve_infection = "ve_infection", dur_V = "dur_V"), function(z){
    purrr::map(seq_along(delta_shift), function(x){
      #produce matrix for age group/vaccine status effiacies
      purrr::map_dbl(seq_along(variant_characteristics$Wild[[z]]), function(y){
        (variant_characteristics$Wild[[z]][y] * (1 - delta_shift[x]) +
           variant_characteristics$Delta[[z]][y] * delta_shift[x]) *
          (1 - omicron_shift[x]) +
          variant_characteristics$Omicron[[z]][y] * omicron_shift[x]
      })
    })
  })

  tt_vaccine_efficacy_change <- as.numeric(date_vaccine_efficacy_change - start_date)
  if(tt_vaccine_efficacy_change[1] > 0){
    tt_vaccine_efficacy_change <- c(0, tt_vaccine_efficacy_change)
  } else {
    tt_vaccine_efficacy_change <- c(tt_vaccine_efficacy_change[1] - 1, tt_vaccine_efficacy_change)
  }
  return(list(
    tt_vaccine_efficacy_change = tt_vaccine_efficacy_change,
    vaccine_efficacy_infection = inputed_values$ve_infection,
    vaccine_efficacy_disease = inputed_values$ve_disease,
    rel_infectiousness_vaccinated = variant_characteristics$Wild$ve_transmission,
    vaccine_coverage_mat = vacc_inputs$vaccine_coverage_mat,
    dur_V = inputed_values$dur_V,
    strategy = vacc_inputs$strategy
  ))
}
