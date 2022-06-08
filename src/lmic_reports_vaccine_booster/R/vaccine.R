get_vaccine_inputs <- function(iso3c){
  vacc_inputs <- readRDS("vacc_inputs.Rds")[[iso3c]]
  vacc_inputs$date_vaccine_change <- vacc_inputs$date_vaccine_change
  vacc_inputs$first_doses <- as.numeric(vacc_inputs$first_doses)
  vacc_inputs$second_doses <- as.numeric(vacc_inputs$second_doses)
  vacc_inputs$booster_doses <- as.numeric(vacc_inputs$booster_doses)
  vacc_inputs$dur_V <- as.numeric(vacc_inputs$dur_V)
  return(vacc_inputs)
}


vaccine_eff_over_time <- function(vacc_inputs, variant_characteristics){
  #check if we need to update dates
  if(variant_characteristics$Omicron$start_date <
     variant_characteristics$Delta$start_date +
     variant_characteristics$Delta$shift_duration){
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
  #calculate variant adjusted effiacis over time
  vaccine_efficacy_infection <- lapply(seq_along(delta_shift), function(x){
    #produce matrix for age group/vaccine status effiacies
    purrr::map_dbl(
      seq_along(variant_characteristics$Wild$ve_infection),
      ~(variant_characteristics$Wild$ve_infection[.x] * (1 - delta_shift[x]) +
        variant_characteristics$Delta$ve_infection[.x] * delta_shift[x]) *
      (1 - omicron_shift[x]) +
      variant_characteristics$Omicron$ve_infection[.x] * omicron_shift[x]
    )
  })

  #scale for break through protection
  variant_characteristics$Wild$ve_disease <- (variant_characteristics$Wild$ve_disease - variant_characteristics$Wild$ve_infection)/
    (1- variant_characteristics$Wild$ve_infection)
  variant_characteristics$Delta$ve_disease <- (variant_characteristics$Delta$ve_disease - variant_characteristics$Delta$ve_infection)/
    (1- variant_characteristics$Delta$ve_infection)
  variant_characteristics$Omicron$ve_disease <- (variant_characteristics$Omicron$ve_disease - variant_characteristics$Omicron$ve_infection)/
    (1- variant_characteristics$Omicron$ve_infection)

  vaccine_efficacy_disease <- lapply(seq_along(delta_shift), function(x){
    #produce matrix for age group/vaccine status effiacies
    purrr::map_dbl(
      seq_along(variant_characteristics$Wild$ve_disease),
      ~(variant_characteristics$Wild$ve_disease[.x] * (1 - delta_shift[x]) +
          variant_characteristics$Delta$ve_disease[.x] * delta_shift[x]) *
        (1 - omicron_shift[x]) +
        variant_characteristics$Omicron$ve_disease[.x] * omicron_shift[x]
    )
  })

  #durations of vaccine protection
  dur_V <- lapply(seq_along(delta_shift), function(x){
    #produce matrix for age group/vaccine status effiacies
    purrr::map_dbl(
      seq_along(variant_characteristics$Wild$dur_V),
      ~(variant_characteristics$Wild$dur_V[.x] * (1 - delta_shift[x]) +
          variant_characteristics$Delta$dur_V[.x] * delta_shift[x]) *
        (1 - omicron_shift[x]) +
        variant_characteristics$Omicron$dur_V[.x] * omicron_shift[x]
    )
  })

  return(list(
    date_vaccine_change = vacc_inputs$date_vaccine_change,
    first_doses = vacc_inputs$first_doses,
    second_doses = vacc_inputs$second_doses,
    booster_doses = vacc_inputs$booster_doses,
    date_vaccine_efficacy_change = date_vaccine_efficacy_change,
    vaccine_efficacy_infection = vaccine_efficacy_infection,
    vaccine_efficacy_disease = vaccine_efficacy_disease,
    rel_infectiousness_vaccinated = variant_characteristics$Wild$ve_transmission,
    vaccine_coverage_mat = vacc_inputs$vaccine_coverage_mat,
    dur_V = dur_V,
    strategy = vacc_inputs$strategy
  ))
}

