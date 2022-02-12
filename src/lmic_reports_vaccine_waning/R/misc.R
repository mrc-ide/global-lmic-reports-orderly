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
  summarised <- lapply(names(Rt_futures), function(trend){
    #for each replicate
    do.call(
      rbind,
      lapply(seq_along(Rt_futures[[trend]]), function(rep){
        #just get the mean
        mean(Rt_futures[[trend]][[rep]])
      })
    )
  })
  #summarise and put into dataframe
  out <- data.frame(temp = 1
  )
  for(i in seq_along(Rt_futures)){
    out <- out %>%
      mutate(
        !! paste0(names(Rt_futures)[i], "_med") := median(summarised[[i]]),
        !! paste0(names(Rt_futures)[i], "_025") := quantile(summarised[[i]], 0.025),
        !! paste0(names(Rt_futures)[i], "_975") := quantile(summarised[[i]], 0.975)
      )
  }
  return(out %>%
           select(!temp))
}

variant_changes_over_time <- function(variant_characteristics,
                                      var){

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
  return(
    list(
      var = variable,
      dates = dates
    )
  )
}

variant_immune_escape <- function(variant_characteristics, dur_R){

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

  #check if we need to drop any due to overlap
  return(
    list(
      variable = 2/c(dur_R_d, dur_R_o),
      dates = c(date_dur_R_change_d, date_dur_R_change_o)
    )
  )
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
  #calculate variant adjusted effiacis over time
  vaccine_efficacy_infection <- ((variant_characteristics$Wild$ve_infection[1] * (1 - delta_shift) +
    variant_characteristics$Delta$ve_infection[1] * delta_shift) *
    (1 - omicron_shift) +
    variant_characteristics$Omicron$ve_infection[1] * omicron_shift) *
    (1 - vacc_inputs$dose_ratio) +
    ((variant_characteristics$Wild$ve_infection[2] * (1 - delta_shift) +
        variant_characteristics$Delta$ve_infection[2] * delta_shift) *
       (1 - omicron_shift) +
       variant_characteristics$Omicron$ve_infection[2] * omicron_shift) *
    vacc_inputs$dose_ratio

  #scale for break through protection
  variant_characteristics$Wild$ve_disease <- (variant_characteristics$Wild$ve_disease - variant_characteristics$Wild$ve_infection)/
    (1- variant_characteristics$Wild$ve_infection)
  variant_characteristics$Delta$ve_disease <- (variant_characteristics$Delta$ve_disease - variant_characteristics$Delta$ve_infection)/
    (1- variant_characteristics$Delta$ve_infection)
  variant_characteristics$Omicron$ve_disease <- (variant_characteristics$Omicron$ve_disease - variant_characteristics$Omicron$ve_infection)/
    (1- variant_characteristics$Omicron$ve_infection)

  vaccine_efficacy_disease <- ((variant_characteristics$Wild$ve_disease[1] * (1 - delta_shift) +
                                    variant_characteristics$Delta$ve_disease[1] * delta_shift) *
                                   (1 - omicron_shift) +
                                   variant_characteristics$Omicron$ve_disease[1] * omicron_shift) *
    (1 - vacc_inputs$dose_ratio) +
    ((variant_characteristics$Wild$ve_disease[2] * (1 - delta_shift) +
        variant_characteristics$Delta$ve_disease[2] * delta_shift) *
       (1 - omicron_shift) +
       variant_characteristics$Omicron$ve_disease[2] * omicron_shift) *
    vacc_inputs$dose_ratio


  return(list(
    date_vaccine_change = vacc_inputs$date_vaccine_change,
    max_vaccine = vacc_inputs$max_vaccine,
    vaccine_efficacy_infection = lapply(vaccine_efficacy_infection,
                                        rep,
                                        17),
    vaccine_efficacy_disease = lapply(vaccine_efficacy_disease,
                                      rep,
                                      17),
    rel_infectiousness_vaccinated = rep(
      variant_characteristics$Wild$ve_transmission, 17
      ),
    vaccine_coverage_mat = vacc_inputs$vaccine_coverage_mat,
    dur_V = vacc_inputs$dur_V,
    dur_vaccine_delay = vacc_inputs$dur_vaccine_delay,
    strategy = vacc_inputs$strategy
  ))
}
