get_vaccine_inputs <- function(iso3c){
  vacc_inputs <- readRDS("vacc_inputs.Rds")[[iso3c]]
  vacc_inputs$date_vaccine_change <- vacc_inputs$date_vaccine_change
  vacc_inputs$max_vaccine <- as.numeric(vacc_inputs$max_vaccine)
  vacc_inputs$vaccine_efficacy_disease <- vacc_inputs$vaccine_efficacy_disease
  vacc_inputs$vaccine_efficacy_infection <- vacc_inputs$vaccine_efficacy_infection
  vacc_inputs$dur_V <- as.numeric(vacc_inputs$dur_V)
  vacc_inputs$dur_vaccine_delay <- as.numeric(vacc_inputs$dur_vaccine_delay)
  return(vacc_inputs)
}
