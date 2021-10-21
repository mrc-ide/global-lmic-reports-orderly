#a function when given a date when vaccinated (2nd dose) and the present date,
#returns the level of protection
get_eff_disease <- function(days_from_vacc, efficacy){
  #for now use Toms RTM values
  reduction_per_day <- 0.9984446
  no_waning_period <- 42
  if_else(days_from_vacc < no_waning_period,
          efficacy,
          efficacy*reduction_per_day^(days_from_vacc - no_waning_period)
          )
}
get_eff_disease_second <- function(days_from_vacc, days_between_doses,
                                     efficacy_first, efficacy_second){
  c(
    get_eff_disease(
      days_from_vacc[days_from_vacc<days_between_doses],
      efficacy_first),
    get_eff_disease(
      days_from_vacc[days_from_vacc>=days_between_doses],
      efficacy_second)
  )
}
get_eff_infection <- function(days_from_vacc, efficacy){
  #for now use Toms RTM values
  reduction_per_day <- 0.9970089
  no_waning_period <- 42
  if_else(days_from_vacc < no_waning_period,
          efficacy,
          efficacy*reduction_per_day^(days_from_vacc - no_waning_period)
  )
}
get_eff_infection_second <- function(days_from_vacc, days_between_doses,
                                     efficacy_first, efficacy_second){
  c(
    get_eff_infection(
      days_from_vacc[days_from_vacc<days_between_doses],
      efficacy_first),
    get_eff_infection(
      days_from_vacc[days_from_vacc>=days_between_doses],
      efficacy_second)
  )
}
#this function controls the dlay in tentry
vaccine_delay <- function(days_from_vacc, delay){
  pgamma(days_from_vacc, rate = 1/delay, shape = 2)
}
#function to calculate the waning efficacy for a country
calculate_waning_eff <- function(
  data, #data frame countaining info
  dates, #name of the variable with the dates
  first_doses, #name of variable of timeseries of first doses
  efficacy_infection_first, #names of variables of the efficacy for each
  efficacy_infection_second, #type of protection/dose
  efficacy_disease_first,
  efficacy_disease_second,
  dur_vaccine_delay, #the level of delay in first dose protection (numeric)
  days_between_doses, #the assumed delay between first dose and second dose
  countries = NULL,#name of the country variables, if null not used
  diagnostic = FALSE #return INTERNAL variables
){
  #inital prep
  if(dur_vaccine_delay == 0){
    delay_func <- function(days_from_first_dose){
      rep(1, length(days_from_first_dose))
    }
  }else{
    delay_func <- function(days_from_first_dose){
      vaccine_delay(days_from_first_dose, dur_vaccine_delay)
    }
  }
  data <-  group_by(data, {{ countries }})
  #calculate effective vaccine efficacy for each time and country
  data <- data %>% mutate(
    #calculate for disease
    INTERNAL = #get all combinations of dates
      expand.grid({{ dates }},
                  {{ dates }}) %>%
      cbind(first_dose =  {{ first_doses }}, #add doses and efficacies
            ve_d_f =  {{ efficacy_disease_first}},
            ve_d_s =  {{ efficacy_disease_second }},
            ve_i_f =  {{ efficacy_infection_first}},
            ve_i_s =  {{ efficacy_infection_second }}) %>%
      mutate( #calculate the days from
        days_from_vacc = as.numeric(Var2 - Var1),
        date = Var2
      )  %>%
      filter(days_from_vacc >= 0) %>% #remove days that are after our date
      arrange(date, days_from_vacc) %>% #arrange by date
      group_by(date) %>% #group by date
      mutate(#adjust doses for the delay
        first_dose_delay = delay_func(days_from_vacc)*first_dose,
        #get the latest efficacy data since this is the one that applies at this
        #date
        ve_d_f = head(ve_d_f, 1),
        ve_d_s = head(ve_d_s, 1),
        ve_i_f = head(ve_i_f, 1),
        ve_i_s = head(ve_i_s, 1)
      ) %>%
      summarise( #calculate the first/second dose efficacy for each type
        vaccine_efficacy_disease_first = sum(
          first_dose_delay*get_eff_disease(days_from_vacc, head(ve_d_f,1))
        )/sum(first_dose_delay),
        vaccine_efficacy_infection_first = sum(
          first_dose_delay*get_eff_infection(days_from_vacc, head(ve_i_f, 1))
        )/sum(first_dose_delay),
        vaccine_efficacy_disease_second = sum(
          first_dose_delay*get_eff_disease_second(days_from_vacc,
                                                  head({{ days_between_doses }}, 1),
                                                  head(ve_d_f, 1),
                                                  head(ve_d_s, 1))
        )/sum(first_dose_delay),
        vaccine_efficacy_infection_second = sum(
          first_dose_delay*get_eff_infection_second(days_from_vacc,
                                                    head({{ days_between_doses }}, 1),
                                                    head(ve_i_f, 1),
                                                    head(ve_i_s, 1))
        )/sum(first_dose_delay),
      ) %>%
      ungroup() %>%
      select(!date)
  ) %>%
    mutate(
      #calculate final values
      vaccine_efficacy_disease = INTERNAL$vaccine_efficacy_disease_first*(1-dose_ratio) +
        INTERNAL$vaccine_efficacy_disease_second*dose_ratio,
      vaccine_efficacy_infection = INTERNAL$vaccine_efficacy_infection_first*(1-dose_ratio) +
        INTERNAL$vaccine_efficacy_infection_second*dose_ratio,
      #nan to 0
      vaccine_efficacy_disease = if_else(
        is.na(vaccine_efficacy_disease) | is.nan(vaccine_efficacy_disease),
        {{ efficacy_disease_first}},
        vaccine_efficacy_disease
      ),
      vaccine_efficacy_infection = if_else(
        is.na(vaccine_efficacy_infection) | is.nan(vaccine_efficacy_infection),
        {{ efficacy_infection_first}},
        vaccine_efficacy_infection
      )
    )
  if(!diagnostic){
    data %>%
      select(!contains("INTERNAL"))
  }
  return(data)
}
