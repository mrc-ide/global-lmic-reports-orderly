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
get_eff_infection <- function(days_from_vacc, efficacy){
  #for now use Toms RTM values
  reduction_per_day <- 0.9970089
  no_waning_period <- 42
  if_else(days_from_vacc < no_waning_period,
          efficacy,
          efficacy*reduction_per_day^(days_from_vacc - no_waning_period)
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
  second_doses, #name of variable of timeseries when people get second doses
  efficacy_infection_first, #names of variables of the efficacy for each
  efficacy_infection_second, #type of protection/dose
  efficacy_disease_first,
  efficacy_disease_second,
  dur_vaccine_delay, #the level of delay in first dose protection (numeric)
  countries = NULL,#name of the country variables, if null not used
  diagnostic = FALSE #return INTERNAL variables
){
  #set up names for easier calling
  data <- data %>%
    mutate(
      INTERNAL_date = .data[[dates]],
      INTERNAL_first_dose = .data[[first_doses]],
      INTERNAL_second_dose = .data[[second_doses]],
      INTERNAL_ve_i_f = .data[[efficacy_infection_first]],
      INTERNAL_ve_i_s = .data[[efficacy_infection_second]],
      INTERNAL_ve_d_f = .data[[efficacy_disease_first]],
      INTERNAL_ve_d_s = .data[[efficacy_disease_second]],
      INTERNAL_country = .data[[countries]]
    )  %>%
    arrange(INTERNAL_country, INTERNAL_date) %>%
    mutate(
      #calculate the mean time of first dose for each time a second dose is received
      INTERNAL_mean_first_dose_date = unlist(
        lapply(
          INTERNAL_date,
          function(date){
            prev_dates <- INTERNAL_date[INTERNAL_date < date]
            prev_first_doses <- INTERNAL_first_dose[INTERNAL_date < date]
            weighted.mean(prev_dates, prev_first_doses)
          }
        )
      )
    )
  #do by country if relevant
  if(!is.null(iso3cs)){
    data <- data %>% group_by(INTERNAL_country)
  }
  #if relevant (dur_vaccine_delay > 0) set up functions
  if(dur_vaccine_delay == 0){
    delay_func <- function(days_from_vacc){
      rep(1, length(days_from_vacc))
    }
    data <- data %>%
      mutate(
        INTERNAL_cum_first_in_comp = cumsum(INTERNAL_first_dose),
        INTERNAL_cum_second_in_comp = cumsum(INTERNAL_second_dose)
      )
  } else {
    delay_func <- function(days_from_vacc){
      vaccine_delay(days_from_vacc, delay = dur_vaccine_delay)
    }
    #calculate entries into first dose compartment with delay
    data <- data %>%
      mutate(
        INTERNAL_cum_first_in_comp = expand.grid(INTERNAL_date,
                                                 INTERNAL_date) %>%
          cbind(dose_values = INTERNAL_first_dose) %>%
          mutate(diff = as.numeric(Var2 - Var1),
                 date = Var2) %>%
          filter(diff >= 0) %>%
          group_by(date) %>%
          summarise(
            value = sum(delay_func(diff) * dose_values)
          ) %>%
          pull(value)
      )
  }
  #some extra things for plotting utility
  if(diagnostic){
    data <- data %>%
      mutate(#people in compartment with second doses
        INTERNAL_cum_second_in_comp = expand.grid(INTERNAL_date,
                                                  INTERNAL_date) %>%
          cbind(dose_values = INTERNAL_second_dose,
                first_dose_date = INTERNAL_mean_first_dose_date) %>%
          mutate(date = Var2,
                 diff_2 = as.numeric(date - first_dose_date),
                 diff_1 = as.numeric(date - Var1),
                 diff_2 = if_else(diff_1 < 0 | is.na(diff_2),
                                0,
                                diff_2)) %>%
          group_by(date) %>%
          summarise(
            value = sum(delay_func(diff_2) * dose_values)
          ) %>%
          pull(value), #inferred dose_ratio
        INTERNAL_dose_ratio_comp = INTERNAL_cum_second_in_comp/INTERNAL_cum_first_in_comp
      )
  }
  #calculate effective vaccine efficacy for each time and country
  data <- data %>%
    mutate(
      #calculate for disease
      value = #get all combinations of dates
        expand.grid(INTERNAL_date,
                      INTERNAL_date) %>%
        cbind(first_dose =  INTERNAL_first_dose, #add doses and efficacies
              second_dose =  INTERNAL_second_dose,
              cum_first_in_comp =  INTERNAL_cum_first_in_comp,
              mean_first_dose_date =  INTERNAL_mean_first_dose_date,
              ve_d_f =  INTERNAL_ve_d_f,
              ve_d_s =  INTERNAL_ve_d_s,
              ve_i_f =  INTERNAL_ve_i_f,
              ve_i_s =  INTERNAL_ve_i_s) %>%
        mutate( #calculate the days from
          days_from_vacc = as.numeric(Var2 - Var1),
          date = Var2
        )  %>%
        filter(days_from_vacc >= 0) %>% #remove days that are after our date
        arrange(date, days_from_vacc) %>% #arrange by date
        group_by(date) %>% #group by date
        mutate(#calculate the things that will be same/used twice in the summarising
          cum_first_in_comp_max = max(cum_first_in_comp),
          days_from_mean_vacc = as.numeric(date - mean_first_dose_date)
        ) %>%
        summarise( #calculate the two parts of the efficacy
          f_ve_d = sum(
            get_eff_disease(days_from_vacc, ve_d_f) *
              delay_func(days_from_vacc) *
              first_dose,
            na.rm = TRUE
          )/head(cum_first_in_comp_max, 1), #this maximum gets the cumulative value at the current time
          s_ve_d_diff = sum(
            #calculate the difference from the average first does time
            (get_eff_disease(days_from_vacc, ve_d_s) -
               get_eff_disease(days_from_mean_vacc, ve_d_f))*
              delay_func(days_from_mean_vacc)*second_dose,
            na.rm = TRUE #this deals with the case where there are no second doses
          )/head(cum_first_in_comp_max, 1), #divide by people with first dose to get
          #dose ratio adjusted addition
          f_ve_i = sum(
            get_eff_disease(days_from_vacc, ve_i_f) *
              delay_func(days_from_vacc) *
              first_dose,
            na.rm = TRUE
          )/head(cum_first_in_comp_max, 1), #this maximum gets the cumulative value at the current time
          s_ve_i_diff = sum(
            #calculate the difference from the average first does time
            (get_eff_disease(days_from_vacc, ve_i_s) -
               get_eff_disease(days_from_mean_vacc, ve_i_f))*
              delay_func(days_from_mean_vacc)*second_dose,
            na.rm = TRUE #this deals with the case where there are no second doses
          )/head(cum_first_in_comp_max, 1), #divide by people with first dose to get
          #dose ratio adjusted addition
        ) %>%
        mutate(#get final efficacy
          vaccine_efficacy_disease = if_else(
            !is.nan(s_ve_d_diff) | !is.na(s_ve_d_diff),
            f_ve_d + s_ve_d_diff,
            f_ve_d
          ),
          vaccine_efficacy_infection = if_else(
            !is.nan(s_ve_i_diff) | !is.na(s_ve_i_diff),
            f_ve_i + s_ve_i_diff,
            f_ve_i
          )
        ) %>%
        #rowwise() %>%
        mutate(
          #merge to extract later
          out = list(cbind(vaccine_efficacy_disease, vaccine_efficacy_infection))
        ) %>% head(1) %>% pull(out),
      vaccine_efficacy_disease = value[[1]][,"vaccine_efficacy_disease"],
      vaccine_efficacy_infection = value[[1]][,"vaccine_efficacy_infection"],
      value = NULL
    )
  if(!diagnostic){
    data %>%
      select(!contains("INTERNAL"))
  }
  return(data)
}
