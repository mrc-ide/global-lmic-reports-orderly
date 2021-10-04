#a function when given a date when vaccinated (2nd dose) and the present date,
#returns the level of protection
get_eff_disease <- function(date, vaccination_date, efficacy){
  #for now use Toms RTM values
  diff_times <- as.numeric(date - vaccination_date)
  reduction_per_day <- 0.9984446
  no_waning_period <- 42
  if_else(diff_times < no_waning_period,
          efficacy,
          efficacy*reduction_per_day^(diff_times - no_waning_period)
          )
}
get_eff_infection <- function(date, vaccination_date, efficacy){
  #for now use Toms RTM values
  diff_times <- as.numeric(date - vaccination_date)
  reduction_per_day <- 0.9970089
  no_waning_period <- 42
  if_else(diff_times < no_waning_period,
          efficacy,
          efficacy*reduction_per_day^(diff_times - no_waning_period)
  )
}
#calculates vaccine efficacy with waning for a country
calculate_waning_eff_country <- function(date, dose_ratio,
                                         first_doses, second_doses,
                                         efficacy_infection_first,
                                         efficacy_infection_second,
                                         efficacy_disease_first,
                                         efficacy_disease_second,
                                         diagnostic = FALSE){
  #calculate efficacy first dose
  if(any(first_doses > 0)){
    #if there are vaccinations
    #when is the first dose?
    start_doses_first <- min(which(first_doses > 0))
    if(any(second_doses > 0)){
      #when is the second dose?
      start_doses_second <- min(which(second_doses > 0))
    } else {
      start_doses_second <- Inf
    }
    do.call(
      rbind, #rbind the list into a df
      lapply(start_doses_first:length(date), function(row){
        #for rows with post vaccinations
        #get the date at the moment
        cur_date <- date[row]
        #get the data in the past
        index_first <- start_doses_first:row
        #calculate efficacies
        ve_i_first <- get_eff_infection(
          cur_date,
          date[index_first],
          efficacy_infection_first[row]
        )
        ve_d_first <- get_eff_disease(
          cur_date,
          date[index_first],
          efficacy_disease_first[row]
        )
        out_df <- data.frame(date_vaccine_change = cur_date,
                             vaccine_efficacy_infection_first = weighted.mean(
                               ve_i_first,
                               first_doses[index_first]
                             ),
                             vaccine_efficacy_disease_first = weighted.mean(
                               ve_d_first,
                               first_doses[index_first]
                             )
        )
        if(start_doses_second <= row){
          #if second doses have occured
          index_second <- start_doses_second:row
          index_cross_over <- which(index_first %in% index_second)
          #calculating effect of second vaccine
          ve_i_second <- get_eff_disease(
            cur_date,
            date[index_second],
            efficacy_infection_second[row]
          )
          ve_d_second <- get_eff_disease(
            cur_date,
            date[index_second],
            efficacy_disease_second[row]
          )
          out_df <- out_df %>% mutate(
            #calculate averages
            vaccine_efficacy_infection_second_diff = weighted.mean(
              ve_i_second,
              second_doses[index_second]
            ) - vaccine_efficacy_infection_first,
            vaccine_efficacy_disease_second_diff = weighted.mean(
              ve_d_second,
              second_doses[index_second]
            ) - vaccine_efficacy_disease_first,
            #calculate overall vaccine efficacy with dose_ratio
            vaccine_efficacy_infection = vaccine_efficacy_infection_first +
            vaccine_efficacy_infection_second_diff * dose_ratio[row],
            vaccine_efficacy_disease = vaccine_efficacy_disease_first +
              vaccine_efficacy_disease_second_diff * dose_ratio[row]
          )

        } else {
          out_df <- out_df %>% mutate(
            vaccine_efficacy_infection_second_diff = NA,
            vaccine_efficacy_disease_second_diff = NA,
            vaccine_efficacy_infection = vaccine_efficacy_infection_first,
            vaccine_efficacy_disease = vaccine_efficacy_disease_first
          )
        }
        if(!diagnostic){
          out_df <- out_df %>%
            select(
              date_vaccine_change,
              vaccine_efficacy_infection,
              vaccine_efficacy_disease
            )
        }
        out_df
      }
      )
    )
    } else {
      if(!diagnostic){
        data.frame(
          date_vaccine_change = date[1],
          vaccine_efficacy_infection = NA,
          vaccine_efficacy_disease = NA
        )
      }
    data.frame(
      date_vaccine_change = date[1],
      vaccine_efficacy_infection_first = NA,
      vaccine_efficacy_infection_second_diff = NA,
      vaccine_efficacy_disease_first = NA,
      vaccine_efficacy_disease_second_diff = NA,
      vaccine_efficacy_infection = NA,
      vaccine_efficacy_disease = NA
    )
  }
}
calculate_waning_eff <- function(data, countries, dose_ratio,
                                 first_doses, second_doses,
                                 efficacy_infection_first,
                                 efficacy_infection_second,
                                 efficacy_disease_first,
                                 efficacy_disease_second){
  data %>%
    full_join(
      do.call(
        rbind,
        lapply(countries, function(country){
          data <- data %>% ungroup() %>%  filter(iso3c == country)
            calculate_waning_eff_country(date = data$date_vaccine_change,
              dose_ratio = data[[dose_ratio]],
              first_doses = data[[first_doses]],
              second_doses = data[[second_doses]],
              efficacy_infection_first = data[[efficacy_infection_first]],
              efficacy_infection_second = data[[efficacy_infection_second]],
              efficacy_disease_first = data[[efficacy_disease_first]],
              efficacy_disease_second = data[[efficacy_disease_second]],
              diagnostic = TRUE
            ) %>%
              mutate(iso3c = country)
        })
      )
    ) #%>% #any values not filled are by design 0
    # mutate(
    #   vaccine_efficacy_infection = if_else(
    #     is.na(vaccine_efficacy_infection),
    #     0,
    #     vaccine_efficacy_infection
    #   ),
    #   vaccine_efficacy_disease = if_else(
    #     is.na(vaccine_efficacy_disease),
    #     0,
    #     vaccine_efficacy_disease
    #   )
    # )
}
