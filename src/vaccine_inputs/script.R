date_0 <- as.Date(date)

### Get the vaccine data from online
#sort out dates from before date_0 and group data by country

vacc_types <- import_vaccine_agreements()

vdm <- import_vaccine_doses_by_manufacturer()

who_vacc <- import_who_vaccination()

who_vacc_meta <- import_who_vaccination_meta()

owid <- import_owid() %>%
  filter(date <= date_0)

## Modifications for single dose Vaccines
#we will treat as only a single dose of the other vaccines, more suitable
#with its lower efficacy
single_dose_df <- get_single_dose_data(
  owid, vacc_types, vdm, who_vacc, who_vacc_meta, single_dose_vaccines =
    c("CanSino", "Johnson&Johnson", "Janssen - Ad26.COV 2-S",
      "Janssen - Ad26.COV 2.5", "CanSino - Convidecia",
      "Gamaleya - Sputnik-Light")
  )

#make adjustments needed to OWID
owid <- owid %>%
  left_join(single_dose_df) %>%
  group_by(iso3c) %>%
  mutate(
    people_fully_vaccinated = if_else(
      (!is.na(people_fully_vaccinated)) & (!is.na(single_doses)),
      people_fully_vaccinated*(1-single_doses/max(people_fully_vaccinated, na.rm = TRUE)),
      people_fully_vaccinated
    ),
    people_fully_vaccinated_per_hundred = if_else(
      !is.na(people_fully_vaccinated) & !is.na(single_doses),
      people_fully_vaccinated_per_hundred*(1-single_doses/max(people_fully_vaccinated, na.rm = TRUE)),
      people_fully_vaccinated_per_hundred
    )) %>% #now to do countries that only had single dose vaccines
  mutate(
    people_fully_vaccinated = if_else(
      identical(single_dose_percentage, 1) & (!is.na(people_fully_vaccinated)),
      0,
      people_fully_vaccinated
    ),
    people_fully_vaccinated_per_hundred = if_else(
      identical(single_dose_percentage, 1) & (!is.na(people_fully_vaccinated_per_hundred)),
      0,
      people_fully_vaccinated_per_hundred
    ),
    total_vaccinations = if_else(
      identical(single_dose_percentage, 1) & (!is.na(total_vaccinations)),
      people_vaccinated,
      total_vaccinations
    ),
    total_vaccinations_per_hundred = if_else(
      identical(single_dose_percentage, 1) & (!is.na(total_vaccinations)),
      people_vaccinated_per_hundred,
      total_vaccinations_per_hundred
    )
  ) %>%
  select(!c(single_dose_percentage, single_doses)) %>%
  ungroup()

#adjustments to who data
who_vacc <- who_vacc %>%
  left_join(single_dose_df) %>%
  mutate(PERSONS_FULLY_VACCINATED = if_else(
    (!is.na(PERSONS_FULLY_VACCINATED)) & (!is.na(single_doses)),
    PERSONS_FULLY_VACCINATED*(1-single_doses/TOTAL_VACCINATIONS),
    as.numeric(PERSONS_FULLY_VACCINATED)
  ),
  PERSONS_FULLY_VACCINATED_PER100 = if_else(
    (!is.na(PERSONS_FULLY_VACCINATED_PER100)) & (!is.na(single_doses)),
    PERSONS_FULLY_VACCINATED_PER100*(1-single_doses/TOTAL_VACCINATIONS),
    as.numeric(PERSONS_FULLY_VACCINATED_PER100)
  )) %>%
  mutate(
    PERSONS_VACCINATED_1PLUS_DOSE = if_else(
      identical(single_dose_percentage, 1) & (!is.na(PERSONS_VACCINATED_1PLUS_DOSE)),
      TOTAL_VACCINATIONS,
      as.numeric(PERSONS_VACCINATED_1PLUS_DOSE)
    ),
    PERSONS_VACCINATED_1PLUS_DOSE_PER100 = if_else(
      identical(single_dose_percentage, 1) & (!is.na(PERSONS_VACCINATED_1PLUS_DOSE_PER100)),
      TOTAL_VACCINATIONS_PER100,
      PERSONS_VACCINATED_1PLUS_DOSE_PER100
    ),
    PERSONS_FULLY_VACCINATED = if_else(
      identical(single_dose_percentage, 1) & (!is.na(PERSONS_FULLY_VACCINATED)),
      0,
      PERSONS_FULLY_VACCINATED
    ),
    PERSONS_FULLY_VACCINATED_PER100 = if_else(
      identical(single_dose_percentage, 1) & (!is.na(PERSONS_FULLY_VACCINATED_PER100)),
      0,
      PERSONS_FULLY_VACCINATED_PER100
    )
  )

### Set up vaccine efficacies

dur_vaccine_delay <- 14
rel_infectiousness_vaccinated <- 0.5
dur_V <- Inf
ve_i_low <- 0.6
ve_i_high <- 0.8
ve_d_low <- 0.8
ve_d_high <- 0.98
days_between_doses <- 14

### Derive effective vaccine efficacies over time

iso3cs <- unique(squire::population$iso3c)

## Use the old function to get dose_ratio and max_vaccine
dose_df <- lapply(iso3cs, function(country){
  country_efficacy <- get_ret_res(owid %>% filter(iso3c == country) %>% select(!iso3c),
                                  who_vacc %>% filter(iso3c == country) %>% select(!iso3c),
                                  who_vacc_meta %>% filter(iso3c == country) %>% select(!iso3c),
                                  date_vaccine_change %>% filter(iso3cs == country) %>% select(!iso3c),
                                  ve_i_low, ve_i_high, ve_d_low, ve_d_high)
  #only interested in dates ratio and max_vaccine, we can recalculate efficacy
  #later on
  data.frame(
    iso3c = country,
    date_vaccine_change = country_efficacy$date_vaccine_change,
    max_vaccine = country_efficacy$max_vaccine,
    dose_ratio = country_efficacy$dose_ratio,
    vaccine_efficacy_infection = unlist(
      lapply(country_efficacy$vaccine_efficacy_infection, function(row){
        if(unique(row) > 1){
          stop("VE varys over age, please adjust code")
        }
        unique(row)
      })
    ),
    vaccine_efficacy_disease = unlist(
      lapply(country_efficacy$vaccine_efficacy_disease, function(row){
        if(unique(row) > 1){
          stop("VE varys over age, please adjust code")
        }
        unique(row)
      })
    )
  )
})

dose_df <- do.call(
  rbind,
  dose_df
  ) %>%
  group_by(iso3c)

#confirm that vaccine change every day
if(length(
  unique(unlist(lapply(iso3cs, function(x){unique(diff(filter(dose_df, iso3c==x)$date_vaccine_change))})))
)>1){
  stop("Code needs adjusting as vaccine changes less than everyday")
}

#make adjustment so that on the first day of vaccinations dose_ratio cannot be greater > 0
dose_df <- dose_df %>%
  mutate( #figure out when the first day vaccinations are
    doses_so_far = cumsum(max_vaccine),
    first_day = if_else(
      doses_so_far > 0 & (
        doses_so_far - lag(doses_so_far) == doses_so_far|
          lag(doses_so_far, default = -888) == -888
        ),
      TRUE,
      FALSE
    ),
    tag = if_else( #is this a problem?
      any(first_day & dose_ratio > 0),
      TRUE,
      FALSE
    )
  ) %>% #we don't actually need days before first day
  filter(doses_so_far > 0) #for those tagged we need to add an extra row
dose_df <- dose_df %>%
  #first and put half the first doses there and set dose_ratio to 0 for that day
  rbind(
    dose_df %>%
      filter(tag == 1, first_day) %>% #get the problematic first days
      mutate(
        date_vaccine_change = date_vaccine_change - 1, #set the date to the day before
        dose_ratio = 0, #dose ratio to 0
        max_vaccine = round(max_vaccine/3), #get 1 third of the vaccinations
        first_day = FALSE, #this is for utlity later
        vaccine_efficacy_infection = ve_i_low,
        vaccine_efficacy_disease = ve_d_low
        )
  ) %>%
  arrange(iso3c, date_vaccine_change) %>%
  mutate(#subtract the moved first doses from the original first date
    max_vaccine = if_else(
      tag == 1 & first_day,
      max_vaccine - lag(max_vaccine),
      max_vaccine
    )
  ) %>% #remove extra variables
  select(!c(tag, first_day, doses_so_far))

#extend data
#we'll assume that max_vaccine stays at its final weekly average
#dose ratio will continue to change at its current rate, until hit 0 or 1
dose_df <- dose_df %>% #add indicator for plotting later
  #extend values to date_0
  mutate(imputed = FALSE)%>%
  #add dates up to current date
  complete(date_vaccine_change = seq(min(date_vaccine_change), date_0, by = 1)) %>%
  left_join( #add averages for max_vaccine
    dose_df %>%
      group_by(iso3c) %>%
      filter(date_vaccine_change >= max(date_vaccine_change) - 7) %>%
      summarise(
        max_vaccine_week_ave = mean(max_vaccine)
      )
  ) %>%
  mutate(
    dose_ratio = if_else(
      is.na(dose_ratio),
      predict(
        lm(y~x, data.frame("x" = seq_along(dose_ratio), "y" = dose_ratio)),
        newdata = data.frame("x" = seq_along(dose_ratio))
      ) %>%
        vapply(min, numeric(1), 1) %>%
        vapply(max, numeric(1), 0),
      dose_ratio
    ),
    max_vaccine = if_else(
      is.na(max_vaccine),
      max_vaccine_week_ave,
      max_vaccine
    ),
    imputed = if_else(
      is.na(imputed),
      TRUE,
      FALSE
    )
  ) %>%
    select(!max_vaccine_week_ave)

## Apply Delta Adjustment if needed
#just set vaccine efficacies to constant values
dose_df <- dose_df %>% mutate(
  ve_i_high = ve_i_high,
  ve_d_high = ve_d_high,
  ve_d_low = ve_d_low,
  ve_i_low = ve_i_low
)
if(adjust_delta){
  dose_df <- add_delta_characteristics(dose_df)
}

if(!waning){
  #use the standard method of waning
  dur_V <- 5000
  #calculate efficacies for extended data
  dose_df <- dose_df %>%
    mutate(
      vaccine_efficacy_disease =
        ve_d_low * (1-dose_ratio) +
        ve_d_high * dose_ratio,
      vaccine_efficacy_infection =
        ve_i_low * (1-dose_ratio) +
        ve_i_high * dose_ratio
    ) %>%
    select(iso3c, date_vaccine_change, max_vaccine, dose_ratio,
           vaccine_efficacy_infection, vaccine_efficacy_disease,
           imputed, any_of(c("shift_start", "shift_end")))
} else {
  gen <- function(lambda, first_doses, per_first_who_get_second, n_dates){
    second_doses <- c(
      rep(0, lambda),
      head(x=per_first_who_get_second, n = n_dates - lambda)
    )
    #calc dose ratio
    if(length(second_doses) != length(first_doses)){
      stop()
    }
    cumsum(second_doses)/cumsum(first_doses)
  }
  #make an estimate of the days between doses for each country
  days_between_doses <- lapply(unique(dose_df$iso3c), function(country){
    this_country <- filter(dose_df, iso3c == country)
    if(nrow(this_country) > 1){
      first_doses <- this_country$max_vaccine
      second_dose_percentage <- tail(this_country$dose_ratio, 1)
      per_first_who_get_second <- first_doses*second_dose_percentage
      n_dates <- length(first_doses)
      err_func <- function(lambda){
        sqrt(sum((this_country$dose_ratio -
                    gen(lambda, first_doses, per_first_who_get_second, n_dates))^2))
      }
      lambdas <- seq(1, length(first_doses), by = 1)
      max(lambdas[which.min(unlist(lapply(lambdas, err_func)))],
          28)
    } else {
      NA
    }
  })
  names(days_between_doses) <- unique(dose_df$iso3c)
  #add to data
  dose_df <- mutate(
    dose_df,
    days_between_doses = days_between_doses[[iso3c[1]]]
  )
  #calculate efficacy on each day for each country with waning and dose adjustment
  dose_df <- dose_df %>%
    arrange(iso3c, date_vaccine_change) %>%
    calculate_waning_eff(dates = date_vaccine_change,
                         first_doses = max_vaccine,
                         efficacy_infection_first = ve_i_low,
                         efficacy_infection_second = ve_i_high,
                         efficacy_disease_first = ve_d_low,
                         efficacy_disease_second = ve_d_high,
                         dur_vaccine_delay = dur_vaccine_delay,
                         days_between_doses = days_between_doses,
                         countries = iso3c,
                         diagnostic = TRUE)
}

### Calibration plot
dir.create("calibration", showWarnings = FALSE)
pdf("calibration/plot.pdf", paper = "a4r")
if(waning){
  #plot efficacy curves for our first and second doses
  print(plot_waning(adjust_delta = adjust_delta, days_between_doses = 28))
}

#plot for each country
for(country in iso3cs){
  this_country <- dose_df %>%
    filter(iso3c == country) %>%
    rename(Date = date_vaccine_change)
  if(nrow(this_country) > 1){
    #if there is data

    #plot of new doses each day
    dose_plot <- plot_doses(this_country)
    #plot dose ratio
    dose_ratio_plot <- plot_dose_ratio(this_country)
    #plot efficacies not accounting for dose_ratio

    #final efficacy plots
    final_eff_plot <- plot_efficacy(this_country, adjust_delta)

    if(waning){
      title_plot <- as_ggplot(text_grob(
        paste0(country, ", Assumed days between doses: ", unique(this_country$days_between_doses)) ,size = 20))
      split_eff <- plot_efficacy_split(this_country, adjust_delta)
      suppressWarnings(
        print(ggarrange(
          title_plot,
          ggarrange(
            ggarrange(
              dose_plot,
              dose_ratio_plot,
              ncol = 1
            ),
            split_eff
          ),
          final_eff_plot,
          ncol = 1,
          heights = c(0.1,1,1.2)
        ))
      )
    } else {
      title_plot <- as_ggplot(text_grob(
        country ,size = 20))
      suppressWarnings(
        print(ggarrange(
          title_plot,
            ggarrange(
              dose_plot,
              dose_ratio_plot
            ),
          final_eff_plot,
          ncol = 1,
          heights = c(0.1,1,1.2)
        ))
      )
    }
  }
}
dev.off()

### Export

#reduce to dates that relevent <date_0 and variables that we need
dose_df <- dose_df %>%
  filter(date_vaccine_change<=date_0) %>%
  select(iso3c, date_vaccine_change, max_vaccine, vaccine_efficacy_disease,
         vaccine_efficacy_infection, imputed) %>% #adjust disease efficacy for infection efficacy
  mutate(vaccine_efficacy_disease = (vaccine_efficacy_disease - vaccine_efficacy_infection)/
           (1- vaccine_efficacy_infection))

#into list format
dose_list <-
  lapply(iso3cs,
         function(country){
           country_df <- dose_df %>%
             filter(iso3c == country)
           list(
             dur_vaccine_delay = dur_vaccine_delay,
             dur_V = dur_V,
             date_vaccine_change = country_df$date_vaccine_change,
             vaccine_efficacy_disease = lapply(country_df$vaccine_efficacy_disease,
                                               rep,
                                               17),
             vaccine_efficacy_infection = lapply(country_df$vaccine_efficacy_infection,
                                                 rep,
                                                 17),
             rel_infectiousness_vaccinated = rep(rel_infectiousness_vaccinated, 17),
             max_vaccine = country_df$max_vaccine,
             imputed = country_df$imputed
           )
         }
)
names(dose_list) <- iso3cs

saveRDS(dose_list,
        "vacc_inputs.Rds"
        )

