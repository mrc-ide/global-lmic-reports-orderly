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
owid <- single_dose_adjust_owid(owid, single_dose_df)

#adjustments to who data
who_vacc <- single_dose_adjust_who_vacc(who_vacc, single_dose_df)

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

#get our doses over time
dose_df <- get_dose_ts(owid, who_vacc, who_vacc_meta, iso3cs, date_0)

## Apply Delta Adjustment if needed
#just set vaccine efficacies to constant values
dose_df <- dose_df %>% mutate(
  ve_i_high = ve_i_high,
  ve_d_high = ve_d_high,
  ve_d_low = ve_d_low,
  ve_i_low = ve_i_low
) %>%
  rename(date_vaccine_change = date)

if(adjust_delta){
  dose_df <- add_delta_characteristics(dose_df)
}

if(!waning){
  #use the standard method of waning
  dur_V <- 446
  #calculate efficacies
  dose_df <- dose_df %>%
    rename(
      max_vaccine = first_dose_per_day
    ) %>%
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

  fitting_df <- dose_df %>%
    mutate(
      first_dose_per_day = if_else(
        imputed,
        as.numeric(NA),
        first_dose_per_day
      ),
      second_dose_per_day = if_else(
        imputed,
        as.numeric(NA),
        second_dose_per_day
      )
    ) %>%
    select(iso3c, date_vaccine_change, first_dose_per_day, second_dose_per_day)

  #make an estimate of the average days between doses for each country
  delays_df <-
    data.frame(
      iso3c = iso3cs
    ) %>%
    mutate(
      delay = unlist(lapply(iso3cs, function(country){
        country_df <- filter(fitting_df, iso3c == country)
        if(nrow(country_df) > 1){
          values <- ccf(country_df$first_dose_per_day, country_df$second_dose_per_day, na.action = na.pass,
                        lag.max = nrow(country_df), plot = FALSE)
          lags <- values$lag[values$lag > 0]
          acfs <- values$acf[values$lag > 0]
          lags[which.max(-acfs)]
        } else {
          NA
        }
      }))
    ) %>%
    filter(!is.na(delay)) %>%
    ungroup() %>%
    transmute(
      iso3c = iso3c,
      days_between_doses = if_else(
        delay < 21,
        mean(delay, na.rm = TRUE),
        delay
      )
    )

  #add to data
  dose_df <- left_join(
    dose_df,
    delays_df,
    by = "iso3c"
  )
  #calculate efficacy on each day for each country with waning and dose adjustment
  dose_df <- dose_df %>%
    arrange(iso3c, date_vaccine_change) %>%
    mutate(max_vaccine = first_dose_per_day) %>%
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
  print(plot_waning(adjust_delta = adjust_delta, days_between_doses = mean(
    delays_df$days_between_doses
  )))
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
        paste0(country, ", Assumed days between doses: ", unique(this_country$days_between_doses)),
        size = 20))
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
             max_vaccine = country_df$max_vaccine
           )
         }
)
names(dose_list) <- iso3cs

saveRDS(dose_list,
        "vacc_inputs.Rds"
        )

