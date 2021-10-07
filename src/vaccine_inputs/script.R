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
    PERSONS_FULLY_VACCINATED
  ),
  PERSONS_FULLY_VACCINATED_PER100 = if_else(
    (!is.na(PERSONS_FULLY_VACCINATED_PER100)) & (!is.na(single_doses)),
    PERSONS_FULLY_VACCINATED_PER100*(1-single_doses/TOTAL_VACCINATIONS),
    PERSONS_FULLY_VACCINATED_PER100
  )) %>%
  mutate(
    PERSONS_VACCINATED_1PLUS_DOSE = if_else(
      identical(single_dose_percentage, 1) & (!is.na(PERSONS_VACCINATED_1PLUS_DOSE)),
      TOTAL_VACCINATIONS,
      PERSONS_VACCINATED_1PLUS_DOSE
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

if(!waning){
  #use the standard method of waning
  dur_V <- 5000
  #extend data
  dose_df <- dose_df %>% #add indicator for plotting later
    #extend values to date_0
    mutate(imputed = FALSE) %>%
    #add dates up to current date
    complete(date_vaccine_change = seq(min(date_vaccine_change), date_0, by = 1)) %>%
    #assume max_vaccine and dose-ratio remain the same
    fill(max_vaccine, dose_ratio, vaccine_efficacy_infection,
         vaccine_efficacy_disease, .direction = "down") %>%
    mutate(imputed = if_else(
      is.na(imputed),
      TRUE,
      FALSE
    ))
  if(adjust_delta){
    #add delta data
    dose_df <- add_delta_characteristics(dose_df) %>%
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
             shift_start, shift_end, imputed)
  }
  #if not waning we just take as is
} else {
  ## Recalculate when people are getting first and second doses

  #we take max_vaccine as people getting first doses

  #confirm that vaccine change every day
  if(length(
    unique(unlist(lapply(iso3cs, function(x){unique(diff(filter(dose_df, iso3c==x)$date_vaccine_change))})))
  )>1){
    stop("Code needs adjusting as vaccine changes less than everyday")
  }

  #extend data
  dose_df <- dose_df %>% #add indicator for plotting later
    #extend values to date_0
    mutate(imputed = FALSE) %>%
    #add dates up to current date
    complete(date_vaccine_change = seq(min(date_vaccine_change), date_0, by = 1)) %>%
    #assume max_vaccine and dose-ratio remain the same
    fill(max_vaccine, dose_ratio, .direction = "down") %>%
    mutate(imputed = if_else(
      is.na(imputed),
      TRUE,
      FALSE
    ))

  dose_df <- dose_df %>%
    mutate(
      people_with_atleast_one_dose = cumsum(max_vaccine),
      #from this we calculate people getting two doses each day
      people_with_two_doses = people_with_atleast_one_dose*dose_ratio,
      #from this we get the number of new second doses each day
      second_doses = c(0, diff(people_with_two_doses))
    )

  #countries to check
  # "ALB" "AUS" "BHS" "BFA" "CHN" "MNG" "NPL" "LCA" "SAU" "VNM"

  #dose_df %>% filter(iso3c == "ALB") %>% View()

  #for now we set negatives to 0
  dose_df <- dose_df %>%
    mutate(
      second_doses = if_else(second_doses <0, 0, second_doses)
    ) %>%
    select(!c(people_with_atleast_one_dose, people_with_two_doses))

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

  #calculate efficacy on each day for each country with waning and dose adjustment
  dose_df <- dose_df %>%
    arrange(iso3c, date_vaccine_change) %>%
    calculate_waning_eff(dates = "date_vaccine_change",
                         first_doses = "max_vaccine",
                         second_doses = "second_doses",
                         efficacy_infection_first = "ve_i_low",
                         efficacy_infection_second = "ve_i_high",
                         efficacy_disease_first = "ve_d_low",
                         efficacy_disease_second = "ve_d_high",
                         dur_vaccine_delay = dur_vaccine_delay,
                         countries = "iso3c",
                         diagnostic = TRUE)
}

### Calibration plot
dir.create("calibration", showWarnings = FALSE)
pdf("calibration/plot.pdf", paper = "a4r")
if(waning){
  #plot efficacy curves for our first and second doses
  plot_waning(adjust_delta = adjust_delta)
}

#plot for each country
for(country in iso3cs){
  this_country <- dose_df %>%
    filter(iso3c == country) %>%
    rename(Date = date_vaccine_change)
  if(nrow(this_country) > 1){
    #if there is data

    #plot of new doses each day
    dose_plot <- plot_doses(this_country, waning)
    if(waning){
      dose_comp_plot <- plot_dose_comp(this_country)
      dose_plot <- suppressWarnings(ggarrange(
        dose_plot,
        dose_comp_plot,
        common.legend = TRUE,
        ncol = 1
      ))
    }
    #plot dose ratio
    dose_ratio_plot <- plot_dose_ratio(this_country, waning)
    #plot mean first does
    if(waning){
      mean_first_plot <- plot_mean_first_date(this_country)
      dose_ratio_plot <- suppressWarnings(ggarrange(
        dose_ratio_plot,
        mean_first_plot,
        common.legend = TRUE,
        ncol = 1
      ))
    }
    #final efficacy plots
    final_eff_plot <- plot_efficacy(this_country, adjust_delta)
    #arrange into a layout
    suppressWarnings(
      print(ggarrange(
        as_ggplot(text_grob(
          country,size = 20)),
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
dev.off()

### Export

#reduce to dates that relevent <date_0 and variables that we need
dose_df <- dose_df %>%
  filter(date_vaccine_change<=date_0) %>%
  select(iso3c, date_vaccine_change, max_vaccine, vaccine_efficacy_disease,
         vaccine_efficacy_infection) %>% #adjust disease efficacy for infection efficacy
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
             vaccine_efficacy_disease = rep(country_df$vaccine_efficacy_disease,
                                            17),
             vaccine_efficacy_infection = rep(country_df$vaccine_efficacy_infection,
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

