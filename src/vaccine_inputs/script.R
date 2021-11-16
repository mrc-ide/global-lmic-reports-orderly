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

#set prioritization + coverage matrix
#standard strategy, might make particular later
strategy <- "HCW, Elderly and High-Risk"
#set vaccine uptake to be 80% or higher if it is in the data
vaccine_uptake <- get_vaccine_uptake(
  iso3cs = iso3cs,
  dose_df = dose_df,
  default_uptake = 0.8
  )
#get the matrices
vaccine_coverage_mats <- get_coverage_mats(
  iso3c = iso3cs,
  strategy = strategy,
  vaccine_uptake = vaccine_uptake
)


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
  #here we assume that doses go to the person vaccinated the earliest
  delays_df <- data.frame(iso3c = iso3cs) %>%
    rowwise() %>%
    mutate(
      days_between_doses = {
        country_df <- dose_df %>%
          rename(country = iso3c) %>%
          filter(country == iso3c)
        if(nrow(country_df) == 1 | sum(country_df$second_dose_per_day) == 0){
          NA
        } else {
          #set up parameters
          delays_df <- data.frame(
            days_between_doses = 1:nrow(country_df),
            frequency = 0
          )
          first_doses_available <- country_df$first_dose_per_day
          #now work though for each day
          for(t in 1:nrow(country_df)){
            second_doses_to_assign <- country_df$second_dose_per_day[t]
            if(second_doses_to_assign>0){
              assigning <- TRUE
              i <- min(which(first_doses_available>0))
              while(assigning){
                assigned <- max(
                  (second_doses_to_assign - first_doses_available[i]),
                  second_doses_to_assign
                )
                #update the delay df
                delays_df[delays_df$days_between_doses == (t - i), "frequency"] <-
                  assigned
                #update parameters
                first_doses_available[i] <- first_doses_available[i] - assigned
                second_doses_to_assign <- second_doses_to_assign - assigned
                i <- i+1
                #check if we have used up all doses to assign
                if(second_doses_to_assign == 0){
                  assigning <- FALSE
                }
              }
            }
          }
          #calculate the median delay
          middle_delay <- (sum(delays_df$frequency)+1)/2
          delays_df$days_between_doses[min(which((cumsum(delays_df$frequency)>middle_delay)))]
        }
      }
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
    delays_df$days_between_doses, na.rm = TRUE
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
             max_vaccine = country_df$max_vaccine,
             vaccine_coverage_mat = vaccine_coverage_mats[[country]],
             strategy = strategy
           )
         }
)
names(dose_list) <- iso3cs

saveRDS(dose_list,
        "vacc_inputs.Rds"
        )

