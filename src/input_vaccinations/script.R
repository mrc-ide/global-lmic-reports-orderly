date_0 <- as.Date(date)

iso3cs <- unique(squire::population$iso3c)

### Get the vaccine data from online
#sort out dates from before date_0 and group data by country

vacc_types <- import_vaccine_agreements()

vdm <- import_vaccine_doses_by_manufacturer()

who_vacc <- import_who_vaccination()

who_vacc_meta <- import_who_vaccination_meta()

owid <- import_owid()

#determine the vaccine platforms used in each country
platforms_df <- get_platforms(iso3cs, vacc_types, vdm, who_vacc, who_vacc_meta)

## Modifications for single dose Vaccines
#we will treat as only a single dose, more suitable
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

#get our doses over time
dose_df <- get_dose_ts(owid, who_vacc, who_vacc_meta, iso3cs, date_0)

dose_df <- dose_df %>%
  rename(date_vaccine_change = date) %>%
  filter(date_vaccine_change <= date_0)

#set prioritization + coverage matrix
#standard strategy, might make particular later
strategy <- "HCW, Elderly and High-Risk"
#set vaccine uptake to be 80% or higher if it is in the data
vaccine_uptake <- get_vaccine_uptake(
  iso3cs = iso3cs,
  dose_df = dose_df,
  default_uptake = 0.8,
  strategy = strategy
  )
#get the matrices
vaccine_coverage_mats <- get_coverage_mats(
  iso3c = iso3cs,
  strategy = strategy,
  vaccine_uptake = vaccine_uptake
)

#into list format
names(iso3cs) <- iso3cs
dose_list <-
  map(iso3cs,
         function(country){
           country_df <- dose_df %>%
             filter(iso3c == country)
           platforms <- platforms_df %>%
             filter(iso3c == country) %>%
             select(!iso3c)
           #estimate the second dose delay
           scaled_first_doses <- cumsum(country_df$first_dose_per_day)/sum(country_df$first_dose_per_day)
           scaled_second_doses <- cumsum(country_df$second_dose_per_day)/sum(country_df$second_dose_per_day)
           delays <- seq(30, 120, by = 1)
           errs <- map_dbl(delays, function(delay){
             mean((na.omit(scaled_first_doses - lead(scaled_second_doses, delay)))^2)
           })
           second_dose_delay <- delays[which.min(errs)]
           list(
             date_vaccine_change = country_df$date_vaccine_change,
             primary_doses = country_df$first_dose_per_day,
             second_dose_delay = second_dose_delay,
             booster_doses = country_df$booster_dose_per_day,
             vaccine_coverage_mat = vaccine_coverage_mats[[country]],
             strategy = strategy,
             platforms = platforms
           )
         }
)

saveRDS(dose_list,
        "vacc_inputs.Rds"
        )

