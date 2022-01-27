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
dur_V <- 446

### Derive effective vaccine efficacies over time

iso3cs <- unique(squire::population$iso3c)

#get our doses over time
dose_df <- get_dose_ts(owid, who_vacc, who_vacc_meta, iso3cs, date_0)

## Apply Delta Adjustment if needed
#just set vaccine efficacies to constant values
dose_df <- dose_df %>%
  rename(date_vaccine_change = date)

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

vaccine_efficacies <- list(
  ve_i_low = c(0.6, 0.5, 0.7),
  ve_i_high = c(0.8, 0.7, 0.85),
  ve_d_low = c(0.8, 0.7, 0.85),
  ve_d_high = c(0.98, 0.95, 0.99),
  ve_i_low_d = c(0.224, 0.1, 0.3),
  ve_i_high_d = c(0.646, 0.55, 0.75),
  ve_d_low_d = c(0.75, 0.65, 0.8),
  ve_d_high_d = c(0.94, 0.9, 0.98)
)

#into list format
dose_list <-
  lapply(iso3cs,
         function(country){
           country_df <- dose_df %>%
             filter(iso3c == country)
           list(
             dur_vaccine_delay = dur_vaccine_delay,
             dur_V = dur_V,
             vaccine_efficacies = vaccine_efficacies,
             date_vaccine_change = country_df$date_vaccine_change,
             max_vaccine = country_df$first_dose_per_day,
             date_vaccine_efficacy = country_df$date_vaccine_change,
             dose_ratio = country_df$dose_ratio,
             rel_infectiousness_vaccinated = rep(rel_infectiousness_vaccinated, 17),
             vaccine_coverage_mat = vaccine_coverage_mats[[country]],
             strategy = strategy
           )
         }
)
names(dose_list) <- iso3cs

saveRDS(dose_list,
        "vacc_inputs.rds"
        )

