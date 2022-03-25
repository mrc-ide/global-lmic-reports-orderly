##Load in fits
vacc_inputs_raw <- readRDS("vacc_inputs.rds")

##Adjust values for excess fitting
vacc_inputs_old_format <- map(vacc_inputs_raw, ~
                                list(
                                  dur_V = 446,
                                  date_vaccine_change = .x$date_vaccine_change,
                                  date_vaccine_efficacy = .x$date_vaccine_change,
                                  max_vaccine = .x$first_doses,
                                  dose_ratio = if_else(.x$first_doses == 0, 0, cumsum(.x$second_doses)/cumsum(.x$first_doses)),
                                  vaccine_coverage_mat = .x$vaccine_coverage_mat,
                                  strategy = .x$strategy,
                                  dur_vaccine_delay = 14,
                                  rel_infectiousness_vaccinated = rep(0.5, 17)
                                )
                          )

##Get dominant vaccine type and others
vaccines_df  <- get_dominant_vaccines(names(vacc_inputs_old_format))

##Define vaccine efficacies (temp)
vaccine_efficacy <- read_csv("vaccine_efficacy_groups.csv")

#Infer missing values
vaccine_efficacy <- vaccine_efficacy %>%
  rbind(
    #for our framework it suffices to assume partial == full for single does vaccines
    vaccine_efficacy %>%
      filter(vaccine_type == "Johnson&Johnson") %>%
      mutate(dose = "Partial")
  )

##Assign to countries
iso3cs <- names(vacc_inputs_old_format)
names(iso3cs) <- iso3cs
vaccine_outputs_multi_vaccine <- map(iso3cs, function(iso){
  inputs <- vacc_inputs_old_format[[iso]]

  #need to get the prior value and the range
  vaccines <- vaccines_df %>%
    filter(iso3c == iso) %>%
    select(!iso3c)
  domiant_values <- vaccine_efficacy %>%
    filter(vaccine_type == vaccines$dominant) %>%
    rename(central = efficacy) %>%
    select(! vaccine_type)
  ranges <- vaccine_efficacy %>%
    filter(vaccine_type  %in% (vaccines %>%
                select(!dominant) %>%
                select(where(~.x)) %>%
                names()
             )
    ) %>%
    group_by(variant, dose, endpoint) %>% #add some variance so that we still have uncertainty in estimates when domiant is a bound or only 1 vaccine type
    summarise(lower = min(efficacy) * 0.8,
              upper = min(max(efficacy) * 1.1, 1),
              .groups = "drop") %>%
    full_join(domiant_values, by = c("variant", "dose", "endpoint")) %>%
    #if dominant is NA we just set the central value to the midpoint of the other, that way it is just uniformly distributed
    mutate(
      cental = if_else(
        is.na(central),
        (lower + upper)/2,
        central
      )
    )

  #assume we are using ves method
  inputs$vaccine_efficacies <- list(
    ve_i_low = ranges %>% filter(dose == "Partial", endpoint == "Infection", variant == "Wild") %>% select(central, lower, upper) %>% as.numeric(),
    ve_i_high = ranges %>% filter(dose == "Full", endpoint == "Infection", variant == "Wild") %>% select(central, lower, upper) %>% as.numeric(),
    ve_d_low = ranges %>% filter(dose == "Partial", endpoint == "Hospitalisation", variant == "Wild") %>% select(central, lower, upper) %>% as.numeric(),
    ve_d_high = ranges %>% filter(dose == "Full", endpoint == "Hospitalisation", variant == "Wild") %>% select(central, lower, upper) %>% as.numeric(),
    ve_i_low_d = ranges %>% filter(dose == "Partial", endpoint == "Infection", variant == "Delta") %>% select(central, lower, upper) %>% as.numeric(),
    ve_i_high_d = ranges %>% filter(dose == "Full", endpoint == "Infection", variant == "Delta") %>% select(central, lower, upper) %>% as.numeric(),
    ve_d_low_d = ranges %>% filter(dose == "Partial", endpoint == "Hospitalisation", variant == "Delta") %>% select(central, lower, upper) %>% as.numeric(),
    ve_d_high_d = ranges %>% filter(dose == "Full", endpoint == "Hospitalisation", variant == "Delta") %>% select(central, lower, upper) %>% as.numeric()
  )

  inputs
})

##Save object
saveRDS(vaccine_outputs_multi_vaccine, "vacc_inputs_multi.rds")
