get_single_dose_data <- function(
  owid, vacc_types, vdm, who_vacc, who_vacc_meta, single_dose_vaccines
){
  vdm %>% #get info on doses of any single dose vaccines from vdm
    filter(date == max(date)) %>%
    group_by(iso3c) %>%
    mutate(single_dose = if_else(
      vaccine %in% single_dose_vaccines,
      1,
      0
    )) %>%
    summarise(
      single_doses = sum(total_vaccinations*single_dose)
    ) %>%
    mutate(single_doses = if_else(
      single_doses == 0,
      as.numeric(NA),
      single_doses
    )) %>%
    full_join( #get info on what single dose vaccines are used, if only single we know
      #all vaccine doses are single doses
      vacc_types %>%
        select(iso3c, vaccine_types) %>%
        rowwise() %>%
        filter(any(single_dose_vaccines %in% vaccine_types)) %>%
        mutate(
          single_dose_percentage = if_else(
            all(vaccine_types %in% single_dose_vaccines),
            1,
            as.numeric(NA)
          )
        ) %>%
        select(iso3c, single_dose_percentage)
    ) %>%
    full_join(
      who_vacc %>%
        rowwise() %>%
        filter(any(VACCINES_USED %in% single_dose_vaccines)) %>%
        mutate(
          single_dose_percentage_1 = if_else(
            all(VACCINES_USED %in% single_dose_vaccines),
            1,
            as.numeric(NA)
          )
        ) %>%
        select(iso3c, single_dose_percentage_1)
    ) %>%
    mutate(
      single_dose_percentage = if_else(
        is.na(single_dose_percentage),
        single_dose_percentage_1,
        single_dose_percentage
      ),
      recieved_single_dose_vaccines = TRUE #indicator for checking later
    ) %>%
    select(iso3c, single_doses, single_dose_percentage, recieved_single_dose_vaccines)
}
