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
        select(iso3c, single_dose_percentage),
      by = "iso3c"
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
        select(iso3c, single_dose_percentage_1),
      by = "iso3c"
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

single_dose_adjust_owid <- function(owid, single_dose_df){
  #only do those with single doses
  modified <- owid %>%
    group_by(iso3c) %>%
    filter(iso3c %in% single_dose_df$iso3c & ( #also only those with vaccine data, so it doesn't throw errors
      any(!is.na(people_fully_vaccinated)) |
        any(!is.na(people_fully_vaccinated_per_hundred)) |
        any(!is.na(total_vaccinations)) |
        any(!is.na(total_vaccinations_per_hundred))
    )) %>%
    left_join(single_dose_df,
              by = "iso3c") %>%
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
  #add back in data
  new <- rbind(
      owid %>%
        filter(!(iso3c %in% modified$iso3c)) %>%
        mutate(
          recieved_single_dose_vaccines = if_else(
            iso3c %in% single_dose_df$iso3c,
            TRUE,
            as.logical(NA)
          )
        ),
      modified
    ) %>%
    arrange(
      iso3c, date
    )
}

single_dose_adjust_who_vacc <- function(who_vacc, single_dose_df){
  who_vacc %>%
    left_join(single_dose_df,
              by = "iso3c") %>%
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
}
