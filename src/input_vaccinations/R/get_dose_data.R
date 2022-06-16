#get first doses/second doses
get_dose_ts <- function(owid, who_vacc, who_vacc_meta, iso3cs, date_0){
  #owid data
  owid_ts <- get_dose_ts_owid(owid, iso3cs)
  who_ts <- get_dose_ts_who(who_vacc, who_vacc_meta, iso3cs)
  #prioritise the owid data then who
  iso_in_owid <- unique(owid_ts$iso3c)
  iso_in_who <- setdiff(unique(who_ts$iso3c), iso_in_owid)
  iso_no_data <- setdiff(iso3cs, c(iso_in_owid, iso_in_who))
  #create blank data for the no data countries
  combined_df <- rbind(
    owid_ts %>%
      filter(iso3c %in% iso_in_owid),
    who_ts%>%
      filter(iso3c %in% iso_in_who),
    data.frame(#create blank data for the no data countries
      date = rep(date_0, length(iso_no_data)),
      iso3c = iso_no_data,
      imputed = FALSE,
      first_dose_per_day = 0,
      second_dose_per_day = 0,
      booster_dose_per_day = 0
    )
  )
  #check that we have a point for every date between min  and max for eahc country
  date_check <- unlist(lapply(iso3cs, function(country){
    country_df <- combined_df %>% filter(iso3c == country)
    all(
      seq(min(country_df$date), max(country_df$date), by = 1) %in%
        country_df$date
    )
  }))
  if(any(!date_check)){
    stop(paste0(
      "Did not generate a observation for each day for countries: ",
      paste0(iso3cs[!date_check], collapse = ", "),
      " Please check get_dose_ts function"
      ))
  }

  #expand up to current date
  #first limit data to this date
  combined_df <- combined_df %>%
    filter(date <= date_0) %>%
    filter(cumsum(first_dose_per_day) > 0 | length(first_dose_per_day) == 1)
  #add dates up to current date
  #we'll assume that vacc_per_day stays at its final weekly average, for each type
  combined_df <- combined_df %>% #add indicator for plotting later
    #extend values to date_0
    #add dates up to current date
    complete(date = seq(min(date), date_0, by = 1)) %>%
    left_join( #add averages for vacc_per_day
      combined_df %>%
        group_by(iso3c) %>%
        filter(date >= max(date) - 7) %>%
        summarise(
          first_dose_per_day_week_ave = mean(first_dose_per_day),
          second_dose_per_day_week_ave = mean(second_dose_per_day),
          booster_dose_per_day_week_ave = mean(booster_dose_per_day)
        ),
      by = "iso3c"
    ) %>%
    mutate(
      first_dose_per_day = if_else(
        is.na(first_dose_per_day),
        first_dose_per_day_week_ave,
        first_dose_per_day
      ),
      second_dose_per_day = if_else(
        is.na(second_dose_per_day),
        second_dose_per_day_week_ave,
        second_dose_per_day
      ),
      booster_dose_per_day = if_else(
        is.na(booster_dose_per_day),
        booster_dose_per_day_week_ave,
        booster_dose_per_day
      ),
      imputed = if_else(is.na(imputed),
                        TRUE,
                        imputed)
    ) %>%
    select(!c(first_dose_per_day_week_ave,
              second_dose_per_day_week_ave,
              booster_dose_per_day_week_ave))

  #any other potential cleaning

  return(combined_df)
}
get_dose_ts_owid <- function(owid, iso3cs){
  #owid data
  owid <- owid %>%
    filter(iso3c %in% iso3cs) %>%
    left_join(
      squire::population %>%
        group_by(iso3c) %>%
        summarise(population = sum(n)),
      by = "iso3c"
    ) %>%
    group_by(iso3c) %>%
    mutate(date = as.Date(date))
  #well rely upon cumulative sums of boosted, second dosed and first dosed
  cumulative_data <- owid %>%
    mutate( #derive from per 100 if not present
      first_dose_cum =
        case_when(
          !is.na(people_vaccinated) ~ people_vaccinated,
          !is.na(people_vaccinated_per_hundred) ~ people_vaccinated*population/100,
          TRUE ~ as.numeric(NA)
        ),
      second_dose_cum =
        case_when(
          !is.na(people_fully_vaccinated) ~ people_fully_vaccinated,
          !is.na(people_fully_vaccinated_per_hundred) ~ people_fully_vaccinated_per_hundred*population/100,
          TRUE ~ as.numeric(NA)
        ),
      booster_dose_cum =
        case_when(
          !is.na(total_boosters) ~ total_boosters,
          !is.na(total_boosters_per_hundred) ~ total_boosters_per_hundred*population/100,
          TRUE ~ as.numeric(NA)
        )
    ) %>%
    select(iso3c, date, first_dose_cum, second_dose_cum, booster_dose_cum) %>%
    filter(!is.na(first_dose_cum) | !is.na(second_dose_cum) | !is.na(booster_dose_cum)) %>%
    mutate(imputed = FALSE) %>%
    #add the first date based on the smoothed new vaccination rate
    left_join(
      owid %>%
        select(iso3c, date, new_vaccinations_smoothed) %>%
        arrange(date) %>%
        filter(!is.na(new_vaccinations_smoothed)) %>%
        filter(new_vaccinations_smoothed > 0) %>%
        summarise(start_date = max(min(date), as_date("2020-12-08")))
    ) %>%
    #drop dates past this date
    filter(date >= start_date) %>%
    #extend datas into past using that start date
    complete(date = seq(
      start_date[1],
      max(date),
      1
    )) %>%
    select(!start_date) %>%
    # complete(date = seq(
    #   max(min(date) - min(30, suppressWarnings(min(first_dose_cum, na.rm = TRUE))), min(as_date("2020-12-08"), start_date)),
    #   max(date),
    #   1
    # )) %>%
    arrange(iso3c, date) %>%
    mutate(
      first_dose_cum = fill_missing_doses(first_dose_cum, "first"),
      second_dose_cum = fill_missing_doses(second_dose_cum, "second"),
      booster_dose_cum = fill_missing_doses(booster_dose_cum, "booster"),
      #ensure there are atleast 30 days between doses then a 15 day buildup
      second_dose_cum = linearly_interpolate(case_when(
        seq_along(second_dose_cum) == length(second_dose_cum) ~ second_dose_cum, #catch so we keep the final value
        seq_along(second_dose_cum) < min(which(first_dose_cum > 0)) + 30 ~ 0,
        seq_along(second_dose_cum) < min(which(first_dose_cum > 0)) + 45 ~ as.numeric(NA),
        TRUE ~ second_dose_cum
      )),
      booster_dose_cum = linearly_interpolate(case_when(
        seq_along(booster_dose_cum) == length(booster_dose_cum) ~ booster_dose_cum, #catch so we keep the final value
        seq_along(booster_dose_cum) < min(which(second_dose_cum > 0)) + 30 ~ 0,
        seq_along(booster_dose_cum) < min(which(second_dose_cum > 0)) + 45 ~ as.numeric(NA),
        TRUE ~ booster_dose_cum
      ))
    ) %>%
    #make sure values are consistent, i.e. increasing and booster<second<first
    mutate(
      #if not increasing set offending values to NA and interpolate
      across(
        c(first_dose_cum, second_dose_cum, booster_dose_cum),
        ~linearly_interpolate(if_else(
          diff(c(0, .x)) < 0,
                 as.numeric(NA),
                 .x
                 ))
      ),
      #to check cosistency between variables
      second_dose_cum = make_consistent(second_dose_cum, booster_dose_cum),
      first_dose_cum = make_consistent(first_dose_cum, second_dose_cum)
  )

  #final check for consistency
  if(any(cumulative_data$second_dose_cum > cumulative_data$first_dose_cum) |
     any(cumulative_data$booster_dose_cum > cumulative_data$second_dose_cum)){
    stop("Consistency adjustments have failed!")
  }

  daily_data <- cumulative_data %>%
    mutate(
      first_dose_per_day  = diff(c(0, first_dose_cum)),
      second_dose_per_day  = diff(c(0, second_dose_cum)),
      booster_dose_per_day  = diff(c(0, booster_dose_cum)),
      imputed = if_else(is.na(imputed), TRUE, imputed)
    ) %>%
    select(iso3c, date, first_dose_per_day, second_dose_per_day, booster_dose_per_day, imputed)

  #now we use this data to split the smooth daily dose data into the dose levels
  daily_data %>%
    full_join(
      owid %>%
        select(date, iso3c, new_vaccinations_smoothed) %>%
        filter(!is.na(new_vaccinations_smoothed))
    ) %>%
    arrange(iso3c, date) %>%
    filter(date >= "2020-12-08") %>% #only keep countries with some non-na data
    filter(any(!is.na(new_vaccinations_smoothed)) | any(!is.na(first_dose_per_day))) %>%
    #interpolate the smoothed vaccination rate
    mutate(
      new_vaccinations_smoothed = if_else(
        date == min(date),
        0,
        new_vaccinations_smoothed
      ),
      new_vaccinations_smoothed = linearly_interpolate(new_vaccinations_smoothed)
    ) %>% #split up first/second/booster based on our calculated values
  mutate(
    temp_sum = first_dose_per_day + second_dose_per_day + booster_dose_per_day,
    across(
      c(first_dose_per_day, second_dose_per_day, booster_dose_per_day),
      ~if_else(temp_sum == 0,
               0,
               new_vaccinations_smoothed * .x /(temp_sum)
      )
    ), #need to set values into the future, we'll just assume that the proportions remain the same
    #assume all na are trailing na
    across(
      c(first_dose_per_day, second_dose_per_day, booster_dose_per_day),
      ~if_else(is.na(.x),
               new_vaccinations_smoothed * tail(na.omit(.x), 1)/tail(na.omit(temp_sum), 1),
               .x
      )
    ),
    imputed = if_else(is.na(imputed), TRUE, imputed)
  ) %>%
  select(iso3c, date, first_dose_per_day, second_dose_per_day, booster_dose_per_day, imputed) %>%
  #remove leading zeros
  filter(cumsum(first_dose_per_day) > 0)

  #probably needs some smoothing at some point since we no longer use the OWID smoothed data
}
get_dose_ts_who <- function(who_vacc, who_vacc_meta, iso3cs){
  #have to impute a time series, only use if owid not available
  #get any useful info out of meta df
  meta_df <- who_vacc_meta %>%
    filter(iso3c %in% iso3cs) %>%
    group_by(iso3c) %>%
    mutate(
      START_DATE = as.Date(START_DATE),
      END_DATE = as.Date(END_DATE),
      AUTHORIZATION_DATE = as.Date(AUTHORIZATION_DATE)
    ) %>%
    summarise(
      start_date = as.Date(ifelse(
        any(!is.na(START_DATE)),
        min(START_DATE, na.rm = TRUE),
        NA
      ), origin = "1970-01-01"),
      end_date = as.Date(ifelse(
        any(!is.na(END_DATE)),
        max(END_DATE, na.rm = TRUE),
        NA
      ), origin = "1970-01-01"),
      aut = as.Date(ifelse(
        any(!is.na(AUTHORIZATION_DATE)),
        min(AUTHORIZATION_DATE, na.rm = TRUE),
        NA
      ), origin = "1970-01-01")
    ) %>%
    mutate(
      start_date = if_else(
        is.na(start_date),
        aut,
        start_date
      )
    ) %>%
    select(iso3c, start_date, end_date) %>%
    #start date can't be before 8th december
    mutate(
      start_date = if_else(start_date < as.Date("2020-12-08"),
                           as.Date("2020-12-08"),
                           start_date)
    )

  who_comb <- full_join(
    who_vacc %>% filter(iso3c %in% iso3cs),
    meta_df,
    by = "iso3c"
  )

  #get the data we'll use to set up vaccinations
  who_comb <- who_comb %>%
    mutate(
      start_date = if_else( #use meta start date if first date is missing
        is.na(as.Date(FIRST_VACCINE_DATE)),
        start_date,
        as.Date(FIRST_VACCINE_DATE)
      ),
      end_date = if_else( #use meta end date if first date is missing
        is.na(as.Date(DATE_UPDATED)),
        end_date,
        as.Date(DATE_UPDATED)
      ),
      single_dose_plus = as.numeric(PERSONS_VACCINATED_1PLUS_DOSE),
      two_dose = as.numeric(PERSONS_FULLY_VACCINATED),
      vaccinations_total = case_when(
        !is.na(as.numeric(TOTAL_VACCINATIONS)) ~ as.numeric(TOTAL_VACCINATIONS),
        (!is.na(single_dose_plus)) & (!is.na(two_dose)) ~ single_dose_plus + two_dose,
        !is.na(single_dose_plus) ~ single_dose_plus,
        !is.na(two_dose) ~ two_dose*2,
        TRUE ~ as.numeric(NA)
      ),
      #update single/second dose if missing
      single_dose_plus = if_else(
        is.na(single_dose_plus),
        vaccinations_total - two_dose,
        single_dose_plus
      ),
      two_dose = if_else(
        is.na(two_dose),
        vaccinations_total - single_dose_plus,
        two_dose
      )
    ) %>%
    select(iso3c, start_date, end_date, vaccinations_total, single_dose_plus, two_dose,
           recieved_single_dose_vaccines)
    #for countries with more two dose than single dose we assume the difference
  # are single dose vaccines and swap them over to better represent one dose efficacy
  who_comb <- who_comb %>%
    mutate(
      recieved_single_dose_vaccines = if_else(
        is.na(recieved_single_dose_vaccines),
        FALSE,
        recieved_single_dose_vaccines
      ),
      single_dose_vacc = if_else(
        #recieved_single_dose_vaccines &
          (single_dose_plus < two_dose),
        two_dose - single_dose_plus,
        0
      ),
      #adjust counts
      single_dose_plus = single_dose_plus + single_dose_vacc,
      two_dose = two_dose - single_dose_vacc
    ) %>% #use the two values to calculate the dose ratio
    mutate(
      dose_ratio_final = two_dose/single_dose_plus
    ) %>%
    select(iso3c, start_date, end_date, vaccinations_total, dose_ratio_final)

  #impute any missing values using regional medians etc
  who_comb <- who_comb %>%
    ungroup() %>%
    mutate(
      #get region/sub_region
      region = countrycode(iso3c, origin = "iso3c", destination = "region"),
      sub_region = countrycode(iso3c, origin = "iso3c", destination = "region23")
    ) %>%
    group_by(sub_region) %>%
    mutate(
      across(
        .cols = c(start_date, end_date, vaccinations_total, dose_ratio_final),
        ~if_else(
          is.na(.x) | is.infinite(.x),
          median(.x, na.rm = TRUE),
          .x
        )
      )
    ) %>%
    group_by(region) %>%
    mutate(
      across(
        .cols = c(start_date, end_date, vaccinations_total, dose_ratio_final),
        ~if_else(
          is.na(.x) | is.infinite(.x),
          median(.x, na.rm = TRUE),
          .x
        )
      )
    ) %>%
    ungroup() %>%
    mutate(
      across(
        .cols = c(start_date, end_date, vaccinations_total, dose_ratio_final),
        ~if_else(
          is.na(.x) | is.infinite(.x),
          median(.x, na.rm = TRUE),
          .x
        )
      )
    ) %>%
    group_by(iso3c)

  #generate TS
  #assume that first 21 days are a build up period (if possible), then vaccinations
  #per day remain steady across that period
  #assume dose ratio is 0 for the first 30 days and then build up to the value
  #for the next 21 days then remain at the same rate
  ts_df <- do.call(
    rbind,
    lapply(1:nrow(who_comb), function(x){
    dates <- seq(who_comb$start_date[x], who_comb$end_date[x], by = 1)
    #calculate the vaccines per day
    if(length(dates) <= 21){
      v_r_b <- who_comb$vaccinations_total[x]/sum(1:length(dates))
      vacc_per_day <- v_r_b*seq(1, length(dates))
    } else {
      v_r_b <- who_comb$vaccinations_total[x]/(sum(1:21) + 21*(length(dates)-21))
      v_r_pb <- 21*v_r_b
      vacc_per_day <- c(v_r_b*seq(1, 21), rep(v_r_pb, length(dates)-21))
    }
    #calculate dose ratio over time
    if(length(dates) <= 30){
      dr_r_0 <- who_comb$dose_ratio_final[x]/length(dates)
      dose_ratio <- c(
        dr_r_0*seq(1, length(dates))
      )
    } else if (length(dates) <= 30 + 21){
      end_0 <- floor((30/(30+21)) * length(dates))
      dr_r_0 <- 0
      dr_r_b <- who_comb$dose_ratio_final[x]/(length(dates) - end_0)
      dose_ratio <- c(
        rep(dr_r_0, end_0),
        dr_r_b*seq(1,length(dates) - end_0)
      )
    } else {
      dr_r_0 <- 0
      dr_r_b <- who_comb$dose_ratio_final[x]/21
      dr_r_pb <- who_comb$dose_ratio_final[x]
      dose_ratio <- c(
        rep(dr_r_0, 30),
        dr_r_b * seq(1, 21),
        rep(dr_r_pb, length(dates) - 21 - 30)
      )
    }
    data.frame(
      iso3c = who_comb$iso3c[x],
      date = dates,
      vacc_per_day = vacc_per_day,
      imputed = TRUE,
      dose_ratio
    )
  })
  )

  #calculate the first/second doses each day (should have no issues)
  ts_df <- ts_df %>%
    group_by(iso3c) %>%
    mutate(
      first_dose_per_day = diff(c(0,
                                  cumsum(vacc_per_day)/(1 + dose_ratio)
                                  )),
      second_dose_per_day = vacc_per_day - first_dose_per_day,
      #possible less than 0s due to round errors so we set this to 0
      second_dose_per_day = if_else(
        abs(second_dose_per_day) < 0.1,
        0,
        second_dose_per_day
      ),
      #assume no boosters
      booster_dose_per_day = 0
    ) %>%
    select(!c(dose_ratio, vacc_per_day))
  #final check
  # problem <- test_for_errors_in_dose_ratio(ts_df, dose_ratio)
  # if(length(problem) > 1){
  #   stop(paste0("Drop with calculating dose ratios in the following countires: ",
  #               paste0(problem, collapse = ", ")
  #               ))
  # }
  return(ts_df)
}
linearly_interpolate <- function(vector){
  if_else(
    is.na({vector}),
    approx(1:length({vector}), {vector}, 1:length({vector}), yright = NA)$y,
    {vector}
  )
}
test_for_errors_in_dose_ratio <- function(df, dose_ratio){
  df %>%
    mutate(
      first_dose = diff(c(0, cumsum(vacc_per_day)/({dose_ratio} + 1))),
      problematic = first_dose < -0.1 | first_dose > vacc_per_day + 0.1
    ) %>%
    filter(problematic) %>%
    pull(iso3c) %>%
    unique()
}
fill_missing_doses <- function(vector, dose_type = "first"){
  #new method to fill missing, less need for consistency as this is handled in the model now
  vector <- if_else(
    #if all NA assume none
    rep(all(is.na(vector )), length(vector )),
    0,
    vector
  )
  vector <- if_else(
    #if all non NA are 0 assume 0
    rep(sum(vector, na.rm = TRUE) == 0, length(vector)),
    0,
    vector
  )
  if(any(vector > 0)){


    #if first dose type we do nothing to it so that it treats the start of the vector as when the dosing started
    #if second dose type we assume it starts 30 days from the start of the vector
    #if its booster we just assume 30 days from when first reported

    if(dose_type == "booster"){

      first_dose_loc <- min(which(vector > 0))
      #catch if we don't actaully have 30 days
      if(first_dose_loc > 30) {
        first_vacc <- vector[first_dose_loc[1]]
        extension <- min(30, first_vacc[1])

        vector[first_dose_loc - extension + seq_len(extension)] <- seq(0, first_vacc, length.out = extension)
      }
    } else if(dose_type == "second"){
      vector[30] <- 0
    }
    #set first value to 0
    vector[1] <- 0
    vector <- linearly_interpolate(vector)
    #sort final values
    if(any(is.na(vector))){
      #assume continous to grow at same rate as previous two weeks
      growth_rate <- mean(diff(tail(na.omit(vector), 14)))
      last_non_na <- if_else(any(is.na(vector)),
                             suppressWarnings(min(which(is.na(vector))) - 1),
                             as.numeric(length(vector))
      )
      vector[is.na(vector)] <- vector[min(which(is.na(vector))) - 1] +
        growth_rate[1] * seq_len(sum(is.na(vector)))
    }
  }
  vector
}
make_consistent <- function(first_cum_ts, second_cum_ts){
  #make daily differences
  first_ts <- diff(c(0, first_cum_ts))
  second_ts <- diff(c(0, second_cum_ts))
  for(t in seq_along(second_ts)[-1]){
    if(cumsum(first_ts)[t-1] < cumsum(second_ts)[t]){
      #have issue scale up first_ts in the previous day to match
      difference <- cumsum(second_ts)[t] - cumsum(first_ts)[t-1]
      first_ts[t-1] <- first_ts[t-1] + difference
    }
  }
  cumsum(first_ts)
}
