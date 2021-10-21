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
      vacc_per_day = 0,
      dose_ratio = 0,
      imputed = FALSE,
      first_dose_per_day = 0,
      second_dose_per_day = 0
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
    filter(cumsum(vacc_per_day) > 0 | length(vacc_per_day) == 1)
  #add dates up to current date
  #we'll assume that vacc_per_day stays at its final weekly average
  #percentage of vacciantions that are second doses will be its final weekly average
  combined_df <- combined_df %>% #add indicator for plotting later
    #extend values to date_0
    mutate(imputed = FALSE)%>%
    #add dates up to current date
    complete(date = seq(min(date), date_0, by = 1)) %>%
    left_join( #add averages for vacc_per_day
      combined_df %>%
        group_by(iso3c) %>%
        filter(date >= max(date) - 7) %>%
        summarise(
          vacc_per_day_week_ave = mean(vacc_per_day)
        ),
      by = "iso3c"
    ) %>%
    left_join( #add averages for percentage second dose
      combined_df %>%
        group_by(iso3c) %>%
        filter(!is.na(second_dose_per_day)) %>%
        filter(date >= max(date) - 7) %>%
        summarise(
          percentage_second_dose_week_ave = mean(second_dose_per_day/vacc_per_day)
        ),
      by = "iso3c"
    ) %>%
    mutate(
      vacc_per_day = if_else(
        is.na(vacc_per_day),
        vacc_per_day_week_ave,
        vacc_per_day
      )
    )
  #now we have to iteratively calculate the first doses, limiting the percentage
  #second doses so that dose ratio does not exceed 1
  missing_df <- combined_df %>% summarise(
    missing = sum(is.na(dose_ratio))
  )
  while(sum(missing_df$missing) > 0){
    combined_df <- mutate(combined_df,
                          percentage_second_dose =
                            pmin(
                              percentage_second_dose_week_ave,
                              (vacc_per_day + lag(cumsum(first_dose_per_day)
                                                  - cumsum(second_dose_per_day), 1))/(2*vacc_per_day)
                            ),
                          second_dose_per_day = if_else(
                            is.na(
                              second_dose_per_day
                            ),
                            vacc_per_day*percentage_second_dose,
                            second_dose_per_day
                          ),
                          first_dose_per_day = if_else(
                            is.na(first_dose_per_day),
                            vacc_per_day - second_dose_per_day,
                            first_dose_per_day
                          ),
                          dose_ratio = if_else(
                            is.na(dose_ratio),
                            cumsum(second_dose_per_day)/cumsum(first_dose_per_day),
                            dose_ratio
                          )
    ) %>%
      select(!percentage_second_dose)


    missing_df <- combined_df %>% summarise(
      missing = sum(is.na(dose_ratio))
    )

    #message(sum(missing_df$missing))
  }
  #final changes
  combined_df <- combined_df %>%
    mutate(
      imputed = if_else(
        is.na(imputed),
        TRUE,
        FALSE
      )
    ) %>%
    select(!c(vacc_per_day_week_ave, percentage_second_dose_week_ave))

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
  owid_m <- owid %>%
    mutate( #derive from per 100 if not present
      vaccinations_cum = if_else(
        rep(all(is.na(total_vaccinations)), length(total_vaccinations)),
        (total_vaccinations_per_hundred*population)/100,
        total_vaccinations
      ),
      #use smoothed if possible
      vaccinations =
        case_when(
          !is.na(new_vaccinations_smoothed) ~ new_vaccinations_smoothed,
          !is.na(new_vaccinations_smoothed_per_million) ~ (new_vaccinations_smoothed_per_million*population)/10^6,
          !is.na(new_vaccinations) ~ new_vaccinations,
          TRUE ~ as.numeric(NA)
        )
    ) %>%
    mutate( #derive from per 100 if not present
      first_dose_cum =
        case_when(
          !is.na(people_vaccinated) ~ people_vaccinated,
          !is.na(people_vaccinated_per_hundred) ~ people_vaccinated*population/100,
          TRUE ~ as.numeric(NA)
        )
    ) %>%
    mutate( #derive from per 100 if not present
      second_dose_cum =
        case_when(
          !is.na(people_fully_vaccinated) ~ people_fully_vaccinated,
          !is.na(people_fully_vaccinated_per_hundred) ~ people_fully_vaccinated_per_hundred*population/100,
          TRUE ~ as.numeric(NA)
        )
    ) %>%
    select(iso3c, date, vaccinations, vaccinations_cum, first_dose_cum, second_dose_cum)
  #extend data so there is on entry each day
  owid_m <- complete(owid_m, date = seq(min(date), max(date), by = 1))
  #filter data and entries
  owid_m <- owid_m %>%
    mutate(aux = !is.na(vaccinations) |
             !is.na(vaccinations_cum) |
             !is.na(first_dose_cum) |
             !is.na(second_dose_cum)) %>% #only keep countries with some data
    filter(any(aux)) %>% #remove trailing zeros
    filter(rev(cumsum(rev(aux))) > 0) %>% #remove leading zeros
    filter(cumsum(aux) > 0) %>%
    select(!aux)
  #if cumulative is greater than 0 at start add two weeks and have it build up
  countries_to_extend <- owid_m %>%
    mutate(start_vacc = case_when(
      !is.na(vaccinations_cum) ~ vaccinations_cum,
      !is.na(first_dose_cum) ~ first_dose_cum,
      !is.na(second_dose_cum) ~ second_dose_cum,
      TRUE ~ as.numeric(NA)
    )) %>%
    filter(!is.na(start_vacc)) %>%
    filter(date == min(date)) %>%
    filter(start_vacc > 0) %>%
    select(iso3c, date, start_vacc) %>% #how many days to extend back
    mutate(extension = min(21, start_vacc))
  #extended data
  owid_m <- owid_m %>%
    mutate(imputed = FALSE) %>%
    rbind(
      do.call(
        rbind,
        lapply(
          1:nrow(countries_to_extend), function(x){
            x <- countries_to_extend[x,]
            df <- data.frame(
              iso3c = x$iso3c,
              date = seq(x$date - x$extension, x$date, by = 1),
              vaccinations = as.numeric(NA),
              vaccinations_cum = seq(0, x$start_vacc, length.out = x$extension + 1),
              first_dose_cum = as.numeric(NA),
              second_dose_cum = as.numeric(NA),
              imputed = TRUE
            )
            df[-nrow(df),]
          }
        )
      )
    ) %>%
    arrange(iso3c, date) %>%
    mutate(#fill in vacciantions
      vaccinations = if_else(
        imputed,
        diff(c(0, vaccinations_cum)),
        vaccinations
      )
    )
  #figure out whether to use daily or cumulative
  #which has more data
  owid_vaccine_per_day <- owid_m %>%
    summarise(daily = sum(is.na(vaccinations)),
              cumulative = sum(is.na(vaccinations_cum))) %>%
    transmute(
      iso3c = iso3c,
      use_daily =
        if_else(
          daily <= cumulative,
          TRUE,
          FALSE
        )
    ) %>%
    right_join(owid_m,
               by = "iso3c") %>%
    group_by(iso3c)

  #get doses each day
  owid_vaccine_per_day <- owid_vaccine_per_day %>% #get rid of any tailing zeros depending on method
    filter(!(use_daily & rev(cumsum(rev(!is.na(vaccinations)))) == 0),
           !((!use_daily) & rev(cumsum(rev(!is.na(vaccinations_cum)))) == 0)
    ) %>%
    filter(!(use_daily & cumsum(!is.na(vaccinations)) == 0),
           !((!use_daily) & cumsum(!is.na(vaccinations_cum)) == 0)
    )  %>%
    mutate(
      #set imputed variable
      imputed = if_else(
        (use_daily & is.na(vaccinations)) |
          (!(use_daily) & is.na(vaccinations_cum)),
        TRUE,
        imputed
      ),
      vacc_per_day =
        if_else(
          use_daily,
          #linaerly interpolate the data
          linearly_interpolate(vaccinations),
          diff(c(0, linearly_interpolate(vaccinations_cum)))
        ),
      before_lst = ((1:length(vacc_per_day)) <= which.max(vaccinations_cum)),
      lt_d = use_daily & before_lst,
      prop_multi = ((max(vaccinations_cum, na.rm = TRUE) -
                      sum(vacc_per_day[before_lst]) +
                      sum(vacc_per_day[lt_d])) /
        sum(vacc_per_day[lt_d])),
       vacc_per_day =
         if_else(
          lt_d,
           #scale so that imputed values add to last cumulative value (if before)
           prop_multi *
             vacc_per_day,
           vacc_per_day
         )
    ) %>%
    select(iso3c, date, vacc_per_day, imputed) %>%
    na.omit() #removes some leading nas

  #now add the second/first doses
  owid_merge <-  owid_vaccine_per_day %>%
    left_join(owid_m %>% select(iso3c, date, first_dose_cum, second_dose_cum) %>%
                filter(!is.na(first_dose_cum) | !is.na(second_dose_cum)),
              by = c("iso3c", "date")
    ) %>% #linearly interpolate assuming that both start at 0 and second
    #dose stays at 0 until 18 days
    group_by(iso3c) %>%
    mutate(
      first_dose_cum = if_else(
        1:length(first_dose_cum) == 1 & is.na(first_dose_cum),
        0,
        first_dose_cum
      ),
      second_dose_cum = if_else(
        1:length(second_dose_cum) <= 18,
        0,
        second_dose_cum
      ),
      first_dose_cum = linearly_interpolate(first_dose_cum),
      second_dose_cum = linearly_interpolate(second_dose_cum),
      #adjust cumulative so they are consistent with the vaccinations per day
      percentage_second = if_else(
        first_dose_cum == 0,
        0,
        diff(c(0, second_dose_cum))/
          (diff(c(0, first_dose_cum)) +
             diff(c(0, second_dose_cum)))
      ),
      first_dose_per_day = vacc_per_day*(1-percentage_second),
      second_dose_per_day = vacc_per_day*percentage_second
      )

  #sort out NAs where there are no new doses in either, we'll just assume values
  #that have no effect on dose ratio
  #due to repition needed we have to do it outside dplyr
  missing_df <- owid_merge %>% summarise(
    missing = sum(is.na(first_dose_per_day))
  )

  while(sum(missing_df$missing) > 0){
    owid_merge <- mutate(owid_merge,
      prev_dose_ratio =
        lag(cumsum(second_dose_per_day), 1)/
        lag(cumsum(first_dose_per_day), 1),
      first_dose_per_day = if_else(
        is.na(first_dose_per_day),
        vacc_per_day/(1 + prev_dose_ratio),
        first_dose_per_day
      ),
      second_dose_per_day = if_else(
        is.na(second_dose_per_day),
        vacc_per_day - first_dose_per_day,
        second_dose_per_day
      )
    ) %>%
      select(!prev_dose_ratio)
      #any values after the first will still be NA due to cumsum return NA,
      #so we need to repeat until they are all filled
    missing_df <- owid_merge %>% summarise(
      missing = sum(is.na(first_dose_per_day))
    )

    #message(sum(missing_df$missing))
  }
  owid_merge <- owid_merge %>%
    mutate(
      first_dose_cum =  cumsum(first_dose_per_day),
      second_dose_cum =  cumsum(second_dose_per_day),
      #correct for if second doses are larger than first,
      #assume that these are single dose vaccines and move them over to first dose
      single_doses_vacc = if_else(
        first_dose_cum < second_dose_cum,
        diff(c(0, second_dose_cum)) - diff(c(0, first_dose_cum)),
        0
      ),
      first_dose_per_day =  diff(c(0, first_dose_cum)) + single_doses_vacc,
      second_dose_per_day =  diff(c(0, second_dose_cum)) - single_doses_vacc,


      first_dose_cum = cumsum(first_dose_per_day),
      second_dose_cum = cumsum(second_dose_per_day),
      all_dose_cum = cumsum(vacc_per_day)
    ) %>%
  select(
    iso3c, date, vacc_per_day, first_dose_per_day, second_dose_per_day, imputed
  )

  #now get the dose ratio
  owid_merge <- owid_merge %>%
    mutate(dose_ratio = if_else(
      cumsum(first_dose_per_day) == 0,
      0,
      cumsum(second_dose_per_day)/cumsum(first_dose_per_day)
      )
    )

  #checl that final dose ratio is not 0.1 off from the final dose ratio in owid
  owid_rescaled <- owid_merge %>%
    arrange(iso3c, date) %>%
    summarise(
      dose_ratio = tail(na.omit(dose_ratio),1)
    ) %>%
    left_join(
      owid_m %>%
        group_by(iso3c) %>%
        summarise(
          final_dose_ratio = max(second_dose_cum, na.rm = TRUE)/
            max(first_dose_cum, na.rm = TRUE)
        ),
      by = "iso3c"
    ) %>%
    mutate(
      diff = abs(dose_ratio - final_dose_ratio)
    ) %>%
    filter(diff > 0.15)
  #slight issue but there is not much we can do

  #final adjustments where change in dose ratio is impossible for given number
  #of doses, can occur because we favour the smoothed vaccinations
  problem <- test_for_errors_in_dose_ratio(owid_merge, dose_ratio)
  if(length(problem) > 1){
    stop(paste0("Drop with calculating dose ratios in the following countires: ",
                paste0(problem, collapse = ", ")
    ))
  }
  return(owid_merge)
}
get_dose_ts_who <- function(who_vacc, who_vacc_meta, iso3cs){
  #have to impute a time series, only use if owid not available
  #get any useful info out of meta df
  meta_df <- who_vacc_meta %>%
    filter(iso3c %in% iso3cs) %>%
    group_by(iso3c) %>%
    summarise(
      start_date = min(as.Date(START_DATE), na.rm = TRUE),
      end_date = max(as.Date(END_DATE), na.rm = TRUE),
      aut = min(as.Date(AUTHORIZATION_DATE), na.rm = TRUE)
    ) %>%
    mutate(
      start_date = if_else(
        is.na(start_date),
        aut,
        start_date
      )
    ) %>%
    select(iso3c, start_date, end_date)

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
  #assume dose ratio is 0 for the first 18 days and then build up to the value
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
    if(length(dates) <= 18){
      dr_r_0 <- who_comb$dose_ratio_final[x]/length(dates)
      dose_ratio <- c(
        dr_r_0*seq(1, length(dates))
      )
    } else if (length(dates) <= 18 + 21){
      end_0 <- floor((18/(18+21)) * length(dates))
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
        rep(dr_r_0, 18),
        dr_r_b * seq(1, 21),
        rep(dr_r_pb, length(dates) - 21 - 18)
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
      )
    )
  #final check
  problem <- test_for_errors_in_dose_ratio(ts_df, dose_ratio)
  if(length(problem) > 1){
    stop(paste0("Drop with calculating dose ratios in the following countires: ",
                paste0(problem, collapse = ", ")
                ))
  }
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
