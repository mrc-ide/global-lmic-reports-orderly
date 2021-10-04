#define a function to calculate vaccine efficientcy for a given set of values
get_ret_res <- function(owid, who_vacc, who_vacc_meta, date_vaccine_change,
                        ve_i_low, ve_i_high, ve_d_low, ve_d_high){
  # does owid have no data
  if(all(is.na(owid$total_vaccinations))) {

    # next then check for who data

    # if this has no data then return 0 for simplicity
    if(nrow(who_vacc) == 0) {

      #must assume dose ratio for future use in counter factuals

      dose_ratio <- 0
      second_dose <- 0

      vaccine_efficacy_infection <- rep((1-dose_ratio)*ve_i_low + dose_ratio*ve_i_high, 1)
      vaccine_efficacy_disease <- rep((1-dose_ratio)*ve_d_low + dose_ratio*ve_d_high, 1)
      vaccine_efficacy_infection <- lapply(vaccine_efficacy_infection, rep, 17)
      vaccine_efficacy_disease <- lapply(vaccine_efficacy_disease, rep, 17)

      ret_res <- list(
        date_vaccine_change = date_0,
        max_vaccine = rep(0, 1),
        vaccine_efficacy_infection = vaccine_efficacy_infection,
        vaccine_efficacy_disease = vaccine_efficacy_disease
      )

    } else if (is.na(who_vacc$TOTAL_VACCINATIONS)) {

      #must assume dose ratio for future use in counter factuals

      dose_ratio <- 0
      second_dose <- 0

      vaccine_efficacy_infection <- rep((1-dose_ratio)*ve_i_low + dose_ratio*ve_i_high, 1)
      vaccine_efficacy_disease <- rep((1-dose_ratio)*ve_d_low + dose_ratio*ve_d_high, 1)
      vaccine_efficacy_infection <- lapply(vaccine_efficacy_infection, rep, 17)
      vaccine_efficacy_disease <- lapply(vaccine_efficacy_disease, rep, 17)

      ret_res <- list(
        date_vaccine_change = date_0,
        max_vaccine = rep(0, 1),
        vaccine_efficacy_infection = vaccine_efficacy_infection,
        vaccine_efficacy_disease = vaccine_efficacy_disease
      )

    } else {

      # first see if they have persons vaccinated or just total.
      if(is.na(who_vacc$PERSONS_VACCINATED_1PLUS_DOSE)) {
        pers_vacc <- who_vacc$TOTAL_VACCINATIONS
      } else {
        pers_vacc <- who_vacc$PERSONS_VACCINATED_1PLUS_DOSE
      }

      # if there is no start date just distribute the vaccines equally over the last 1 month for now
      if(nrow(who_vacc_meta) == 0) {

        date_vaccine_change <- seq.Date(as.Date(who_vacc$DATE_UPDATED) - 30, as.Date(who_vacc$DATE_UPDATED), 1)

      } else if (nrow(who_vacc_meta) >= 1 && all(is.na(who_vacc_meta$START_DATE))) {

        date_vaccine_change <- seq.Date(as.Date(who_vacc$DATE_UPDATED) -30, as.Date(who_vacc$DATE_UPDATED), 1)

      } else {

        date_vaccine_change <- seq.Date(min(as.Date(who_vacc_meta$START_DATE), na.rm=TRUE), as.Date(who_vacc$DATE_UPDATED), 1)

      }

      # assume took 1 month if still only one data point
      if(length(unique(date_vaccine_change)) == 1) {
        date_vaccine_change <- seq.Date(as.Date(date_vaccine_change[1])-30, as.Date(date_vaccine_change[1]), 1)
      }

      max_vaccine <- round(pers_vacc/length(date_vaccine_change))
      max_vaccine <- rep(max_vaccine, length(date_vaccine_change))


      # ratio of 1st to 2nd doses given
      if(is.na(who_vacc$PERSONS_VACCINATED_1PLUS_DOSE)) {
        dose_ratio <- 0.5
      } else {
        second_dose <- who_vacc$TOTAL_VACCINATIONS - who_vacc$PERSONS_VACCINATED_1PLUS_DOSE
        dose_ratio <- second_dose/who_vacc$PERSONS_VACCINATED_1PLUS_DOSE
        if(is.nan(dose_ratio)) {
          dose_ratio <- 0
        }
      }

      # format for odin
      vaccine_efficacy_infection <- (1-dose_ratio)*ve_i_low + dose_ratio*ve_i_high
      vaccine_efficacy_disease <- (1-dose_ratio)*ve_d_low + dose_ratio*ve_d_high
      vaccine_efficacy_infection <- c(rep(ve_i_low,28),  rep(vaccine_efficacy_infection, length.out = max(length(max_vaccine) - 28, 1)))
      vaccine_efficacy_disease <- c(rep(ve_d_low,28),  rep(vaccine_efficacy_disease, length.out = max(length(max_vaccine) - 28, 1)))
      vaccine_efficacy_infection <- lapply(vaccine_efficacy_infection, rep, 17)
      vaccine_efficacy_disease <- lapply(vaccine_efficacy_disease, rep, 17)

      ret_res <- list(
        date_vaccine_change = date_vaccine_change,
        max_vaccine = max_vaccine,
        vaccine_efficacy_infection = vaccine_efficacy_infection,
        vaccine_efficacy_disease = vaccine_efficacy_disease
      )
    }

  } else if (sum(!is.na(owid$people_vaccinated)) == 1) {

    # first is there any additional data in the WHO
    if(nrow(who_vacc) == 1) {

      if(!is.na(who_vacc$PERSONS_VACCINATED_1PLUS_DOSE)) {

        # if we have an extra data point here then add it in
        date_vaccine_change <- c(who_vacc$DATE_UPDATED, owid$date[!is.na(owid$people_vaccinated)])
        max_vaccine <- c(who_vacc$PERSONS_VACCINATED_1PLUS_DOSE, owid$people_vaccinated[!is.na(owid$people_vaccinated)])

        # and is there a start date
        if(nrow(who_vacc_meta) >= 1) {
          if(any(!is.na(who_vacc_meta$START_DATE))) {
            date_vaccine_change <- c(min(as.Date(who_vacc_meta$START_DATE), na.rm = TRUE), date_vaccine_change)
            max_vaccine <- c(0, max_vaccine)
          }
        }

        # order them
        max_vaccine <- max_vaccine[order(date_vaccine_change)]
        date_vaccine_change <- date_vaccine_change[order(date_vaccine_change)]

      } else {

        date_vaccine_change <- owid$date[!is.na(owid$people_vaccinated)]
        max_vaccine <- owid$people_vaccinated[!is.na(owid$people_vaccinated)]

        # and is there a start date
        if(nrow(who_vacc_meta) >= 1) {
          if(any(!is.na(who_vacc_meta$START_DATE))) {
            date_vaccine_change <- c(min(as.Date(who_vacc_meta$START_DATE), na.rm = TRUE), date_vaccine_change)
            max_vaccine <- c(0, max_vaccine)
          }
        }
      }

    } else {

      date_vaccine_change <- owid$date[!is.na(owid$people_vaccinated)]
      max_vaccine <- owid$people_vaccinated[!is.na(owid$people_vaccinated)]

      # and is there a start date
      if(nrow(who_vacc_meta) >= 1) {
        if(any(!is.na(who_vacc_meta$START_DATE))) {
          date_vaccine_change <- c(min(as.Date(who_vacc_meta$START_DATE), na.rm = TRUE), date_vaccine_change)
          max_vaccine <- c(0, max_vaccine)
        }
      }

    }

    # assume took 1 month if still only one data point
    if(length(unique(date_vaccine_change)) == 1) {
      max_vaccine <- c(0, max(max_vaccine))
      date_vaccine_change <- c(as.Date(date_vaccine_change[1])-30, as.Date(date_vaccine_change[1]))
    }

    peeps <- interp_diffs(date_vacc = date_vaccine_change, tot = max_vaccine)
    date_vaccine_change <- peeps$date
    max_vaccine <- peeps$max

    # ratio of 1st to 2nd doses given
    if (all(is.na(owid$people_vaccinated))) {
      if(nrow(who_vacc) == 1) {
        if(is.na(who_vacc$PERSONS_VACCINATED_1PLUS_DOSE)) {
          dose_ratio <- 0.5
        } else {
          second <- who_vacc$TOTAL_VACCINATIONS - who_vacc$PERSONS_VACCINATED_1PLUS_DOSE
          dose_ratio <- second/who_vacc$PERSONS_VACCINATED_1PLUS_DOSE
          dose_ratio[seq_len(min(length(dose_ratio), 28))] <- 0
        }
      } else {
        dose_ratio <- 0.5
      }
    } else {
      dose_ratio <- mean(1 - (owid$people_vaccinated[!is.na(owid$people_vaccinated)]/owid$total_vaccinations[!is.na(owid$people_vaccinated)]))
    }

    # format for odin
    vaccine_efficacy_infection <- (1-dose_ratio)*ve_i_low + dose_ratio*ve_i_high
    vaccine_efficacy_disease <- (1-dose_ratio)*ve_d_low + dose_ratio*ve_d_high
    vaccine_efficacy_infection <- seq(ve_i_low, vaccine_efficacy_infection, length.out = length(max_vaccine))
    vaccine_efficacy_disease <- seq(ve_d_low, vaccine_efficacy_disease, length.out = length(max_vaccine))
    vaccine_efficacy_infection <- lapply(vaccine_efficacy_infection, rep, 17)
    vaccine_efficacy_disease <- lapply(vaccine_efficacy_disease, rep, 17)

    ret_res <- list(
      date_vaccine_change = date_vaccine_change,
      max_vaccine = max_vaccine,
      vaccine_efficacy_infection = vaccine_efficacy_infection,
      vaccine_efficacy_disease = vaccine_efficacy_disease
    )



  } else {

    # if no indication on split in dosng assume all 1st dose
    if(all(is.na(owid$people_vaccinated))) {

      if(sum(!is.na(owid$total_vaccinations)) == 1) {

        if(nrow(who_vacc) == 1) {

          if(!is.na(who_vacc$TOTAL_VACCINATIONS)) {

            # if we have an extra data point here then add it in
            date_vaccine_change <- c(who_vacc$DATE_UPDATED, owid$date[!is.na(owid$total_vaccinations)])
            max_vaccine <- c(who_vacc$TOTAL_VACCINATIONS, owid$total_vaccinations[!is.na(owid$total_vaccinations)])

            # and is there a start date
            if(nrow(who_vacc_meta) >= 1) {
              if(any(!is.na(who_vacc_meta$START_DATE))) {
                if(any(who_vacc_meta$START_DATE != "")) {
                  date_vaccine_change <- c(min(as.Date(who_vacc_meta$START_DATE), na.rm = TRUE), date_vaccine_change)
                  max_vaccine <- c(0, max_vaccine)
                }
              }
            }

            # order them
            max_vaccine <- max_vaccine[order(date_vaccine_change)]
            date_vaccine_change <- date_vaccine_change[order(date_vaccine_change)]

            # assume took 1 month.
            if(length(unique(date_vaccine_change)) == 1) {
              max_vaccine <- c(0, max(max_vaccine))
              date_vaccine_change <- c(as.Date(date_vaccine_change[1])-30, as.Date(date_vaccine_change[1]))
            }

          }

        }

        peeps <- interp_diffs(date_vacc = date_vaccine_change, tot = max_vaccine)

      } else {

        peeps <- interp_diffs(date_vacc = owid$date, tot = owid$total_vaccinations)

      }
      # peeps vaccinations given out per day
      date_vaccine_change <- peeps$date
      max_vaccine <- peeps$max
      dose_ratio <- rep(0, length(date_vaccine_change))

    } else {

      # interpolate the total vaccinations. Some though have most recnt data only for
      # people vaccinated in which case do that first and use the dates from that for the
      # total vaccination interpolation
      min_date <- min(c(owid$date[!is.na(owid$people_vaccinated)], owid$date[!is.na(owid$total_vaccinations)]))
      max_date <- max(c(owid$date[!is.na(owid$people_vaccinated)], owid$date[!is.na(owid$total_vaccinations)]))
      dates_to_interp <- seq.Date(as.Date(min_date), as.Date(max_date), 1)
      tots <- interp_for_dates(date_vacc = owid$date, tot = owid$total_vaccinations, dates = dates_to_interp)
      peeps <- interp_for_dates(date_vacc = owid$date, tot = owid$people_vaccinated, dates = dates_to_interp)

      # peeps vaccinations given out per day
      date_vaccine_change <- peeps$date
      max_vaccine <- peeps$max

      # now for doses
      firsts <- tots$max
      if (all(is.na(owid$people_fully_vaccinated))) {
        dose_ratio <- rep(0, length(date_vaccine_change))
      } else {

        seconds <- tots$max - peeps$max
        seconds[seconds < 0] <- 0
        firsts <- cumsum(firsts)
        seconds <- cumsum(seconds)
        dose_ratio <- vapply(seconds/firsts, min, numeric(1), 1.0)
        dose_ratio[is.na(dose_ratio)] <- 0
        dose_ratio <- dose_ratio[which(tots$date %in% peeps$date)]

        # and now let sort it out to be more realistic
        # ratio should start at 0 for the first 28 days de facto
        dose_ratio[seq_len(min(length(dose_ratio), 28))] <- 0

      }

    }

    # possibly best to interpolate

    # now to work out the efficacy
    vaccine_efficacy_infection <- (1-dose_ratio)*ve_i_low + dose_ratio*ve_i_high
    vaccine_efficacy_disease <- (1-dose_ratio)*ve_d_low + dose_ratio*ve_d_high
    vaccine_efficacy_infection <- lapply(vaccine_efficacy_infection, rep, 17)
    vaccine_efficacy_disease <- lapply(vaccine_efficacy_disease, rep, 17)

    ret_res <- list(
      date_vaccine_change = date_vaccine_change,
      max_vaccine = max_vaccine,
      vaccine_efficacy_infection = vaccine_efficacy_infection,
      vaccine_efficacy_disease = vaccine_efficacy_disease
    )
  }
  ret_res$dose_ratio <- dose_ratio #append this for future use
  return(ret_res)
}
