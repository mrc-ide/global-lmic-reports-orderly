# function to remove over estimates of vaccine, i.e. where total vaccinations have gone down
remove_overestimates <- function(tot) {

  tots <- tot[!is.na(tot)]
  if(any(diff(tots) < 0)) {

    off <- which(diff(tots) < 0)
    tot[which(!is.na(tot))[off]] <- NA
    return(tot)
  } else {
    return(tot)
  }

}

# function to interpolate missing vaccine dates
interp_diffs <- function(date_vacc, tot) {

  tot <- remove_overestimates(tot)

  # if there are breaks then we need to interpolate what has happened for ease
  tot <- approx(
    x = as.Date(date_vacc),
    y = tot,
    xout = seq.Date(min(as.Date(date_vacc)), max(as.Date(date_vacc)), 1)
  )$y

  # filter to just dates now with vaccines given
  date_vaccine_change <- seq.Date(min(as.Date(date_vacc)), max(as.Date(date_vacc)), 1)[!is.na(tot)]

  # we should also linearly extend back to 0 vaccinations
  early_vaccs <- predict(
    lm(y~x, data = data.frame(y = c(head(na.omit(tot), min(length(tot), 7))), x = seq_len(min(length(na.omit(tot)), 7)))),
    newdata = data.frame(x = -30:0),
    type = "response")
  early_vaccs <- round(early_vaccs[early_vaccs > 0])

  if(length(early_vaccs) > 0) {

    # corresponding dates
    back <- as.numeric(names(early_vaccs)) - max(as.numeric(names(early_vaccs)))
    extra_dates <- as.Date(date_vaccine_change[1]) - 1 + back

    # add the early vaccinations
    max_vaccine <- c(vapply(early_vaccs, min, numeric(1), na.omit(tot)[1]), na.omit(tot))
    date_vaccine_change <- c(extra_dates, date_vaccine_change)

  } else {

    max_vaccine <- na.omit(tot)

  }

  # create our vacc inputs
  max_vaccine <- as.integer(c(max_vaccine[1], diff(max_vaccine)))
  max_vaccine[max_vaccine < 0] <- 0

  return(list(date = date_vaccine_change, max = max_vaccine))

}

interp_for_dates <- function(date_vacc, tot, dates) {

  tot <- remove_overestimates(tot)

  # if there are breaks then we need to interpolate what has happened for ease
  tot <- approx(
    x = as.Date(date_vacc),
    y = tot,
    xout = seq.Date(min(as.Date(date_vacc)), max(as.Date(date_vacc)), 1)
  )$y

  # filter to just dates now with vaccines given
  date_vaccine_change <- seq.Date(min(as.Date(date_vacc)), max(as.Date(date_vacc)), 1)[!is.na(tot)]
  tot <- tot[!is.na(tot)]

  # we should also linearly extend back to 0 vaccinations
  if(length(dates[dates < date_vaccine_change[1]]) > 0) {

    early_vaccs <- predict(
      lm(y~x, data = data.frame(y = head(tot, min(7, length(tot))), x = seq_len(min(7, length(tot))))),
      newdata = data.frame(x = as.integer(-((as.Date(date_vaccine_change[1]) - as.Date(dates[dates < date_vaccine_change[1]]))-1))),
      type = "response")
    early_vaccs[early_vaccs < 0] <- 0
    early_vaccs <- round(early_vaccs[early_vaccs >= 0])

  } else {
    early_vaccs <- numeric(0)
  }

  if(length(early_vaccs) > 0) {

    # corresponding dates
    back <- as.numeric(names(early_vaccs)) - max(as.numeric(names(early_vaccs)))
    extra_dates <- as.Date(date_vaccine_change[1]) - 1 + back

    # add the early vaccinations
    max_vaccine <- c(early_vaccs, tot)
    date_vaccine_change <- c(extra_dates, date_vaccine_change)

  } else {

    max_vaccine <- tot

  }

  # now we need to check for linear interpolation forward in time as well
  if(any(!dates %in% date_vaccine_change)) {

    y_in <- tail(tot, min(7, length(max_vaccine)))
    x_in <- seq_len(min(7, length(max_vaccine)))
    x_in <- x_in - max(x_in)


    late_vaccs <- predict(
      lm(y~x, data = data.frame(y = y_in, x = x_in)),
      newdata = data.frame(x = seq_along(dates[!dates %in% date_vaccine_change])),
      type = "response")
    late_vaccs[late_vaccs < 0] <- 0
    late_vaccs <- round(late_vaccs[late_vaccs >= 0])

    max_vaccine <- c(max_vaccine, late_vaccs)
    date_vaccine_change <- dates

  }

  # create our vacc inputs
  max_vaccine <- as.integer(c(max_vaccine[1], diff(max_vaccine)))
  max_vaccine[max_vaccine < 0] <- 0

  return(list(date = date_vaccine_change, max = max_vaccine))

}
