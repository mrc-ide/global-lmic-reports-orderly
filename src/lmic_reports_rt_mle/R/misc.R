r_list_format <- function(out, date_0) {

  df <- nimue_format(out,
                      var_select = c("infections","deaths","hospital_demand",
                                     "ICU_demand", "D", "hospital_incidence","ICU_incidence"),
                      date_0 = date_0)

  pr <- nimue_format(out, var_select = c("S","R","D"), date_0 = date_0) %>%
    na.omit %>%
    pivot_wider(names_from = compartment, values_from = y) %>%
    mutate(y = sum(out$parameters$population)-D-R-S,
           compartment = "prevalence") %>%
    select(replicate, compartment, t, y, date)

  rbind(df, pr)

}
named_list <- function(...) {
  get <- as.character(match.call())
  l <- list(...)
  names(l) <- utils::tail(get, -1)
  return(l)
}

summarise_rt_futures <- function(Rt_futures){
  #for each trend type
  summarised <- lapply(names(Rt_futures), function(trend){
    #for each replicate
    do.call(
      rbind,
      lapply(seq_along(Rt_futures[[trend]]), function(rep){
        #just get the mean
        mean(Rt_futures[[trend]][[rep]])
      })
    )
  })
  #summarise and put into dataframe
  out <- data.frame(temp = 1
  )
  for(i in seq_along(Rt_futures)){
    out <- out %>%
      mutate(
        !! paste0(names(Rt_futures)[i], "_med") := median(summarised[[i]]),
        !! paste0(names(Rt_futures)[i], "_025") := quantile(summarised[[i]], 0.025),
        !! paste0(names(Rt_futures)[i], "_975") := quantile(summarised[[i]], 0.975)
      )
  }
  return(out %>%
           select(!temp))
}


variant_changes_over_time <- function(variant_characteristics,
                                      var, start_date) {

  #delta
  dates <- seq(
    variant_characteristics$Delta$start_date,
    variant_characteristics$Delta$start_date +
      variant_characteristics$Delta$shift_duration,
    by = 1
  )
  variable <- c(seq(variant_characteristics$Wild[[var]],
                    variant_characteristics$Delta[[var]],
                    length.out = length(dates)+1))
  #omicron
  dates_omicron <- seq(
    variant_characteristics$Omicron$start_date,
    variant_characteristics$Omicron$start_date +
      variant_characteristics$Omicron$shift_duration,
    by = 1
  )
  variable <- c(
    variable,
    seq(variant_characteristics$Delta[[var]],
        variant_characteristics$Omicron[[var]],
        length.out = length(dates_omicron)+1)[-1]
  )
  dates <- c(dates, dates_omicron)
  tt <- as.numeric(dates - start_date)
  if(tt[1] > 0){
    tt <- c(0, tt)
  } else {
    tt <- c(tt[1] - 1, tt)
  }
  return(
    list(
      var = variable,
      tt = tt
    )
  )
}

variant_immune_escape <- function(variant_characteristics, dur_R, start_date) {

  #delta
  dur_R_d <- c(dur_R, 1 / (
    (variant_characteristics$Delta$shift_duration / dur_R - log(1 - variant_characteristics$Delta$immune_escape)) /
      variant_characteristics$Delta$shift_duration
  ), dur_R)
  date_dur_R_change_d <- c(variant_characteristics$Delta$start_date,
                           variant_characteristics$Delta$start_date +
                             variant_characteristics$Delta$shift_duration)
  #omicron
  dur_R_o <- c( 1 / (
    (variant_characteristics$Omicron$shift_duration / dur_R - log(1 - variant_characteristics$Omicron$immune_escape)) /
      variant_characteristics$Omicron$shift_duration
  ), dur_R)
  date_dur_R_change_o <- c(variant_characteristics$Omicron$start_date,
                           variant_characteristics$Omicron$start_date +
                             variant_characteristics$Omicron$shift_duration)

  tt <- as.numeric(c(date_dur_R_change_d, date_dur_R_change_o) - start_date)
  if(tt[1] > 0){
    tt <- c(0, tt)
  } else {
    tt <- c(tt[1] - 1, tt)
  }
  return(
    list(
      variable = c(dur_R_d, dur_R_o),
      tt = tt
    )
  )
}
