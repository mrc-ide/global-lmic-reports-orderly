init_state_nimue <- function(deaths_removed, iso3c, seeding_cases = 5,
                             vaccinated_already = 0) {

  # get an initial
  pop <- squire::get_population(iso3c = iso3c, simple_SEIR = FALSE)
  init <- nimue:::init(pop$n, seeding_cases = seeding_cases)

  if(deaths_removed > 0) {
    # work out how many deaths and where
    probs <- (squire:::probs$prob_hosp * squire:::probs$prob_severe * squire:::probs$prob_severe_death_treatment) +
      (squire:::probs$prob_hosp * (1-squire:::probs$prob_severe * squire:::probs$prob_non_severe_death_treatment))
    probs <- probs*pop$n
    probs <- probs/sum(probs)
    deaths <- as.numeric(t(rmultinom(1, deaths_removed, probs)))
    # approximate IFR for income group
    wb_metadata <- read.csv("gdp_income_group.csv", fileEncoding="UTF-8-BOM", stringsAsFactors = TRUE)
    income <- wb_metadata$income_group[match(iso3c, wb_metadata$country_code)]
    ifrs <- data.frame("income" = c("Low income", "Lower middle income", "Upper middle income", "High income"),
                       "ifr" = c(0.17, 0.31, 0.51, 1.02))
    ifr <- ifrs$ifr[ifrs$income == income]
    R <- rpois(1, deaths_removed*1/ifr/0.01)
    R <- as.numeric(t(rmultinom(1, R, rep(1/length(probs), length(probs)))))
    R <- R - deaths
    # and update the inital to reflect
    init$D_0[,1] <- deaths
    init$S_0[,1] <- init$S_0[,1] - R - deaths
    init$R1_0[,1] <- R
  }

  if(vaccinated_already > 0) {

    S_0 <- init$S_0
    prop <- vaccinated_already / sum(tail(S_0[,1],-3))
    S_0[-(1:3),3] <- round(S_0[-(1:3),1] * prop)
    S_0[-(1:3),1] <- S_0[-(1:3),1] - S_0[-(1:3),3]
    init$S_0 <- S_0
  }

  #drop last column
  init <- purrr::map(init, ~.x[,-6])

  return(init)
}

post_lockdown_date_relative <- function(x, above = 1.1, max_date, min_date) {

  if(nrow(x)==0) {

    return(NA)

  } else {

    if(any(x$observed)) {
      m <- predict(loess(C~as.numeric(date), data=x, span = 0.2), type = "response")
    } else {
      m <- predict(loess(C~as.numeric(date), data=x, span = 0.2), type = "response")
      #m <- x$C
    }
    min_mob <- min(m)
    diff <- 1 - min_mob

    pl <- NA
    attempts <- 50
    while(is.na(pl) && attempts == 0) {
      above15 <- which(m >= ((above-1)*diff)+min_mob)
      pl <- above15[which(above15>which.min(m))[1]]
      above <- above*0.99
      attempts <- attempts -1
    }

    # if still NA then just return max date
    if(is.na(pl)) {
      return(max_date)
    }

    # if past max date then take the last local  minimum and grow by 4 days
    if(x$date[pl] > max_date) {
      min_f <- which(diff(sign(diff(m)))==2)+1
      pl <- min_f[tail(which(x$date[min_f] < max_date),1)] + 4
    }

    dat <- max(min_date, as.Date(x$date[pl]))

    return(dat)

  }

}


