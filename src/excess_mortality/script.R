RhpcBLASctl::blas_set_num_threads(1L)
RhpcBLASctl::omp_set_num_threads(1L)

system(paste0("echo Excess Fit for  ",iso3c))
date_0 <- as.Date(date)

version_min <- "0.6.9"
if(packageVersion("squire") < version_min) {
  stop("squire needs to be updated to at least v", version_min)
}

version_min <- "0.1.22"
if(packageVersion("nimue") < version_min) {
  stop("nimue needs to be updated to at least v", version_min)
}

## -----------------------------------------------------------------------------
## 1. GET INPUT DATA
## -----------------------------------------------------------------------------

## a. Get data from local files
## -----------------------------------------------------------------------------

# get data from file
data <- readRDS("excess_deaths.Rds")
df <- data[data$iso3c == iso3c, ]

#add a check for brunei, a long gap of no deaths makes this is impossible to fit
if(iso3c == "BRN"){
  #proceed if there's a lot of space with less than 1 tenth of the deaths
  if(sum(df[df$week_start < "2021-06-20",]$deaths) < sum(df$deaths)/10){
    df <- df[df$week_start >= "2021-06-20",]
  }
}

## b. Sort out what is to be our death time series
## -----------------------------------------------------------------------------

# here I have just taken the maximum of either excess or covid on each day.
# but this is probably not the best idea as it likely overestimates covid deaths

# I'm not sure this does over estimate, if excess mortality is disconnected to
# covid deaths then we are essentially using reported deaths (an under-estimate)
# Only case I can think of for over-estimate is when excess deaths spikes due to
# lack of treatment etc, whilst covid deaths are well reported and tested for.

#if we have no deaths then we do not proceed

if(nrow(df) == 0 | sum(df$deaths) == 0){
  saveRDS(NULL, "res.rds")
  ggsave("fitting.pdf",width=12, height=12,
         NULL)
} else{
  ## c. Any other parameters needed to be worked out
  ## -----------------------------------------------------------------------------

  # check that we have this iso3c in squire
  if(!(iso3c %in% squire::population$iso3c)) {
    stop("iso3c not found in squire")
  }
  country <- squire::population$country[match(iso3c, squire::population$iso3c)]
  pop <- squire::get_population(country)$n

  if(model == "SQUIRE"){
    adjust_delta <- FALSE
  }

  if(adjust_delta){
    #open data from covariants
    delta_characteristics <- readRDS("delta_characteristics.Rds") %>% ungroup()

    #get data or impute if not there
    if (iso3c %in% delta_characteristics$iso3c) {
      this_iso3c <- iso3c
      delta_characteristics <- delta_characteristics %>%
        filter(iso3c == this_iso3c) %>%
        select(where(~is.numeric(.x) | is.Date(.x)))
    } else if (
      countrycode::countrycode(iso3c, origin = "iso3c",
                               destination = "un.regionsub.name") %in%
      delta_characteristics$sub_region
    ) {
      this_sub_region <- countrycode::countrycode(iso3c,
                                                  origin = "iso3c", destination = "un.regionsub.name")
      #we then use the median values of all countries in that sub region
      delta_characteristics <- delta_characteristics %>%
        filter(
          sub_region == this_sub_region
        ) %>%
        summarise(across(
          where(~is.numeric(.x) | is.Date(.x)),
          ~median(.x, na.rm = T)
        ))
    } else if (
      countrycode::countrycode(iso3c, origin = "iso3c",
                               destination = "continent") %in%
      delta_characteristics$continent
    ) {
      this_continent <- countrycode::countrycode(iso3c, origin = "iso3c",
                                                 destination = "continent")
      #we then use the median values of all countries in that contient
      delta_characteristics <- delta_characteristics %>%
        filter(
          continent == this_continent
        ) %>%
        summarise(across(
          where(~is.numeric(.x) | is.Date(.x)),
          ~median(.x, na.rm = T)
        ))
    } else{
      #else we use the median values for the world
      delta_characteristics <- delta_characteristics %>%
        summarise(across(
          where(~is.numeric(.x) | is.Date(.x)),
          ~median(.x, na.rm = T)
        ))
    }
  } else{
    #these settings should lead to no adjustment
    delta_characteristics <- data.frame(
    )
  }

  ## -----------------------------------------------------------------------------
  ## 2. Fit Model
  ## -----------------------------------------------------------------------------

  # fit model
  res <- fit_spline_rt(
    data = df,
    country = country,
    pop = pop,
    n_mcmc = as.numeric(n_mcmc),
    replicates = as.numeric(replicates),
    model = model,
    delta_characteristics = delta_characteristics
  )

  ## -----------------------------------------------------------------------------
  ## 3. Summarise model for ease of viewing outputs and goodness of fit
  ## -----------------------------------------------------------------------------

  # remove the output for memory and ease
  output <- res$output
  res$output <- NULL

  # save output without output for memory
  saveRDS(res, "res.rds")
  res$output <- output

  # make a series of quick plots so we can check fits easily afterwards
  if (model == "SQUIRE") {
    rtp <- rt_plot_immunity(res)
  } else {
    rtp <- rt_plot_immunity_vaccine(res)
  }

  #just for this
  res$pmcmc_results$inputs$data$date <- res$pmcmc_results$inputs$data$week_start

  dp <- dp_plot(res)
  cdp <- cdp_plot(res)
  ar <- ar_plot(res)

  res$pmcmc_results$inputs$data$date <- NULL

  ggsave("fitting.pdf",width=12, height=12,
         cowplot::plot_grid(rtp$plot + ggtitle(country),
                            dp, cdp, ar, ncol = 1))

}
