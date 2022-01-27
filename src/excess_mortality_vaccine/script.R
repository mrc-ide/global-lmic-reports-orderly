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

#this step removes deaths that were likely due to importations that lead nowhere
#or small contained epidemics simply to help the model fit better
df <- preprocess_for_fitting(df)
removed_deaths <- df[[2]] #save to store later
df <- df[[1]]

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


    #open data from covariants, we won't use omicron adjustments
    delta_characteristics <- readRDS("variant_characteristics.Rds") %>%
      ungroup() %>%
      rename(iso3c_ = iso3c) %>%
      filter(iso3c_ == iso3c) %>%
      select(where(~is.numeric(.x) | is.Date(.x))) %>%
      select(!contains("omicron"))

    vaccine_inputs <- readRDS("vacc_inputs.Rds")[[iso3c]]


  #use the poisson distribution if we need to close in more quickly,
  #should not be used for final fits
  if(iso3c %in% c(
    # "COL", "TZA", "TCD", "UGA", "CAN", "SWE", "AFG", "KEN", "ZAF",
    # "GIN", "MOZ", "NER", "SAU",
    # "CUB", "NGA", "SDN"
                  )){
    version <- "Poisson"
  } else {
    version <- "Negative Binomial-Cumulative"
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
    n_chains = as.numeric(n_chains),
    replicates = as.numeric(replicates),
    delta_characteristics = delta_characteristics,
    vaccine_inputs = vaccine_inputs,
    likelihood_version = version
  )

  ## -----------------------------------------------------------------------------
  ## 3. Summarise model for ease of viewing outputs and goodness of fit
  ## -----------------------------------------------------------------------------

  # remove the output for memory and ease
  output <- res$output
  res$output <- NULL

  #add removed deaths to interventions list
  res$interventions$pre_epidemic_isolated_deaths <- removed_deaths

  # save output without output for memory
  saveRDS(res, "res.rds")
  res$output <- output

  # make a series of quick plots so we can check fits easily afterwards
  rtp <- rt_plot_immunity(res, vaccine = TRUE, Rt_plot = TRUE)
  dp <- dp_plot(res)
  cdp <- cdp_plot(res)
  ar <- ar_plot(res)

  ggsave("fitting.pdf",width=12, height=12,
         cowplot::plot_grid(rtp$plot + ggtitle(country),
                            dp, cdp, ar, ncol = 1))

}

plot(res$pmcmc_results$results$ves, type = "l")
plot(res$pmcmc_results$results$delta_dur_R, type = "l")
cov2cor(res$pmcmc_results$covariance_matrix[[1]])

plot(res$pmcmc_results$chains$chain2$results$ves, type = "l")
plot(res$pmcmc_results$chains$chain2$results$delta_dur_R, type = "l")
res$pmcmc_results$chains$chain1$covariance_matrix

