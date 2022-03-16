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
data <- readRDS("combined_data.Rds") %>%
  mutate(date = dateRep) %>%
  filter(countryterritoryCode == iso3c) %>%
  rename(iso3c = countryterritoryCode)

#this step removes deaths that were likely due to importations that lead nowhere
#or small contained epidemics simply to help the model fit better
data <- preprocess_for_fitting(data, iso3c, date_0)
removed_deaths <- data[[2]] #save to store later
data <- data[[1]]
#drop iso3c code
data$iso3c <- NULL

## b. Sort out what is to be our death time series
## -----------------------------------------------------------------------------

# here I have just taken the maximum of either excess or covid on each day.
# but this is probably not the best idea as it likely overestimates covid deaths

# I'm not sure this does over estimate, if excess mortality is disconnected to
# covid deaths then we are essentially using reported deaths (an under-estimate)
# Only case I can think of for over-estimate is when excess deaths spikes due to
# lack of treatment etc, whilst covid deaths are well reported and tested for.

#if we have no deaths then we do not proceed

if(nrow(data) == 0 | sum(data$deaths) == 0 | is.na(sum(data$deaths))){
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
  delta_characteristics <- readRDS("variant_characteristics.Rds")[[iso3c]]$Delta

  vaccine_inputs <- readRDS("vacc_inputs.Rds")[[iso3c]]

  likelihood_version <- "Negative Binomial-Cumulative"

  ## -----------------------------------------------------------------------------
  ## 2. Fit Model
  ## -----------------------------------------------------------------------------

  # fit model
  res <- fit_spline_rt(
    data = data,
    country = country,
    pop = pop,
    n_mcmc = as.numeric(n_mcmc),
    n_chains = as.numeric(n_chains),
    n_burnin = as.numeric(n_burnin),
    replicates = as.numeric(replicates),
    delta_characteristics = delta_characteristics,
    vaccine_inputs = vaccine_inputs,
    likelihood_version = likelihood_version
  )

  ## -----------------------------------------------------------------------------
  ## 3. Summarise model for ease of viewing outputs and goodness of fit
  ## -----------------------------------------------------------------------------

  # remove the output for memory and ease
  output <- res$output
  res$output <- NULL

  #also remove drjacoby results, too large and we shouldn't need them
  drjacoby_out <- res$pmcmc_results$drjacoby_out
  res$pmcmc_results$drjacoby_out <- NULL

  #add removed deaths to interventions list
  res$interventions$pre_epidemic_isolated_deaths <- removed_deaths

  # save output without output for memory
  saveRDS(res, "res.rds")
  res$output <- output

  # make a series of quick plots so we can check fits easily afterwards
  #rtp <- rt_plot_immunity(res, vaccine = TRUE, Rt_plot = TRUE)
  dp <- dp_plot(res)
  cdp <- cdp_plot(res)
  #ar <- ar_plot(res)
  coupling <- drjacoby::plot_mc_acceptance(drjacoby_out)
  #multiplot
  parameters <- setdiff(names(drjacoby_out$output),
                        c("chain", "phase", "iteration", "logprior", "loglikelihood")
  )

  log_plots <- ggplot(data = drjacoby_out$output %>%
           pivot_longer(cols = all_of(parameters),
                        names_to = "parameter",
                        values_to = "value") %>%
             filter(loglikelihood != -.Machine$double.xmax), aes(x = value, y = logprior + loglikelihood,
                                                  colour = as.character(chain))) +
    geom_point(alpha = 0.25) +
    facet_wrap(vars(parameter), scales = "free_x") +
    theme_pubclean() +
    theme(legend.position = "none") +
    labs(x = "", y = "Log-Posterior") +
    scale_y_continuous(labels = scales::scientific)

  ggsave("fitting.pdf",width=18, height=12,
         cowplot::plot_grid(
           log_plots,
           cowplot::plot_grid(dp + ggtitle(country),
                            cdp, coupling, ncol = 1),
           nrow = 1)
  )

}

