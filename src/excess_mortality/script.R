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
data <- read.csv("export_country.csv")
data2 <- data[data$iso3c == iso3c, ] %>%
  mutate(date = as.Date(date)) %>%
  arrange(date) %>%
  select(date, estimated_daily_excess_deaths, daily_covid_deaths) %>%
  rename(deaths = estimated_daily_excess_deaths,
         covid = daily_covid_deaths) %>%
  complete(date = seq.Date(min(date), max(date), 1)) %>%
  fill(deaths, covid,  .direction = "down") %>%
  mutate(only_pos = deaths) %>%
  mutate(only_pos = replace(only_pos, which(only_pos < 0), 0)) %>%
  rowwise() %>%
  mutate(new = max(c(deaths, covid))) %>%
  filter(date <= date_0)

# TODO: Figure out what our actual timeseries of daily deaths is and how to create it:
# e.g. currently the following plot will show isues with say using the max of either deaths or covid above:
# ggplot(data2, aes(date, cumsum(deaths))) + geom_line() +
# geom_line(aes(y = cumsum(covid)), color = "red") +
# geom_line(aes(y = cumsum(new)), color = "blue") +
# geom_line(aes(y = cumsum(only_pos)), color = "green")

## b. Sort out what is to be our death time series
## -----------------------------------------------------------------------------

# here I have just taken the maximum of either excess or covid on each day.
# but this is probably not the best idea as it likely overestimates covid deaths
df <- data2 %>% select(date, new) %>% rename(deaths = new)
df$deaths <- as.integer(df$deaths)

## c. Any other parameters needed to be worked out
## -----------------------------------------------------------------------------

# check that we have this iso3c in squire
if(!(iso3c %in% squire::population$iso3c)) {
  stop("iso3c not found in squire")
}
country <- squire::population$country[match(iso3c, squire::population$iso3c)]
pop <- squire::get_population(country)$n


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
  pars_obs_dur_R = as.numeric(dur_R),
  pars_obs_prob_hosp_multiplier = as.numeric(prob_hosp_multiplier),
  pars_obs_delta_start_date = as.Date(delta_start_date),
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

dp <- dp_plot(res)
cdp <- cdp_plot(res)
ar <- ar_plot(res)

ggsave("fitting.pdf",width=12, height=12,
       cowplot::plot_grid(rtp$plot + ggtitle(country),
                          dp, cdp, ar, ncol = 1))

