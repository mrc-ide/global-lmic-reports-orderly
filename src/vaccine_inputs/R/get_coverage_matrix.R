#vectorised version of the old function, removed the scale_cov_mat since that
#doesn't seem to do anything to the model results, essentially vaccines
# are only limited by vaccine uptake and max_vaccine
get_coverage_mats <- function(iso3c, strategy, vaccine_uptake) {
  #get the strategy matrices
  strategies <- readRDS("coverage_strategies.Rds")
  #set up values for loop
  if(length(strategy) == 1){
    strategy <- rep(strategy, length(iso3c))
  }
  if(length(vaccine_uptake) == 1){
    vaccine_uptake <- rep(vaccine_uptake, length(iso3c))
  }
  if(length(strategy) != length(iso3c) |
     length(vaccine_uptake) != length(iso3c)){
    stop("Error: iso3c, strategy and vaccine_uptake must be the same length or have a length of 1.")
  }
  #loop through each country to get the scaled covareg matrix
  cov_mats <- lapply(seq_along(iso3c), function(index){
    if(strategy[index] == "HCW and Elderly") {
      cov_mat <- strategies[[iso3c[index]]]$whoPriority * vaccine_uptake[index]
    } else if (strategy[index] == "HCW, Elderly and High-Risk") {
      cov_mat <- strategies[[iso3c[index]]]$etagePriority * vaccine_uptake[index]
    } else if (strategy[index] == "Elderly") {
      cov_mat <- nimue::strategy_matrix("Elderly", max_coverage = vaccine_uptake[index], 0)
    } else if (strategy[index] == "All") {
      cov_mat <- nimue::strategy_matrix("All", max_coverage = vaccine_uptake[index], 0)
    } else {
      stop('Incorrect strategy. Must be one of "HCW and Elderly", "HCW, Elderly and High-Risk", "Elderly", "All"')
    }
  })
  #give names
  names(cov_mats) <- iso3cs
  return(cov_mats)
}
#A function to calculate the vaccine uptakes for each country
get_vaccine_uptake <- function(iso3cs, dose_df, default_uptake, strategy){
  if(length(strategy) != length(iso3c)){
    stop("Error: iso3c and strategy must be the same length or have a length of 1.")
  }
  #set up defaults
  uptakes <- rep(default_uptake, length(iso3cs))
  #calculate relevant vaccination population
  pop_df <- tibble(
    iso3c = iso3cs,
    strat = strategy
  ) %>%
    left_join(
      squire::population,
      by = "iso3c"
    ) %>%
    group_by(iso3c) %>%
    summarise(
      pop = if_else(
        strat %in% c("All", "Elderly"),
        #use all pop for these strategies
        sum(n),
        #else remove non-adults
        sum(if_else(
          age_group %in% c("0-4", "5-9", "10-14", "15-19"),
          0,
          n
        ))
      )
    ) %>%
    ungroup()
  #check if higher uptake in data
  higher_iso3cs <- dose_df %>%
    group_by(iso3c) %>%
    summarise(vaccinated = sum(first_dose_per_day)) %>%
    left_join(
      pop_df,
      by = "iso3c"
    ) %>%
    mutate( #note some of these are >1 likely to inaccuracy in populations, but these are likely still above 80%
      vaccine_uptake = vaccinated/pop
    ) %>%
    filter(vaccine_uptake > default_uptake) %>%
    pull(iso3c)

  #if uptakes are larger we set the uptake to 95%, the highest we can correctly
  #model or expect to see
  uptakes[iso3cs %in% higher_iso3cs] <- 0.95
  return(uptakes)
}
