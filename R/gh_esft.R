update_gh_esft <- function(origin, destination){
  iso3cs <- list.files(origin, pattern = "[A-Z]")
  iso3cs <- iso3cs[nchar(iso3cs) == 3]

  format_func <- function(df){
    df %>%
      dplyr::filter(
        .data$scenario %in% c(
          "Maintain Status Quo", "Additional 50% Reduction", "Relax Interventions 50%"
        ),
        .data$compartment %in% c(
          "cumulative_infections", "infections", "hospital_incidence", "ICU_incidence", "hospital_demand", "ICU_demand"
        )
      ) %>%
      dplyr::select(.data$scenario, .data$compartment, .data$date,
                    .data$death_calibrated, .data$y_mean) %>%
      suppressMessages()
  }

  purrr::walk(iso3cs, function(iso3c){
    #get projections
    proj <- readr::read_csv(
      paste0(origin, "/", iso3c, "/projections.csv")
    ) %>%
      format_func()
    #for now split up into their groups and save
    dir.create(paste0(destination, "/", iso3c))
    proj <- split(proj, proj$scenario)
    for(x in  names(proj)){
      saveRDS(
        proj[[x]] %>%
          dplyr::select(!.data$scenario),
        paste0(destination, "/", iso3c, "/", x, ".Rds")
      )
    }
  })
  #alternative
  large <- purrr::map_dfr(iso3cs, function(iso3c){
    #get projections
    readr::read_csv(
      paste0(origin, "/", iso3c, "/projections.csv")
    ) %>%
      format_func %>%
      dplyr::mutate(iso3c = iso3c)
  })
  saveRDS(
    large,
    paste0(destination, "/all.Rds")
  )
}
