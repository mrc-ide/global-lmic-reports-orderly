gather_report_ids <- function(task){
  #copy reports over
  repo <- here::here()
  destination <- file.path(
    repo, "gh-pages"
  )
  #gather latest fit for each country
  db <- orderly::orderly_db("destination",
                            root = repo)

  #just get the latest rt_optimise report, up to this date, for each country
  #I cannot be bothered to relearn SQL so we'll just extract the relevant data
  #and use tidyverse
  DBI::dbGetQuery(db, paste0(
    'SELECT
                                   report_version.id as id,
                                   name as parameter,
                                   value as value
                                 FROM report_version
                                 JOIN parameters
                                   ON report_version.id = parameters.report_version
                                 WHERE report = "', task, '"')) %>%
    dplyr::filter(.data$parameter %in% c("iso3c", "date")) %>%
    tidyr::pivot_wider(names_from = .data$parameter, values_from = .data$value) %>%
    dplyr::mutate(date = lubridate::as_date(.data$date)) %>%
    dplyr::group_by(.data$iso3c) %>%
    dplyr::filter(.data$date == max(.data$date)) %>%
    purrr::transpose()
}
