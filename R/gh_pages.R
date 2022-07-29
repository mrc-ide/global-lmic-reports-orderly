

update_gh_pages <- function(ids, date){
  repo <- here::here()
  destination <- file.path(
    repo, "gh-pages"
  )

  file_copy <- function(from, to) {
    ok <- file.copy(from, to, overwrite = TRUE, recursive = TRUE)
    if (any(!ok)) {
      stop("There was an error copying files")
    }
  }

 #copy reports over
  is_HIC <- purrr::map_chr(ids, ~.x$iso3c) %>%
    squire.page::get_income_group() %>%
    `%in%`(c("HIC", "UMIC"))
  copy <- c("index.html",
            "projections.csv",
            "index.pdf")
  report_origins <- file.path(repo, "archive", task, purrr::map_chr(ids, ~.x$id))
  report_destinations <- file.path(destination, purrr::map_chr(ids, ~paste0(.x$iso3c, "/", as.character(lubridate::as_date(.x$date)))))
  made_report <- file.size(file.path(report_origins, "index.pdf")) > 0
  purrr::walk(seq_along(report_origins), function(i){
    message(sprintf("Copying %s (%s)", ids[[i]]$iso3c, ids[[i]]$id))
    dir.create(report_destinations[i], FALSE, TRUE)
    if(!made_report[i]){
      copy <- setdiff(copy, c("index.html",
                              "index.pdf"))
    }
    to_copy <- file.path(report_origins[i], copy)
    file_copy(to_copy, report_destinations[i])
    #update main files too
    if(is_HIC[i]){
      to_copy <- to_copy[copy %in% "projections.csv"]
    }
    file_copy(to_copy, dirname(report_destinations[i]))
  })
  #combined reports (LMIC/LIC only)
  pdf_input <- file.path(report_destinations[!is_HIC & made_report], "index.pdf")
  message(sprintf("Building combined pdf from %d files", length(pdf_input)))
  qpdf::pdf_combine(pdf_input, file.path(destination, "combined_reports.pdf"))
  #copy data to regional page
  message("Copying summaries to pages_regional_page")
  summaries <- do.call(rbind,
                       lapply(file.path(report_origins[!is_HIC & made_report], "summary_df.rds"), readRDS))
  saveRDS(summaries, file.path(repo, "src", "pages_regional_page", "summaries.rds"))
  #combine projections
  #legacy projections
  message("Copying legacy projection summaries")
  filename <- file.path(destination, "data", paste0(date,"_v8.csv"))
  purrr::map_dfr(file.path(report_origins, "projections.csv"), ~readr::read_csv(.x, progress = FALSE, show_col_types = FALSE) %>% filter(scenario %in% c(
    "Maintain Status Quo",
    "Relax Interventions 50%",
    "Additional 50% Reduction",
    "Surged Maintain Status Quo",
    "Surged Additional 50% Reduction",
    "Surged Relax Interventions 50%"
  ))) %>%
    dplyr::mutate(version = "v8") %>%
    dplyr::filter(.data$date > min(.data$date) + 250) %>%
    dplyr::ungroup() %>%
    readr::write_csv(filename)
  zip(paste0(filename, ".zip"), filename, extras = '-j')
  file.remove(filename)

  message("Copying projection summaries")
  filename <- file.path(destination, "data", paste0(date,"_v9.csv"))
  projections <- purrr::map_dfr(file.path(report_origins, "projections.csv"), ~readr::read_csv(.x, progress = FALSE, show_col_types = FALSE) %>% filter(scenario %in% c(
    "Central",
    "Optimistic",
    "Pessimistic",
    "Surged Central",
    "Surged Optimistic",
    "Surged Pessimistic"
  ))) %>%
    dplyr::mutate(version = "v9") %>%
    dplyr::filter(.data$date > min(.data$date) + 250) %>%
    dplyr::ungroup()
  readr::write_csv(projections, filename)
  zip(paste0(filename, ".zip"), filename, extras = '-j')
  file.remove(filename)
  #projections to regional page
  saveRDS(projections %>%
            filter(iso3c %in% purrr::map_chr(ids[!is_HIC & made_report], ~.x$iso3c)),
  file.path(repo, "src", "pages_regional_page", "all_data.rds"))
  rm(projections)
  #run orderly tasks
  message("Running pages_ orderly tasks")
  pages_index_page_id <- suppressMessages(orderly::orderly_run("pages_index_page", list(date = date), echo = FALSE))
  orderly::orderly_commit(pages_index_page_id)
  pages_methods_id <- suppressMessages(orderly::orderly_run("pages_methods", list(date = date), echo = FALSE))
  orderly::orderly_commit(pages_methods_id)
  pages_404_id <- suppressMessages(orderly::orderly_run("pages_404", list(date = date), echo = FALSE))
  orderly::orderly_commit(pages_404_id)
  pages_FAQ_id <- suppressMessages(orderly::orderly_run("pages_FAQ", list(date = date), echo = FALSE))
  orderly::orderly_commit(pages_FAQ_id)
  pages_news_id <- suppressMessages(orderly::orderly_run("pages_news", list(date = date), echo = FALSE))
  orderly::orderly_commit(pages_news_id)
  pages_regional_page_ids <- purrr::map_chr(c(Africa = "Africa", Asia = "Asia", Americas = "Americas", Europe = "Europe",
                   Oceania = "Oceania"),
              function(continent){
                id <- suppressMessages(orderly::orderly_run("pages_regional_page", list(date=date, continent = continent)))
                orderly::orderly_commit(id)
                id
              })
  #write data schema
  message("Writing data schema")
  lines <- c(
    "## Data Schema",
    "",
    "* **date** - ISO Date for the predicted number of infection/deaths/hospital burden",
    "* **compartment** - One of deaths, infections, hospital demand and ICU demand or Rt, Reff. Cumulatives are also given as well as prevalence - which is all active infections, which includes hospitilised and non-hopsitilised cases",
    "* **y_025, y_25, y_median, y_mean, y_75, y_975** - Summary statistics for the compartment. E.g. y_25 is the 25% quantile",
    "* **scenario** - The intervention scenario explored. One of Maintain Status Quo, Optimistic, Pessimistic. Surged Maintain Status Quo also present in countries estimated to pass capacity in next 28 days. Additionally old scenarios still output in V8: One of Maintain Status Quo, Additional 50% Reduction, Relax Interventions 50%. Surged Maintain Status Quo also present in countries estimated to pass capacity in next 28 days.",
    "* **country** - Country name",
    "* **iso3c** - Country ISO3C letter",
    "* **report_date** - ISO Date at which the reports were generated, i.e. what is the current date in the dataset",
    "* **version** - Report version. If not present means they were run with verion 1. See Methods page for more details.",
    "* **death_calibrated** - Scenario calibrated to death data. In countries with no deaths to date, we provide scenario projections of what would happen if 5 imported cases occurred today with no future intervention changes. We assume three different R0s. For our upper estimate, we assume that R0 is equal to the mean R0 of countries in the same income classifcation (High, Upper Middle, Lower Middle, Low Income). For our lower estimate, we assume that R0 is equal to either 1.2 or the current Rt of countries in the same income classifcation (whichever is highest). Lastly, our central estimate is half way between the upper and lower estimate. These 3 scenarios are labelled as Relax Interventions 50%, Additional 50% Reduction and Maintain Status Quo in order to maintain consistency with the calibrated scenarios."
  )
  schema_loc <- file.path(destination, "data", "schema.md")
  file.create(schema_loc, showWarnings = FALSE)
  writeLines(lines, schema_loc)
  #copy files over
  message(sprintf("Copying index (%s), methods (%s), 404 (%s), FAQ (%s) and News (%s) pages",
                  pages_index_page_id, pages_methods_id, pages_404_id, pages_FAQ_id, pages_news_id))
  src_index <- file.path(repo, "archive", "pages_index_page", pages_index_page_id,"index.html")
  src_params <- file.path(repo, "archive", "pages_methods", pages_methods_id, "parameters.html")
  src_404 <- file.path(repo, "archive", "pages_404", pages_404_id, "404.html")
  src_FAQ <- file.path(repo, "archive", "pages_FAQ", pages_FAQ_id, "FAQ.html")
  src_News <- file.path(repo, "archive", "pages_news", pages_news_id, "News.html")
  file_copy(c(src_index, src_params, src_404, src_FAQ, src_News), destination)
  message("Copying regional pages")
  src_regional <- file.path(repo, "archive", "pages_regional_page", pages_regional_page_ids)
  dest_regional <- file.path(destination, names(pages_regional_page_ids), date)
  copy_regional <- c("index.html",
            "index.pdf")
  purrr::walk(seq_along(src_regional), function(x){
    message(sprintf("Copying %s (%s)", names(pages_regional_page_ids)[x], pages_regional_page_ids[x]))
    dir.create(dest_regional[x], FALSE, TRUE)
    file_copy(file.path(src_regional[x], copy_regional), dest_regional[x])
    file_copy(file.path(src_regional[x], copy_regional), dirname(dest_regional[x]))
  })
}
