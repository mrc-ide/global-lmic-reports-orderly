#this repo
repo <- here::here()
if(stringr::str_detect(task, "lmic_reports")){
  excess_mortality <- FALSE
  destination <- "standard"
  filename <- "grid_out.Rds"
} else if (stringr::str_detect(task, "excess_mortality")) {
  excess_mortality <- TRUE
  destination <- "excess"
  filename <- "res.Rds"
} else if (task == "reported_deaths_vaccine") {
  excess_mortality <- TRUE
  destination <- "reported"
  filename <- "res.Rds"
}
destination <- file.path(
  repo, "gh-fits",
  destination
)
#get the fits
message("Gathering Fits")
fits <- squire.page::get_fits(repo = repo, date = date, iso3cs = NULL, excess = excess_mortality)

#optionally use task ids
# library(stringr)
# ids <- map_chr(str_split(tasks, "[\\\\.]"), ~tail(.x,2)[1])
# fits <- map(ids, ~readRDS(file.path("archive", task, .x, filename)))
# names(fits) <- names(bundles)


#upload to folder replacing existing files
dir.create(destination, recursive = TRUE, showWarnings = FALSE)
for (iso in names(fits)) {
  saveRDS(fits[[iso]],
          file.path(
            destination,
            paste0(iso, ".Rds")
          ))
}

#update pars inits
message("Updating pars_inits.Rds")
#get old conditions
pars_init <- readRDS(
  file.path("src", task, "pars_init.rds")
  )
for(iso3c in names(fits)) {
  out <- fits[[iso3c]]
  #check if real results (relevant for LMIC fits)
  if("pmcmc_results" %in% names(out)) {
    if("chains" %in% names(out$pmcmc_results)) {
      mc <- do.call(rbind, lapply(out$pmcmc_results$chains, "[[", "results"))
    } else {
      mc <- out$pmcmc_results$results
    }
    best <- mc[which.max(mc$log_posterior),]
    best <- best[, setdiff(colnames(best), c("log_prior", "log_likelihood", "log_posterior"))]
    rownames(best) <- NULL
    best$iso3c <- iso3c
    if(excess_mortality) {
      best$start_date <- as.character(squire:::offset_to_start_date(out$pmcmc_results$inputs$data$week_start[1], round(best$start_date)))

    } else {
      best$start_date <- as.character(squire:::offset_to_start_date(out$pmcmc_results$inputs$data$date[1], round(best$start_date)))

      best$date_Meff_change <- out$pmcmc_results$inputs$Rt_args$date_Meff_change
      best$Rt_shift_duration <- out$pmcmc_results$inputs$Rt_args$Rt_shift_duration
      best$Rt_rw_duration <- out$pmcmc_results$inputs$Rt_args$Rt_rw_duration
    }

    ## Don't use this as drjacoby
    # if("chains" %in% names(out$pmcmc_results)) {
    #   best$covariance_matrix <- out$pmcmc_results$chains$chain1$covariance_matrix[1]
    #   best$scaling_factor <- mean(
    #     c(tail(na.omit(out$pmcmc_results$chains$chain1$scaling_factor),1),
    #       tail(na.omit(out$pmcmc_results$chains$chain2$scaling_factor),1),
    #       tail(na.omit(out$pmcmc_results$chains$chain3$scaling_factor),1))
    #   )
    # } else {
    #   best$covariance_matrix <- out$pmcmc_results$covariance_matrix[1]
    #   best$scaling_factor <- tail(na.omit(out$pmcmc_results$scaling_factor),1)
    # }

    # update pars inits
    pars_init[[iso3c]] <- best
  }
}

saveRDS(
  pars_init,
  file.path("src", task, "pars_init.rds")
)
