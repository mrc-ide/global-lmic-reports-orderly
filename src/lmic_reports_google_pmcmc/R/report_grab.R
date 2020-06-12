reports_day <- function(date = NULL) {
  
  db <- orderly::orderly_db("destination")
  if (is.null(date)) {
    date <- as.character(Sys.Date())
  }
  
  ## First find the id corresponding to the ecdc report with data.  If
  ## there are more than one, it's not totally clear what you want to
  ## do as you might want to take the earliest or the latest.
  ## Probably we want to take *all* and do the join over that, which
  ## is easy enough to do if you replace the '= $1' and replace with
  ## 'IN (%s)' and interpolate 'paste(sprintf('"%s"', id), collapse = ", ")'
  sql <- 'SELECT report_version.id
            FROM report_version
            JOIN parameters
              ON parameters.report_version = report_version.id
           WHERE report_version.report = "ecdc"
             AND parameters.value = $1'
  id <- DBI::dbGetQuery(db, sql, date)$id
  if (length(id) == 0L) {
    message(sprintf("No 'ecdc' report for '%s'", as.character(date)))
    return(NULL)
  } else if (length(id) > 1) {
    message(sprintf("Multiple 'ecdc' reports for '%s'", as.character(date)))
  }
  
  ## Then find all lmic_reports reports that use files from this ecdc
  ## report.  This is a bit awful and I might add direct link or a
  ## view to make this easier at some point.
  sql <- 'SELECT report_version.id, parameters.value as country
            FROM report_version_artefact
            JOIN file_artefact
              ON file_artefact.artefact = report_version_artefact.id
            JOIN depends
              ON depends.use = file_artefact.id
            JOIN report_version
              ON report_version.id = depends.report_version
            JOIN parameters
              ON parameters.report_version = report_version.id
           WHERE report_version_artefact.report_version IN (%s)
             AND report = "lmic_reports_google"
             AND parameters.name = "iso3c"
           ORDER BY country, report_version.id'
  sql <- sprintf(sql, paste(sprintf('"%s"', id), collapse = ", "))
  reports <- DBI::dbGetQuery(db, sql)
  
  if(nrow(reports) == 0) {
    DBI::dbDisconnect(db)
    return(NULL)
  } else {
  
  if (any(duplicated(reports$country))) {
    keep <- tapply(seq_len(nrow(reports)), reports$country, max)
    reports <- reports[keep, ]
    rownames(reports) <- NULL
  }
  
  reports$date <- as.character(date)
  DBI::dbDisconnect(db)
  return(reports)
  
  }
}


generate_draws <- function(scan_results, squire_model, replicates, n_particles, forecast,
                           country, population, interventions, data) {
  
  # carry out sims drawn from the grid
  if (is.null(scan_results$z)) {
    res <- squire:::sample_grid_scan(scan_results = scan_results,
                                     n_sample_pairs = replicates,
                                     n_particles = n_particles,
                                     forecast_days = forecast ,
                                     full_output = TRUE)
  } else {
    res <- squire:::sample_3d_grid_scan(scan_results = scan_results,
                                        n_sample_pairs = replicates,
                                        n_particles = n_particles,
                                        forecast_days = forecast ,
                                        full_output = TRUE)
  }
  
  
  # recreate model output for each type of model(ish)
  if (inherits(squire_model, "stochastic")) {
    
    # create a fake run object and fill in the required elements
    r <- squire_model$run_func(country = country,
                               contact_matrix_set = scan_results$inputs$model_params$contact_matrix_set,
                               tt_contact_matrix = scan_results$inputs$model_params$tt_matrix,
                               hosp_bed_capacity = scan_results$inputs$model_params$hosp_bed_capacity,
                               tt_hosp_beds = scan_results$inputs$model_params$tt_hosp_beds,
                               ICU_bed_capacity = scan_results$inputs$model_params$ICU_bed_capacity,
                               tt_ICU_beds = scan_results$inputs$model_params$tt_ICU_beds,
                               population = population,
                               replicates = 1,
                               time_period = nrow(res$trajectories))
    
    # first let's create the output
    names(res)[names(res) == "trajectories"] <- "output"
    dimnames(res$output) <- list(dimnames(res$output)[[1]], dimnames(r$output)[[2]], NULL)
    r$output <- res$output
    
    # and adjust the time as before
    full_row <- match(0, apply(r$output[,"time",],2,function(x) { sum(is.na(x)) }))
    saved_full <- r$output[,"time",full_row]
    for(i in seq_len(replicates)) {
      na_pos <- which(is.na(r$output[,"time",i]))
      full_to_place <- saved_full - which(rownames(r$output) == as.Date(max(data$date))) + 1L
      if(length(na_pos) > 0) {
        full_to_place[na_pos] <- NA
      }
      r$output[,"time",i] <- full_to_place
    }
    
  } else if (inherits(squire_model, "deterministic")) {
    r <- list("output" = res$trajectories)
    r <- structure(r, class = "squire_simulation")
  }
  
  # second let's recreate the output
  r$model <- res$inputs$squire_model$odin_model(
    user = res$inputs$model_params, unused_user_action = "ignore"
  )
  
  # we will add the interventions here so that we now what times are needed for projection
  r$interventions <- interventions
  
  # as well as adding the scan_results so it's easy to draw from the scan again in the future
  r$scan_results <- scan_results
  
  # and add the parameters that changed between each simulation, i.e. drawn from gris
  r$replicate_parameters <- res$param_grid
  
  # and fix the replicates
  r$parameters$replicates <- replicates
  r$parameters$time_period <- as.numeric(diff(as.Date(range(rownames(r$output)))))
  
  return(r)
  
}


generate_draws_pmcmc <- function(pmcmc, burnin, n_chains, squire_model, replicates, n_particles, forecast,
                           country, population, interventions, data) {
  

#--------------------------------------------------------
# Section 3 of pMCMC Wrapper: Sample PMCMC Results
#--------------------------------------------------------
pmcmc_samples <- squire:::sample_pmcmc(pmcmc_results = pmcmc,
                              burnin = burnin,
                              n_chains = n_chains,
                              n_trajectories = replicates,
                              n_particles = n_particles,
                              forecast_days = forecast)

#--------------------------------------------------------
# Section 4 of pMCMC Wrapper: Tidy Output
#--------------------------------------------------------

  # create a fake run object and fill in the required elements
  r <- squire_model$run_func(country = country,
                             contact_matrix_set = pmcmc$inputs$model_params$contact_matrix_set,
                             tt_contact_matrix = pmcmc$inputs$model_params$tt_matrix,
                             hosp_bed_capacity = pmcmc$inputs$model_params$hosp_bed_capacity,
                             tt_hosp_beds = pmcmc$inputs$model_params$tt_hosp_beds,
                             ICU_bed_capacity = pmcmc$inputs$model_params$ICU_bed_capacity,
                             tt_ICU_beds = pmcmc$inputs$model_params$tt_ICU_beds,
                             population = population,
                             day_return = TRUE,
                             replicates = 1,
                             time_period = nrow(pmcmc_samples$trajectories))
  
  # and add the parameters that changed between each simulation, i.e. posterior draws
  r$replicate_parameters <- pmcmc_samples$sampled_PMCMC_Results
  
  # as well as adding the pmcmc chains so it's easy to draw from the chains again in the future
  r$pmcmc_results <- pmcmc
  
  # then let's create the output that we are going to use
  names(pmcmc_samples)[names(pmcmc_samples) == "trajectories"] <- "output"
  dimnames(pmcmc_samples$output) <- list(dimnames(pmcmc_samples$output)[[1]], dimnames(r$output)[[2]], NULL)
  r$output <- pmcmc_samples$output
  
  # and adjust the time as before
  full_row <- match(0, apply(r$output[,"time",],2,function(x) { sum(is.na(x)) }))
  saved_full <- r$output[,"time",full_row]
  for(i in seq_len(replicates)) {
    na_pos <- which(is.na(r$output[,"time",i]))
    full_to_place <- saved_full - which(rownames(r$output) == as.Date(max(data$date))) + 1L
    if(length(na_pos) > 0) {
      full_to_place[na_pos] <- NA
    }
    r$output[,"time",i] <- full_to_place
  }

# second let's recreate the output
r$model <- pmcmc_samples$inputs$squire_model$odin_model(
  user = pmcmc_samples$inputs$model_params, unused_user_action = "ignore"
)

# we will add the interventions here so that we know what times are needed for projection
r$interventions <- interventions

# and fix the replicates
r$parameters$replicates <- replicates
r$parameters$time_period <- as.numeric(diff(as.Date(range(rownames(r$output)))))
r$parameters$dt <- pmcmc$inputs$model_params$dt

return(r)

}

named_list <- function(...) {
  get <- as.character(match.call())
  l <- list(...)
  names(l) <- tail(get, -1)
  return(l)
}


rt_creation <- function(out, date_0, max_date) {
  
  date_0 <- as.Date(date_0)
  
  rts <- lapply(seq_len(nrow(out$replicate_parameters)), function(y) {
    
    tt <- squire:::intervention_dates_for_odin(dates = out$interventions$date_R0_change, 
                                               change = out$interventions$R0_change, 
                                               start_date = out$replicate_parameters$start_date[y],
                                               steps_per_day = 1/out$parameters$dt)
    
    df <- data.frame(
      "Rt" = c(out$replicate_parameters$R0[y], 
               vapply(tt$change, out$pmcmc_results$inputs$Rt_func, numeric(1), 
                      R0 = out$replicate_parameters$R0[y], Meff = out$replicate_parameters$Meff[y])),
      "date" = c(as.character(out$replicate_parameters$start_date[y]), 
                 as.character(out$interventions$date_R0_change[match(tt$change, out$interventions$R0_change)])),
      rep = y,
      stringsAsFactors = FALSE)
    
    if("projection_args" %in% names(out)) {
      
      extra <- data.frame("Rt" = tail(df$Rt, 1) * out$projection_args$R0_change,
                          "date" = as.character(date_0 + 1 + out$projection_args$tt_R0),
                          "rep" = y)
      
      df <- rbind(df, extra)
      
    }
    
    df$pos <- seq_len(nrow(df))
    return(df)
  } )
  
  rt_all <- do.call(rbind, rts)
  
  rt_all$date <- as.Date(rt_all$date)
  rt_all <- rt_all[,c(3,2,1,4)]
  
  new_rt_all <- rt_all %>%
    dplyr::group_by(rep) %>% 
    dplyr::arrange(date) %>% 
    tidyr::complete(date = seq.Date(min(rt_all$date), max_date, by = "days")) 
  
  column_names <- colnames(new_rt_all)[-c(1,2)]
  new_rt_all <- tidyr::fill(new_rt_all, tidyselect::all_of(column_names), .direction = c("down"))
  new_rt_all <- tidyr::fill(new_rt_all, tidyselect::all_of(column_names), .direction = c("up"))
  
  sum_rt <- dplyr::group_by(new_rt_all, date) %>% 
    dplyr::summarise(compartment = "Rt",
                     y_025 = quantile(Rt, 0.025),
                     y_25 = quantile(Rt, 0.25),
                     y_median = median(Rt),
                     y_mean = mean(Rt),
                     y_75 = quantile(Rt, 0.75),
                     y_975 = quantile(Rt, 0.975)) 
  
  head(tail(sum_rt, -1),-1)
}

post_lockdown_date <- function(x, above = 1.1, max_date) {
  
  if(nrow(x)==0) {
    
    return(NA)
    
  } else {
    
    if(any(x$observed)) {
    m <- predict(loess(C~as.numeric(date), data=x, span = 0.2), type = "response")
    } else {
      m <- x$C
    }
    min_mob <- min(m)
    above15 <- which(m >= above*min_mob)
    pl <- above15[which(above15>which.min(m))[1]]
    
    if(x$date[pl] > max_date) {
        min_f <- which(diff(sign(diff(m)))==2)+1
        pl <- min_f[tail(which(x$date[min_f] < max_date),1)]
    }
    
    return(as.Date(x$date[pl])-3)
    
  }

}
