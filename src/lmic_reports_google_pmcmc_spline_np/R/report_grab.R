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

generate_draws_pmcmc_fitted <- function(out, pmcmc, burnin, n_chains, squire_model, replicates, n_particles, forecast,
                                        country, population, interventions, data) {
  
  
  #--------------------------------------------------------
  # Section 1 # what is our predicted gradient
  #--------------------------------------------------------
  
  infections <- format_output(out, "infections", date_0 = max(data$date))
  infections <- infections %>% filter(date > (max(data$date) - 30) & date < (max(data$date))) %>% 
    group_by(date) %>% summarise(y = median(y))
  
  get_grad <- function(x) {
    lm(y~x, data = data.frame(y = x, x = seq_along(x)))$coefficients[2]
  }
  
  pred_grad <- get_grad(infections$y)
  des_grad <- get_grad(tail(data$cases, 30))
  
  index <- squire:::odin_index(out$model)
  
  if(sign(pred_grad) != sign(des_grad)) {
  
  # do we need to go up or down
  if(des_grad <= pred_grad) {
    alters <- seq(0, 0.4, 0.05)
  } else {
    alters <- seq(0, -0.2, -0.05) # more conservative on increasing Rt
  }
  
  # store our grads
  ans <- alters
  last_rw <- ncol(out$pmcmc_results$chains$chain1$results) - 3
  
  # for later
  all_case_compartments <- unlist(
    index[c("IMild", "ICase1", "ICase2", "IOxGetLive1", "IOxGetLive2",
            "IOxGetDie1", "IOxGetDie2", "IOxNotGetLive1", "IOxNotGetLive2",
            "IOxNotGetDie1", "IOxNotGetDie2", "IMVGetLive1", "IMVGetLive2",
            "IMVGetDie1", "IMVGetDie2", "IMVNotGetLive1", "IMVNotGetLive2",
            "IMVNotGetDie1", "IMVNotGetDie2", "IRec1", "IRec2", "R", "D")])
  
  #--------------------------------------------------------
  # Section 2 # # find best grad correction
  #--------------------------------------------------------
  
  for(alt in seq_along(alters)) {
    
    for(ch in seq_along(out$pmcmc_results$chains)) {
      out$pmcmc_results$chains[[ch]]$results[,last_rw] <- out$pmcmc_results$chains[[ch]]$results[,last_rw] + alters[alt]
    }
    
    pmcmc_samples <- squire:::sample_pmcmc(pmcmc_results = out$pmcmc_results,
                                           burnin = burnin,
                                           n_chains = n_chains,
                                           n_trajectories = replicates,
                                           n_particles = n_particles,
                                           forecast_days = forecast)
    
    nt <- nrow(pmcmc_samples$trajectories)
    
    # assign the infections
    for(i in seq_along(out$parameters$population)) {
      collect <- vapply(1:dim(pmcmc_samples$trajectories)[3], function(j) {
        pos <- seq(i,length(all_case_compartments), by = length(out$parameters$population))
        pos <- all_case_compartments[pos]
        diff(rowSums(pmcmc_samples$trajectories[,pos,j]))
      }, FUN.VALUE = numeric(nt-1))
      pmcmc_samples$trajectories[1+seq_len(nt-1),index$n_E2_I[i],] <- collect
    }
    
    
    this_infs <- as.numeric(rowMeans(tail(matrix(unlist(lapply(seq_len(replicates), function(i) {
      rowSums(pmcmc_samples$trajectories[,index$n_E2_I,i])
    })), ncol = replicates),30)))
    
    ans[alt] <- get_grad(this_infs)
    
    
    for(ch in seq_along(out$pmcmc_results$chains)) {
      out$pmcmc_results$chains[[ch]]$results[,last_rw] <- out$pmcmc_results$chains[[ch]]$results[,last_rw] - alters[alt]
    }
    
  }
  
  
  # adapt our whole last chain accordingly
  alts <- which.min(abs(ans-des_grad))
  for(ch in seq_along(out$pmcmc_results$chains)) {
    out$pmcmc_results$chains[[ch]]$results[,last_rw] <- out$pmcmc_results$chains[[ch]]$results[,last_rw] + alters[alts]
  }
  
  }
  
  # set up now to do the stochastic draws
  out$pmcmc_results$inputs$squire_model <- squire:::explicit_model()
  out$pmcmc_results$inputs$model_params$dt <- 0.02
  
  #--------------------------------------------------------
  # Section 3 of pMCMC Wrapper: Sample PMCMC Results
  #--------------------------------------------------------
  pmcmc_samples <- squire:::sample_pmcmc(pmcmc_results = out$pmcmc_results,
                                         burnin = burnin,
                                         n_chains = n_chains,
                                         n_trajectories = replicates,
                                         n_particles = n_particles,
                                         forecast_days = forecast)
  
  #--------------------------------------------------------
  # Section 4 of pMCMC Wrapper: Tidy Output
  #--------------------------------------------------------
  
  # create a fake run object and fill in the required elements
  r <- out$pmcmc_results$inputs$squire_model$run_func(country = out$parameters$country,
                             contact_matrix_set = out$pmcmc_results$inputs$model_params$contact_matrix_set,
                             tt_contact_matrix = out$pmcmc_results$inputs$model_params$tt_matrix,
                             hosp_bed_capacity = out$pmcmc_results$inputs$model_params$hosp_bed_capacity,
                             tt_hosp_beds = out$pmcmc_results$inputs$model_params$tt_hosp_beds,
                             ICU_bed_capacity = out$pmcmc_results$inputs$model_params$ICU_bed_capacity,
                             tt_ICU_beds = out$pmcmc_results$inputs$model_params$tt_ICU_beds,
                             population = out$pmcmc_results$inputs$population,
                             replicates = 1,
                             day_return = TRUE,
                             time_period = nrow(pmcmc_samples$trajectories))
  
  # and add the parameters that changed between each simulation, i.e. posterior draws
  r$replicate_parameters <- pmcmc_samples$sampled_PMCMC_Results

  # as well as adding the pmcmc chains so it's easy to draw from the chains again in the future
  r$pmcmc_results <- out$pmcmc_results
  
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
  r$parameters$dt <- out$pmcmc_results$inputs$model_params$dt
  
  return(r)
  
}


named_list <- function(...) {
  get <- as.character(match.call())
  l <- list(...)
  names(l) <- tail(get, -1)
  return(l)
}


rt_creation <- function(out, date_0, max_date) {
  
  iso3c <- squire::get_population(out$parameters$country)$iso3c[1]
  
  if("pmcmc_results" %in% names(out)) {
    wh <- "pmcmc_results"
  } else {
    wh <- "scan_results"
  }
  
  date <- max(as.Date(out$pmcmc_results$inputs$data$date))
  date_0 <- date
  
  # impact of immunity ratios
  ratios <- get_immunity_ratios(out)
  
  # create the Rt data frame
  rts <- lapply(seq_len(length(out$replicate_parameters$R0)), function(y) {
    
    tt <- squire:::intervention_dates_for_odin(dates = out$interventions$date_R0_change, 
                                               change = out$interventions$R0_change, 
                                               start_date = out$replicate_parameters$start_date[y],
                                               steps_per_day = 1/out$parameters$dt)
    
    if(wh == "scan_results") {
      Rt <- c(out$replicate_parameters$R0[y], 
              vapply(tt$change, out[[wh]]$inputs$Rt_func, numeric(1), 
                     R0 = out$replicate_parameters$R0[y], Meff = out$replicate_parameters$Meff[y])) 
    } else {
      Rt <- squire:::evaluate_Rt_pmcmc(
        R0_change = tt$change, 
        date_R0_change = tt$dates, 
        R0 = out$replicate_parameters$R0[y], 
        pars = as.list(out$replicate_parameters[y,]),
        Rt_args = out$pmcmc_results$inputs$Rt_args) 
    }
    
    df <- data.frame(
      "Rt" = Rt,
      "Reff" = Rt*tail(na.omit(ratios[[y]]),length(Rt)),
      "date" = tt$dates,
      "iso" = iso3c,
      rep = y,
      stringsAsFactors = FALSE)
    df$pos <- seq_len(nrow(df))
    return(df)
  } )
  
  rt <- do.call(rbind, rts)
  rt$date <- as.Date(rt$date)
  
  rt <- rt[,c(3,4,1,2,5,6)]
  
  new_rt_all <- rt %>%
    group_by(iso, rep) %>% 
    arrange(date) %>% 
    complete(date = seq.Date(min(rt$date), date_0, by = "days")) 
  
  column_names <- colnames(new_rt_all)[-c(1,2,3)]
  new_rt_all <- fill(new_rt_all, all_of(column_names), .direction = c("down"))
  new_rt_all <- fill(new_rt_all, all_of(column_names), .direction = c("up"))
  
  sum_rt <- dplyr::group_by(new_rt_all, date) %>% 
    dplyr::summarise(compartment = "Rt",
                     y_025 = quantile(Rt, 0.025),
                     y_25 = quantile(Rt, 0.25),
                     y_median = median(Rt),
                     y_mean = mean(Rt),
                     y_75 = quantile(Rt, 0.75),
                     y_975 = quantile(Rt, 0.975)) 
  sum_rt <- head(sum_rt,-1)
  
  sum_reff <- dplyr::group_by(new_rt_all, date) %>% 
    dplyr::summarise(compartment = "Reff",
                     y_025 = quantile(Reff, 0.025),
                     y_25 = quantile(Reff, 0.25),
                     y_median = median(Reff),
                     y_mean = mean(Reff),
                     y_75 = quantile(Reff, 0.75),
                     y_975 = quantile(Reff, 0.975)) 
  sum_reff <- head(sum_reff,-1)
  
  return(rbind(sum_rt, sum_reff))
}

post_lockdown_date <- function(x, above = 1.1, max_date, min_date) {
  
  if(nrow(x)==0) {
    
    return(NA)
    
  } else {
    
    if(any(x$observed)) {
      m <- predict(loess(C~as.numeric(date), data=x, span = 0.2), type = "response")
    } else {
      m <- x$C
    }
    min_mob <- min(m)
    
    pl <- NA
    while(is.na(pl)) {
      above15 <- which(m >= above*min_mob)
      pl <- above15[which(above15>which.min(m))[1]]
      above <- above*0.99
    }
    
    # if past max date then take the last local  minimum and grow by 4 days
    if(x$date[pl] > max_date) {
      min_f <- which(diff(sign(diff(m)))==2)+1
      pl <- min_f[tail(which(x$date[min_f] < max_date),1)] + 4
    } 
    
    dat <- max(min_date, as.Date(x$date[pl]))
    
    return(dat)
    
  }
  
}

post_lockdown_date_relative <- function(x, above = 1.1, max_date, min_date) {
  
  if(nrow(x)==0) {
    
    return(NA)
    
  } else {
    
    if(any(x$observed)) {
      m <- predict(loess(C~as.numeric(date), data=x, span = 0.2), type = "response")
    } else {
      m <- predict(loess(C~as.numeric(date), data=x, span = 0.2), type = "response")
      #m <- x$C
    }
    min_mob <- min(m)
    diff <- 1 - min_mob
    
    pl <- NA
    while(is.na(pl)) {
      above15 <- which(m >= ((above-1)*diff)+min_mob)
      pl <- above15[which(above15>which.min(m))[1]]
      above <- above*0.99
    }
    
    # if past max date then take the last local  minimum and grow by 4 days
    if(x$date[pl] > max_date) {
      min_f <- which(diff(sign(diff(m)))==2)+1
      pl <- min_f[tail(which(x$date[min_f] < max_date),1)] + 4
    } 
    
    dat <- max(min_date, as.Date(x$date[pl]))
    
    return(dat)
    
  }
  
}

ede <- function (x, y, index) {
  n = length(x)
  if (index == 1) {
    y = -y
  }
  ifelse(n >= 4, {
    LF = y - (y[1] + (y[n] - y[1]) * (x - x[1])/(x[n] - x[1]))
    jf1 = which.min(LF)
    xf1 = x[jf1]
    jf2 = which.max(LF)
    xf2 = x[jf2]
    ifelse(jf2 < jf1, {
      xfx <- NaN
    }, {
      xfx <- 0.5 * (xf1 + xf2)
    })
  }, {
    jf1 = NaN
    jf2 = NaN
    xfx = NaN
  })
  out = matrix(c(jf1, jf2, xfx), nrow = 1, ncol = 3, byrow = TRUE)
  rownames(out) = "EDE"
  colnames(out) = c("j1", "j2", "chi")
  return(out)
}



mobility_decline_date <- function(x, above = 1.1, max_date, min_date) {
  
  message(x$iso3c[1])
  
  if(nrow(x)==0) {
    
    return(NA)
    
  } else {
    
    max_date <- as.Date(max_date)
    min_date <- as.Date(min_date)
    
    x <- x[x$date > min_date,]
    
    if(any(x$observed)) {
      
      get <- x[,c("date","C")]
      names(get) <- c("time","intensity")
      
    } else {
      
      m <- predict(loess(C~as.numeric(date), data=x, span = 0.2), type = "response")
      get <- data.frame("time" = x$date, "intensity" = m)
      
    }
    
    get$intensity <- (get$intensity - min(get$intensity))/(max(get$intensity)-min(get$intensity))
    get$time <- as.numeric(get$time)-as.numeric(get$time[1])
    
    fitObj <- sicegar::fitAndCategorize(get,
                                        threshold_minimum_for_intensity_maximum = 0.3,
                                        threshold_intensity_range = 0.2,
                                        threshold_t0_max_int = 0.2)
    
    
    dat <- max(min_date, x$date[round(fitObj$doubleSigmoidalModel$startDeclinePoint_x)])
    
    return(dat)
    
  }
  
}

