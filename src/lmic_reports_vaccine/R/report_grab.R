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

generate_draws_pmcmc_fitted <- function(out, n_particles = 10, grad_dur = 21) {
  
  pmcmc <- out$pmcmc_results
  n_chains <- length(out$pmcmc_results$chains)
  burnin <- round(out$pmcmc_results$inputs$n_mcmc/10)
  squire_model <- out$pmcmc_results$inputs$squire_model
  replicates <- dim(out$output)[3]
  forecast <- 0
  country <- out$parameters$country
  population <- out$parameters$population
  interventions <- out$interventions
  data <- out$pmcmc_results$inputs$data
  rw_dur <- out$pmcmc_results$inputs$Rt_args$Rt_rw_duration
  
  #--------------------------------------------------------
  # Section 1 # what is our predicted gradient
  #--------------------------------------------------------
  
  # first what is the model predicted infections
  infections <- format_output(out, "infections", date_0 = max(data$date))
  infections_end <- infections %>% filter(date > (max(data$date) - grad_dur) & date <= (max(data$date))) %>% 
    group_by(date) %>% summarise(y = median(y))
  
  infections_pre_end <- infections %>% 
    filter(date > (max(data$date) - grad_dur - rw_dur) & date <= (max(data$date) - grad_dur) ) %>% 
    group_by(date) %>% summarise(y = median(y))
  
  # and the observed cases
  cases_end <- tail(data$cases, grad_dur)
  cases_pre_end <- head(tail(data$cases, grad_dur+rw_dur), rw_dur)
  
  # get these gradients
  get_grad <- function(x) {
    x[x==0] <- NA
    # only get a grad if more than half data is there
    if(sum(is.na(x)) > 0.5*length(x)) {
      return(NA)
    } else {
    lm(log(y)~x, data = data.frame(y = x, x = seq_along(x)))$coefficients[2]
    }
  }
  
  pred_grad_end <- get_grad(infections_end$y)
  pred_grad_pre_end <- get_grad(infections_pre_end$y)
  
  des_grad_end <- get_grad(cases_end)
  des_grad_pre_end <- get_grad(cases_pre_end)
  
  # if the cases are just not good enough then don't
  if(!is.na(des_grad_end) && !is.na(des_grad_pre_end) && country != "Indonesia") {
  
  if(sign(pred_grad_pre_end) == sign(des_grad_pre_end)) {
  
    ca_grad_frac <-  pred_grad_pre_end / des_grad_pre_end
    
  } else {
    
    infections_pre_end <- infections %>% 
      filter(date <= (max(data$date) - grad_dur) ) %>% 
      group_by(date) %>% summarise(y = median(y))
    
    cases_pre_end <- data$cases[match(infections_pre_end$date, data$date)]
    
    pos <- seq_len(floor(length(cases_pre_end)/rw_dur)*rw_dur)
    
    breaks <- split(pos,
          sort(unlist(replicate(floor(length(cases_pre_end)/rw_dur),
                                seq_len(rw_dur),simplify = FALSE))))
    
    na_to_0 <- function(x) {x[is.na(x)] <- 0; return(x)}
    
    inf_grads <- lapply(breaks, function(x) {get_grad(na_to_0(infections_pre_end$y[x]))})
    case_grads <- lapply(breaks, function(x) {get_grad(na_to_0(cases_pre_end[x]))})
    ca_grad_frac <- median((unlist(inf_grads)/unlist(case_grads)), na.rm = TRUE)
    
    # if we still can't get a gradient then don't adjust from the fitting (is risky)
    if(is.na(ca_grad_frac) | is.infinite(ca_grad_frac)) {
      ca_grad_frac <- NA
    }
    
  }
  
  # desired model predictd final gradient
  wanted_grad <- des_grad_end * ca_grad_frac
  
  # if actual gradient available and there are at least half non 0s
  if(!is.nan(wanted_grad) || !is.na(wanted_grad) || !is.infinite(wanted_grad) &&
     (sum(cases_end == 0)/length(cases_end)) < 0.5) {
  
  index <- squire:::odin_index(out$model)
  index$n_E2_I <- seq(tail(unlist(index),1)+1, tail(unlist(index),1)+length(index$S),1)
  index$delta_D <- seq(tail(unlist(index),1)+1, tail(unlist(index),1)+length(index$S),1)
  
  # do we need to go up or down
  if(wanted_grad < pred_grad_end) {
    alters <- seq(0.025, 0.175, 0.025)
  } else {
    alters <- seq(-0.025, -0.125, -0.025) # more conservative on the way up
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
    
    message(alt)  
    
    for(ch in seq_along(out$pmcmc_results$chains)) {
      out$pmcmc_results$chains[[ch]]$results[,last_rw] <- out$pmcmc_results$chains[[ch]]$results[,last_rw] + alters[alt]
    }
    
    pmcmc_samples <- squire:::sample_pmcmc(pmcmc_results = out$pmcmc_results,
                                           burnin = burnin,
                                           n_chains = n_chains,
                                           n_trajectories = replicates,
                                           n_particles = n_particles,
                                           forecast_days = forecast)
    
    dimnms <- dimnames(pmcmc_samples$trajectories)
    
    # make e2_i space
    dimnms[[2]] <- c(dimnms[[2]], paste0("n_E2_I[", seq_len(length(index$S)),"]"))
    new_data <- array(data = 0,
                      dim = c(dim(pmcmc_samples$trajectories) + c(0, length(index$n_E2_I), 0)),
                      dimnames = dimnms)
    
    new_data[, seq_len(dim(pmcmc_samples$trajectories)[2]), ] <- pmcmc_samples$trajectories
    pmcmc_samples$trajectories <- new_data
    
    # make D space
    dimnms <- dimnames(pmcmc_samples$trajectories)
    dimnms[[2]] <- c(dimnms[[2]], paste0("delta_D[", seq_len(length(index$S)),"]"))
    new_data <- array(data = 0,
                      dim = c(dim(pmcmc_samples$trajectories) + c(0, length(index$delta_D), 0)),
                      dimnames = dimnms)
    
    new_data[, seq_len(dim(pmcmc_samples$trajectories)[2]), ] <- pmcmc_samples$trajectories
    pmcmc_samples$trajectories <- new_data
    nt <- nrow(pmcmc_samples$trajectories)
    
    # are the steps not 1 apart? if so we need to sum the incident variables (infecions/deaths)
      if (out$parameters$day_return || !squire:::odin_is_discrete(out$model)) {
        
        # assign the infections
        for(i in seq_along(out$parameters$population)) {
          collect <- vapply(1:out$parameters$replicates, function(j) {
            pos <- seq(i, length(index$cum_infs), by = length(out$parameters$population))
            pos <- index$cum_infs[pos]
            diff(pmcmc_samples$trajectories[,pos,j])
          }, FUN.VALUE = numeric(nt-1))
          pmcmc_samples$trajectories[1+seq_len(nt-1),index$n_E2_I[i],] <- collect
        }
        
        # assign the deaths
        for(i in seq_along(out$parameters$population)) {
          collect <- vapply(1:out$parameters$replicates, function(j) {
            pos <- seq(i, length(index$D), by = length(out$parameters$population))
            pos <- index$D[pos]
            diff(pmcmc_samples$trajectories[,pos,j])
          }, FUN.VALUE = numeric(nt-1))
          pmcmc_samples$trajectories[1+seq_len(nt-1),index$delta_D[i],] <- collect
        }
        
      }
    
    this_infs <- as.numeric(rowMeans(tail(matrix(unlist(lapply(seq_len(replicates), function(i) {
      rowSums(pmcmc_samples$trajectories[,index$n_E2_I,i])
    })), ncol = replicates),grad_dur)))
    
    ans[alt] <- get_grad(this_infs)
    
    
    for(ch in seq_along(out$pmcmc_results$chains)) {
      out$pmcmc_results$chains[[ch]]$results[,last_rw] <- out$pmcmc_results$chains[[ch]]$results[,last_rw] - alters[alt]
    }
    
  }
  
  
  # adapt our whole last chain accordingly
  alts <- which.min(abs(ans-wanted_grad))
  for(ch in seq_along(out$pmcmc_results$chains)) {
    out$pmcmc_results$chains[[ch]]$results[,last_rw] <- out$pmcmc_results$chains[[ch]]$results[,last_rw] + alters[alts]
  }
  
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
                                                      time_period = nrow(pmcmc_samples$trajectories),
                                                      dur_R = out$pmcmc_results$inputs$model_params$dur_R)
  
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

generate_draws_pmcmc_case_fitted <- function(out, n_particles = 10, grad_dur = 21) {
  
  pmcmc <- out$pmcmc_results
  n_chains <- length(out$pmcmc_results$chains)
  burnin <- round(out$pmcmc_results$inputs$n_mcmc/10)
  squire_model <- out$pmcmc_results$inputs$squire_model
  replicates <- dim(out$output)[3]
  forecast <- 0
  country <- out$parameters$country
  population <- out$parameters$population
  interventions <- out$interventions
  data <- out$pmcmc_results$inputs$data
  rw_dur <- out$pmcmc_results$inputs$Rt_args$Rt_rw_duration
  
  #--------------------------------------------------------
  # Section 1 # what is our predicted gradient
  #--------------------------------------------------------
  
  # first what is the model predicted infections
  infections <- format_output(out, "infections", date_0 = max(data$date))
  infections_end <- infections %>% filter(date > (max(data$date) - grad_dur) & date <= (max(data$date))) %>% 
    group_by(date) %>% summarise(y = median(y))
  
  infections_pre_end <- infections %>% 
    filter(date > (max(data$date) - grad_dur - rw_dur) & date <= (max(data$date) - grad_dur) ) %>% 
    group_by(date) %>% summarise(y = median(y))
  
  # and the observed cases
  cases_end <- tail(data$cases, grad_dur)
  cases_pre_end <- head(tail(data$cases, grad_dur+rw_dur), rw_dur)
  
  # get these gradients
  get_infs <- function(x) {
    sum(x, na.rm = TRUE)
  }
  
  pred_infs_end <- get_infs(infections_end$y)
  pred_infs_pre_end <- get_infs(infections_pre_end$y)
  
  des_infs_end <- get_infs(cases_end)
  des_infs_pre_end <- get_infs(cases_pre_end)
  
  # if there are less than 100 cases in both windowns then don't bother
  if(des_infs_end > 100 && des_infs_pre_end > 100) {
      
    ca_infs_frac <-  pred_infs_pre_end / des_infs_pre_end
      
    # desired model predictd final infs
    wanted_infs <- des_infs_end * ca_infs_frac
    
    # if actual infs available
    if(!is.nan(wanted_infs) || !is.na(wanted_infs) || !is.infinite(wanted_infs)) {
      
      index <- squire:::odin_index(out$model)
      index$n_E2_I <- seq(tail(unlist(index),1)+1, tail(unlist(index),1)+length(index$S),1)
      index$delta_D <- seq(tail(unlist(index),1)+1, tail(unlist(index),1)+length(index$S),1)
      
      # do we need to go up or down
      if(wanted_infs < pred_infs_end) {
        alters <- seq(0.025, 0.175, 0.025)
      } else {
        alters <- seq(-0.025, -0.125, -0.025) # more conservative on the way up
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
        
        message(alt)  
        
        for(ch in seq_along(out$pmcmc_results$chains)) {
          out$pmcmc_results$chains[[ch]]$results[,last_rw] <- out$pmcmc_results$chains[[ch]]$results[,last_rw] + alters[alt]
        }
        
        pmcmc_samples <- squire:::sample_pmcmc(pmcmc_results = out$pmcmc_results,
                                               burnin = burnin,
                                               n_chains = n_chains,
                                               n_trajectories = replicates,
                                               n_particles = n_particles,
                                               forecast_days = forecast)
        
        dimnms <- dimnames(pmcmc_samples$trajectories)
        
        # make e2_i space
        dimnms[[2]] <- c(dimnms[[2]], paste0("n_E2_I[", seq_len(length(index$S)),"]"))
        new_data <- array(data = 0,
                          dim = c(dim(pmcmc_samples$trajectories) + c(0, length(index$n_E2_I), 0)),
                          dimnames = dimnms)
        
        new_data[, seq_len(dim(pmcmc_samples$trajectories)[2]), ] <- pmcmc_samples$trajectories
        pmcmc_samples$trajectories <- new_data
        
        # make D space
        dimnms <- dimnames(pmcmc_samples$trajectories)
        dimnms[[2]] <- c(dimnms[[2]], paste0("delta_D[", seq_len(length(index$S)),"]"))
        new_data <- array(data = 0,
                          dim = c(dim(pmcmc_samples$trajectories) + c(0, length(index$delta_D), 0)),
                          dimnames = dimnms)
        
        new_data[, seq_len(dim(pmcmc_samples$trajectories)[2]), ] <- pmcmc_samples$trajectories
        pmcmc_samples$trajectories <- new_data
        nt <- nrow(pmcmc_samples$trajectories)
        
        # are the steps not 1 apart? if so we need to sum the incident variables (infecions/deaths)
        if (out$parameters$day_return || !squire:::odin_is_discrete(out$model)) {
          
          # assign the infections
          for(i in seq_along(out$parameters$population)) {
            collect <- vapply(1:out$parameters$replicates, function(j) {
              pos <- seq(i, length(index$cum_infs), by = length(out$parameters$population))
              pos <- index$cum_infs[pos]
              diff(pmcmc_samples$trajectories[,pos,j])
            }, FUN.VALUE = numeric(nt-1))
            pmcmc_samples$trajectories[1+seq_len(nt-1),index$n_E2_I[i],] <- collect
          }
          
          # assign the deaths
          for(i in seq_along(out$parameters$population)) {
            collect <- vapply(1:out$parameters$replicates, function(j) {
              pos <- seq(i, length(index$D), by = length(out$parameters$population))
              pos <- index$D[pos]
              diff(pmcmc_samples$trajectories[,pos,j])
            }, FUN.VALUE = numeric(nt-1))
            pmcmc_samples$trajectories[1+seq_len(nt-1),index$delta_D[i],] <- collect
          }
          
        }
        
        this_infs <- as.numeric(rowMeans(tail(matrix(unlist(lapply(seq_len(replicates), function(i) {
          rowSums(pmcmc_samples$trajectories[,index$n_E2_I,i])
        })), ncol = replicates),grad_dur)))
        
        ans[alt] <- get_infs(this_infs)
        
        
        for(ch in seq_along(out$pmcmc_results$chains)) {
          out$pmcmc_results$chains[[ch]]$results[,last_rw] <- out$pmcmc_results$chains[[ch]]$results[,last_rw] - alters[alt]
        }
        
      }
      
      
      # adapt our whole last chain accordingly
      alts <- which.min(abs(ans-wanted_infs))
      for(ch in seq_along(out$pmcmc_results$chains)) {
        out$pmcmc_results$chains[[ch]]$results[,last_rw] <- out$pmcmc_results$chains[[ch]]$results[,last_rw] + alters[alts]
      }
      
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
                                                      time_period = nrow(pmcmc_samples$trajectories),
                                                      dur_R = out$pmcmc_results$inputs$model_params$dur_R)
  
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
  } else if("scan_results" %in% names(out)) {
    wh <- "scan_results"
  } else {
    wh <- "simple"
  }
  
  if (wh != "simple") {
  
  date <- max(as.Date(out$pmcmc_results$inputs$data$date))
  date_0 <- date
  
  # impact of immunity ratios
  ratios <- get_immunity_ratios(out, max_date)
  
  # create the Rt data frame
  rts <- lapply(seq_len(length(out$replicate_parameters$R0)), function(y) {
    
    dates <- c(out$interventions$date_R0_change, 
               seq.Date(as.Date(tail(out$interventions$date_R0_change,1)) + 1,
                        as.Date(max_date),
                        1))
    
    change <- c(out$interventions$R0_change, 
                rep(tail(out$interventions$R0_change,1), length(dates)-length(out$interventions$R0_change)))
    
    tt <- squire:::intervention_dates_for_odin(dates = dates, 
                                               change = change, 
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
  
  sum_reff <- dplyr::group_by(new_rt_all, date) %>% 
    dplyr::summarise(compartment = "Reff",
                     y_025 = quantile(Reff, 0.025),
                     y_25 = quantile(Reff, 0.25),
                     y_median = median(Reff),
                     y_mean = mean(Reff),
                     y_75 = quantile(Reff, 0.75),
                     y_975 = quantile(Reff, 0.975)) 
  
  ret <- rbind(sum_rt, sum_reff)
  
  } else {
    ret <- rt_creation_simple_out(out, date_0, max_date)
  }
  
  return(ret)
}

rt_creation_simple_out <- function(out, date_0, max_date) {
  
  iso3c <- squire::get_population(out$parameters$country)$iso3c[1]
  
  get_immunity_ratios_simple <- function(out, max_date) {
    
    mixing_matrix <- squire:::process_contact_matrix_scaled_age(
      out$parameters$contact_matrix_set[[1]],
      out$parameters$population
    )
    
    dur_ICase <- out$parameters$dur_ICase
    dur_IMild <- out$parameters$dur_IMild
    prob_hosp <- out$parameters$prob_hosp
    
    # assertions
    squire:::assert_single_pos(dur_ICase, zero_allowed = FALSE)
    squire:::assert_single_pos(dur_IMild, zero_allowed = FALSE)
    squire:::assert_numeric(prob_hosp)
    squire:::assert_numeric(mixing_matrix)
    squire:::assert_square_matrix(mixing_matrix)
    squire:::assert_same_length(mixing_matrix[,1], prob_hosp)
    
    if(sum(is.na(prob_hosp)) > 0) {
      stop("prob_hosp must not contain NAs")
    }
    
    if(sum(is.na(mixing_matrix)) > 0) {
      stop("mixing_matrix must not contain NAs")
    }
    
    index <- squire:::odin_index(out$model)
    pop <- out$parameters$population
    if(is.null(max_date)) {
      max_date <- max(out$pmcmc_results$inputs$data$date)
    }
    t_now <- which(as.Date(rownames(out$output)) == max_date)
    prop_susc <- lapply(seq_len(dim(out$output)[3]), function(x) {
      t(t(out$output[seq_len(t_now), index$S, x])/pop)
    } )
    
    relative_R0_by_age <- prob_hosp*dur_ICase + (1-prob_hosp)*dur_IMild
    
    adjusted_eigens <- lapply(prop_susc, function(x) {
      
      unlist(lapply(seq_len(nrow(x)), function(y) {
        if(any(is.na(x[y,]))) {
          return(NA)
        } else {
          Re(eigen(mixing_matrix*x[y,]*relative_R0_by_age)$values[1])
        }
      }))
      
    })
    
    
    mat <- squire:::process_contact_matrix_scaled_age(out$parameters$contact_matrix_set[[1]],
                                                      out$parameters$population)
    
    betas <- lapply(rep(out$parameters$R0, dim(out$output)[3]), function(x) {
      squire:::beta_est_explicit(dur_IMild = dur_IMild, 
                                 dur_ICase = dur_ICase, 
                                 prob_hosp = prob_hosp, 
                                 mixing_matrix = mat, 
                                 R0 = x)
    })
    
    ratios <- lapply(seq_along(betas), function(x) {
      (betas[[x]] * adjusted_eigens[[x]]) / out$parameters$R0
    })
    
    return(ratios)
  }
  
  # impact of immunity ratios
  ratios <- get_immunity_ratios_simple(out, max_date)
  
  # create the Rt data frame
  rts <- lapply(seq_len(dim(out$output)[3]), function(y) {
    
    Rt <- rep(out$parameters$R0, length(ratios[[y]]))
    
    df <- data.frame(
      "Rt" = Rt,
      "Reff" = Rt*tail(na.omit(ratios[[y]]),length(Rt)),
      "date" = head(rownames(out$output),length(Rt)),
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
  
  sum_reff <- dplyr::group_by(new_rt_all, date) %>% 
    dplyr::summarise(compartment = "Reff",
                     y_025 = quantile(Reff, 0.025),
                     y_25 = quantile(Reff, 0.25),
                     y_median = median(Reff),
                     y_mean = mean(Reff),
                     y_75 = quantile(Reff, 0.75),
                     y_975 = quantile(Reff, 0.975)) 
  
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
    attempts <- 50
    while(is.na(pl) && attempts == 0) {
      above15 <- which(m >= ((above-1)*diff)+min_mob)
      pl <- above15[which(above15>which.min(m))[1]]
      above <- above*0.99
      attempts <- attempts -1
    }
    
    # if still NA then just return max date
    if(is.na(pl)) {
      return(max_date)
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

extend_df_for_covidsim <- function(df, out, ext = 240) {
  
  betas <- df$beta_set
  tt_R0 <- df$tt_beta
  dates <- df$date
  
  # get an initial with the same seeds as before
  population <- squire::get_population(country = out$parameters$country, simple_SEIR = FALSE)
  init <- squire:::init_check_explicit(NULL, population$n, seeding_cases = 5)
  init$S <- init$S + init$E1
  init$E1 <- c(0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0)
  init$S <- init$S - c(0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0)
  
  # run model
  det_out <- squire::run_deterministic_SEIR_model(
    country = out$parameters$country,
    beta_set = betas,
    walker_params = FALSE,
    init = init,
    day_return = TRUE,
    tt_R0 = tt_R0+1,
    R0 = tt_R0,
    time_period = length(dates) + ext)
  
  # summarise the changes on the susceptibles
  index <- squire:::odin_index(det_out$model)
  
  # get the ratios
  mixing_matrix <- squire:::process_contact_matrix_scaled_age(
    out$pmcmc_results$inputs$model_params$contact_matrix_set[[1]],
    out$pmcmc_results$inputs$model_params$population
  )
  dur_ICase <- out$parameters$dur_ICase
  dur_IMild <- out$parameters$dur_IMild
  prob_hosp <- out$parameters$prob_hosp
  pop <- out$parameters$population
  
  prop_susc <- t(t(det_out$output[, index$S, 1])/pop)
  relative_R0_by_age <- prob_hosp*dur_ICase + (1-prob_hosp)*dur_IMild
  
  adjusted_eigens <-  unlist(lapply(seq_len(nrow(prop_susc)), function(y) {
    if(any(is.na(prop_susc[y,]))) {
      return(NA)
    } else {
      Re(eigen(mixing_matrix*prop_susc[y,]*relative_R0_by_age)$values[1])
    }
  }))
  
  betas <- squire:::beta_est(squire_model = out$pmcmc_results$inputs$squire_model, 
                             model_params = out$pmcmc_results$inputs$model_params, 
                             R0 = df$Rt[1])
  
  ratios <- (betas * adjusted_eigens) / df$Rt[1]
  
  # extend df
  df2 <- df %>% 
    complete(tt_beta = seq(0, max(df$tt_beta) + ext, 1)) %>% 
    mutate(date = seq.Date(min(df$date), max(df$date)+ext,1)) %>% 
    fill(c("beta_set",4:10), .direction = "down") %>% 
    mutate(Reff = ratios*Rt)
  
  return(df2)
  
}
