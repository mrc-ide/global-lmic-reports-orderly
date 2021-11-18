get_immunity_ratios <- function(out, max_date = NULL) {

  mixing_matrix <- squire:::process_contact_matrix_scaled_age(
    out$pmcmc_results$inputs$model_params$contact_matrix_set[[1]],
    out$pmcmc_results$inputs$model_params$population
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

  betas <- lapply(out$replicate_parameters$R0, function(x) {
    squire:::beta_est(squire_model = out$pmcmc_results$inputs$squire_model,
                      model_params = out$pmcmc_results$inputs$model_params,
                      R0 = x)
  })

  ratios <- lapply(seq_along(betas), function(x) {
    (betas[[x]] * adjusted_eigens[[x]]) / out$replicate_parameters$R0[[x]]
  })

  return(ratios)
}


rt_plot_immunity <- function(out) {

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
      "R0" = na.omit(Rt)[1]*tail(na.omit(ratios[[y]]),length(Rt)),
      "date" = tt$dates,
      rep = y,
      stringsAsFactors = FALSE)
    df$pos <- seq_len(nrow(df))
    return(df)
  } )

  rt <- do.call(rbind, rts)
  rt$date <- as.Date(rt$date)

  rt <- rt[,c(5,4,1,2,3,6)]

  new_rt_all <- rt %>%
    group_by(rep) %>%
    arrange(date) %>%
    complete(date = seq.Date(min(rt$date), date_0, by = "days"))

  column_names <- colnames(new_rt_all)[-c(1,2,3)]
  new_rt_all <- fill(new_rt_all, all_of(column_names), .direction = c("down"))
  new_rt_all <- fill(new_rt_all, all_of(column_names), .direction = c("up"))

  suppressMessages(sum_rt <- group_by(new_rt_all, date) %>%
                     summarise(Rt_min = quantile(Rt, 0.025,na.rm=TRUE),
                               Rt_q25 = quantile(Rt, 0.25,na.rm=TRUE),
                               Rt_q75 = quantile(Rt, 0.75,na.rm=TRUE),
                               Rt_max = quantile(Rt, 0.975,na.rm=TRUE),
                               Rt_median = median(Rt,na.rm=TRUE),
                               Rt = mean(Rt,na.rm=TRUE),
                               R0_min = quantile(R0, 0.025,na.rm=TRUE),
                               R0_q25 = quantile(R0, 0.25,na.rm=TRUE),
                               R0_q75 = quantile(R0, 0.75,na.rm=TRUE),
                               R0_max = quantile(R0, 0.975,na.rm=TRUE),
                               R0_median = median(R0),
                               R0 = mean(R0),
                               Reff_min = quantile(Reff, 0.025,na.rm=TRUE),
                               Reff_q25 = quantile(Reff, 0.25,na.rm=TRUE),
                               Reff_q75 = quantile(Reff, 0.75,na.rm=TRUE),
                               Reff_max = quantile(Reff, 0.975,na.rm=TRUE),
                               Reff_median = median(Reff,na.rm=TRUE),
                               Reff = mean(Reff,na.rm=TRUE)))

  min_date <- min(as.Date(out$replicate_parameters$start_date))

  country_plot <- function(vjust = -1.2) {
    ggplot(sum_rt %>% filter(
      date > min_date & date <= as.Date(as.character(date_0+as.numeric(lubridate::wday(date_0)))))) +
      geom_ribbon(mapping = aes(x=date, ymin=R0_min, ymax = R0_max), fill = "#8cbbca") +
      geom_ribbon(mapping = aes(x = date, ymin = R0_q25, ymax = R0_q75), fill = "#3f8da7") +
      geom_ribbon(mapping = aes(x=date, ymin=Reff_min, ymax = Reff_max), fill = "#96c4aa") +
      geom_ribbon(mapping = aes(x = date, ymin = Reff_q25, ymax = Reff_q75), fill = "#48996b") +
      geom_line(mapping = aes(x = date, y = Reff_median), color = "#48996b") +
      geom_hline(yintercept = 1, linetype = "dashed") +
      geom_hline(yintercept = sum_rt$R0_median[1], linetype = "dashed") +
      theme_bw() +
      theme(axis.text = element_text(size=12)) +
      xlab("") +
      ylab("Reff") +
      scale_x_date(breaks = "2 weeks",
                   limits = as.Date(c(as.character(min_date),
                                      as.character(date_0+as.numeric(lubridate::wday(date_0))))),
                   date_labels = "%d %b",
                   expand = c(0,0)) +
      theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1, colour = "black"),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            panel.border = element_blank(),
            panel.background = element_blank(),
            axis.line = element_line(colour = "black")
      )
  }


  res <- list("plot" = suppressWarnings(country_plot()), "rts" = sum_rt)
  return(res)
}

sero_plot <- function(res, sero_df) {

  # seroconversion data from brazeay report 34
  sero_det <- res$pmcmc_results$inputs$pars_obs$sero_det

  # additional_functions for rolling
  roll_func <- function(x, det) {
    l <- length(det)
    ret <- rep(0, length(x))
    for(i in seq_along(ret)) {
      to_sum <- tail(x[seq_len(i)], length(det))
      ret[i] <- sum(rev(to_sum)*head(det, length(to_sum)), na.rm=TRUE)
    }
    return(ret)
  }


  # get symptom onset data
  date_0 <- max(res$pmcmc_results$inputs$data$date)
  inf <- nim_sq_format(res, c("infections"), date_0 = date_0) %>%
    rename(symptoms = y) %>%
    left_join(nim_sq_format(res, "S", date_0 = date_0),
              by = c("replicate", "t", "date")) %>%
    rename(S = y) %>%
    select(replicate, t, date, S, symptoms)

  inf <- inf %>% group_by(replicate) %>%
    mutate(sero_positive = roll_func(symptoms, sero_det),
           sero_perc = sero_positive/max(S,na.rm = TRUE)) %>%
    group_by(date) %>%
    summarise(sero_perc_med = median(sero_perc, na.rm=TRUE),
              sero_perc_min = quantile(sero_perc, 0.025, na.rm=TRUE),
              sero_perc_max = quantile(sero_perc, 0.975, na.rm=TRUE))

  gg <- ggplot(inf, aes(date, sero_perc_med, ymin = sero_perc_min, ymax = sero_perc_max)) +
    geom_line() +
    geom_ribbon(alpha = 0.2) +
    geom_point(aes(x = date_start + (date_end-date_start)/2, y = sero),
               sero_df, inherit.aes = FALSE) +
    geom_errorbar(aes(x = date_start + (date_end-date_start)/2,
                      ymin = sero_min, ymax = sero_max),
                  sero_df, inherit.aes = FALSE, width = 0) +
    geom_errorbarh(aes(y = sero, xmin = date_start, xmax = date_end),
                   sero_df, inherit.aes = FALSE, height = 0) +
    scale_x_date(date_labels = "%b %Y", date_breaks = "3 months") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    theme_bw() + ylab("Seroprevalence") + xlab("")
  gg

  return(gg)

}

ar_plot <- function(res) {

  S_tot <- sum(res$pmcmc_results$inputs$model_params$population)
  date_0 <- max(res$pmcmc_results$inputs$data$date)
  inf <- nim_sq_format(res, "infections", date_0 = date_0) %>%
    mutate(infections = as.integer(y)) %>%
    select(replicate, t, date, infections) %>%
    group_by(replicate) %>%
    mutate(infections = lag(cumsum(replace_na(infections, 0)), 5, default = 0))


  g2 <- ggplot(inf, aes(date, infections/S_tot, group = replicate)) + geom_line() +
    scale_x_date(date_labels = "%b %Y", date_breaks = "3 months") +
    ylab("Attack Rate") + xlab("") + theme_bw()

  return(g2)
}

cdp_plot <- function(res) {

  suppressWarnings(
    cdp <- plot(res, "D", date_0 = max(res$pmcmc_results$inputs$data$date), x_var = "date") +
      theme_bw() +
      theme(legend.position = "none", axis.title.x = element_blank()) +
      ylab("Cumulative Deaths") +
      scale_x_date(date_labels = "%b %Y", date_breaks = "3 months") +
      xlab("")
  )

  cdp
}

dp_plot <- function(res) {

  suppressWarnings(
    dp <- plot(res, particle_fit = TRUE) +
      theme_bw() +
      theme(legend.position = "none", axis.title.x = element_blank()) +
      scale_x_date(date_labels = "%b %Y", date_breaks = "3 months") +
      xlab("")
  )

  dp

}

get_projections <- function(res, vaccine_inputs, state) {

  ## Functions for working out the relative changes in R0 for given scenarios
  date_0 <- max(res$pmcmc_results$inputs$data$date)
  time_period <- as.integer(as.Date("2022-01-01") - date_0)

  ## We need to know work out vaccine doses and efficacy going forwards
  model_user_args <- extend_vaccine_inputs(vaccine_inputs, time_period, res)
  model_user_args <- lapply(model_user_args, function(x) {
    x$prob_hosp_multiplier <- res$pmcmc_results$inputs$pars_obs$prob_hosp_multiplier
    return(x)
  })

  state <- res$parameters$state
  proj1 <- squire::projections(res, time_period = time_period, model_user_args = model_user_args)
  proj2 <- squire::projections(res, time_period = time_period, R0_change = 1.1, model_user_args = model_user_args)
  proj3 <- squire::projections(res, time_period = time_period, R0_change = 1.25, model_user_args = model_user_args)

  # get data
  get_ret <- function(proj, r_increase) {
    inf <- nim_sq_format(proj, c("infections"), date_0 = date_0) %>%
      rename(symptoms = y) %>%
      left_join(nim_sq_format(proj, "S", date_0 = date_0),
                by = c("replicate", "t", "date")) %>%
      rename(S = y) %>%
      select(replicate, t, date, S, symptoms)

    deaths <- nim_sq_format(proj, c("deaths"), date_0 = date_0) %>%
      rename(deaths = y) %>%
      select(replicate, t, date, deaths)

    ret <- left_join(inf, deaths)
    ret$state <- state
    ret$r_increase <- r_increase

    return(ret)
  }

  rbind(get_ret(proj1, "No Change"),
        get_ret(proj2, "10% Increase"),
        get_ret(proj3, "25% Increase"))

}
