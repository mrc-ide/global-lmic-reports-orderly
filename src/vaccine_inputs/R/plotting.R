plot_waning <- function(adjust_delta, days = 365*2,
                        days_between_doses = days_between_doses,
                        ve_i_first = ve_i_low, ve_i_second = ve_i_high,
                        ve_d_first = ve_d_low, ve_d_second = ve_d_high){
  #main data
  waning_data <- tibble(
    `Days from Vaccination` = seq(1, 365*2)
  ) %>%
    mutate(
      Delta = FALSE,
      `Protection against Infection,\nfirst dose` =
        get_eff_infection(`Days from Vaccination`, ve_i_first),
      `Protection against Infection,\nsecond dose` =
        get_eff_infection_second(`Days from Vaccination`, days_between_doses,
                               ve_i_first, ve_i_second),
      `Protection against Disease,\nfirst dose` =
        get_eff_disease(`Days from Vaccination`, ve_d_first),
      `Protection against Disease,\nsecond dose` =
        get_eff_disease_second(`Days from Vaccination`, days_between_doses,
                               ve_d_first, ve_d_second)
    )
  if(adjust_delta){
    delta_adjustments <- readRDS("delta_characteristics.Rds") %>%
      select(contains("ve")) %>%
      unique()
    waning_data <- waning_data %>%  rbind(
      tibble(`Days from Vaccination` = waning_data$`Days from Vaccination`) %>%
        mutate(
          Delta = TRUE,
          `Protection against Infection,\nfirst dose` =
            get_eff_infection(`Days from Vaccination`, delta_adjustments$ve_i_low_d),
          `Protection against Infection,\nsecond dose` =
            get_eff_infection_second(`Days from Vaccination`, days_between_doses,
                                   delta_adjustments$ve_i_low_d,
                                   delta_adjustments$ve_i_high_d),
          `Protection against Disease,\nfirst dose` =
            get_eff_disease(`Days from Vaccination`, delta_adjustments$ve_d_low_d),
          `Protection against Disease,\nsecond dose` =
            get_eff_disease_second(`Days from Vaccination`, days_between_doses,
                                   delta_adjustments$ve_d_low_d,
                                   delta_adjustments$ve_d_high_d)
        )
    )
    waning_plots <- lapply(names(waning_data)[-c(1,2)], function(var){
      ggplot(waning_data) +
        geom_line(
          aes(x = `Days from Vaccination`, y = .data[[var]], colour = Delta),
          alpha = 0.75
        ) + theme_pubclean()
    })
  } else {
    waning_plots <- lapply(names(waning_data)[-c(1,2)], function(var){
      ggplot(waning_data) +
        geom_line(
          aes(x = `Days from Vaccination`, y = .data[[var]]),
          alpha = 0.75
        ) + theme_pubclean()
    })
  }
  ggarrange(
    plotlist = waning_plots,
    common.legend = TRUE,
    legend = "right"
  )
}
plot_doses <- function(data){
  ggplot(data %>%
                    rename(`First Doses given` = max_vaccine
                    ),
                  aes(x = Date,
                      y = `First Doses given`,
                      linetype = imputed
                      )) +
    geom_line() + theme_pubclean() +
    scale_linetype(guide = "none")
}
plot_dose_ratio <- function(data){
  ggplot(data %>%
             rename(`Dose Ratio` = dose_ratio),
           aes(x = Date, y = `Dose Ratio`, linetype = imputed))+
    geom_line() +
    theme_pubclean() + ylim(c(0,1)) +
    scale_linetype(guide = "none")
}
plot_mean_first_date <- function(data){
  ggplot(data %>%
           mutate(`Assumed days\nbetween doses` =as.numeric(Date - INTERNAL_mean_first_dose_date)),
         aes(x = Date, y = `Assumed days\nbetween doses`,
             linetype = imputed)) +
    geom_point() +
    theme_pubclean() +
    scale_linetype(guide = "none")
}
plot_efficacy <- function(data, adjust_delta){
  out <- ggplot(data %>%
           rename(`Infection`=vaccine_efficacy_infection,
                  `Disease`=vaccine_efficacy_disease) %>%
           pivot_longer(c(`Infection`, `Disease`), names_to = "Protection:",
                        values_to = "Effective VE"),
         aes(x = Date, y = `Effective VE`,
             linetype = imputed, colour = `Protection:`)) +
    geom_step() +
    theme_pubclean() +
    scale_linetype(guide = "none") + ylim(c(0,1))
  if(adjust_delta){
    delta_values <- data %>%
      ungroup() %>%
      select(delta_shift_start, delta_shift_end) %>%
      unique() %>%
      pivot_longer(c(delta_shift_start, delta_shift_end), values_to = "x")
    delta_values <- delta_values %>%
      rbind(delta_values) %>%
      arrange(x) %>%
      mutate(y = c(Inf, -Inf, -Inf, Inf))
    out <- out +
      geom_polygon(inherit.aes = FALSE, data = delta_values, aes(
        x = x,
        y = y
      ), alpha = 0.05) +
      geom_vline(data = delta_values, aes(xintercept = x),
                 linetype = "dashed")
  }
  out
}
plot_efficacy_split <- function(data, adjust_delta){
  plot_1 <- ggplot(this_country %>%
                     mutate(
                       `First` = INTERNAL$vaccine_efficacy_infection_first,
                       `Second` = INTERNAL$vaccine_efficacy_infection_second
                     ) %>%
                     pivot_longer(c(`First`, `Second`),
                                  names_to = "Dose:",
                                  values_to = "Protection against\nInfection"),
                   aes(x = Date, y = `Protection against\nInfection`, colour = `Dose:`)) +
    geom_line() +
    theme_pubclean() +
    scale_linetype(guide = "none") + ylim(c(0,1))
  plot_2 <-
    ggplot(this_country %>%
             mutate(
               `First` = INTERNAL$vaccine_efficacy_disease_first,
               `Second` = INTERNAL$vaccine_efficacy_disease_second
             ) %>%
             pivot_longer(c(`First`, `Second`),
                          names_to = "Dose:",
                          values_to = "Protection against\nDisease"),
           aes(x = Date, y = `Protection against\nDisease`, colour = `Dose:`)) +
    geom_line() +
    theme_pubclean() +
    scale_linetype(guide = "none") + ylim(c(0,1))
  if(adjust_delta){
      delta_values <- data %>%
        ungroup() %>%
        select(shift_start, shift_end) %>%
        unique() %>%
        pivot_longer(c(shift_start, shift_end), values_to = "x")
      delta_values <- delta_values %>%
        rbind(delta_values) %>%
        arrange(x) %>%
        mutate(y = c(Inf, -Inf, -Inf, Inf))
      plot_1 <- plot_1 +
        geom_polygon(inherit.aes = FALSE, data = delta_values, aes(
          x = x,
          y = y
        ), alpha = 0.05) +
        geom_vline(data = delta_values, aes(xintercept = x),
                   linetype = "dashed")
      plot_2 <- plot_2 +
        geom_polygon(inherit.aes = FALSE, data = delta_values, aes(
          x = x,
          y = y
        ), alpha = 0.05) +
        geom_vline(data = delta_values, aes(xintercept = x),
                   linetype = "dashed")
    }
  suppressWarnings(ggarrange(
    plot_1, plot_2,
    ncol = 1,
    common.legend = TRUE
  ))
}
