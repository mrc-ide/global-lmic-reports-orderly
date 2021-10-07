plot_waning <- function(adjust_delta, days = 365*2,
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
        get_eff_infection(`Days from Vaccination`, ve_i_second),
      `Protection against Disease,\nfirst dose` =
        get_eff_infection(`Days from Vaccination`, ve_d_first),
      `Protection against Disease,\nsecond dose` =
        get_eff_infection(`Days from Vaccination`, ve_d_second)
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
            get_eff_infection(`Days from Vaccination`, delta_adjustments$ve_i_high_d),
          `Protection against Disease,\nfirst dose` =
            get_eff_infection(`Days from Vaccination`, delta_adjustments$ve_d_low_d),
          `Protection against Disease,\nsecond dose` =
            get_eff_infection(`Days from Vaccination`, delta_adjustments$ve_d_high_d)
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
plot_doses <- function(data, waning){
  if(waning){
    out <- ggplot(data %>%
                    rename(`First Dose` = max_vaccine,
                           `Second Dose` = second_doses
                    ) %>%
                    pivot_longer(c(`First Dose`,
                                   `Second Dose`),
                                 values_to = "Doses given each day",
                                 names_to = "Dose:"),
                  aes(x = Date,
                      y = `Doses given each day`,
                      colour = `Dose:`,
                      linetype = imputed))
  } else {
    out <- ggplot(data %>%
                    rename(`First Doses given` = max_vaccine
                    ),
                  aes(x = Date,
                      y = `First Doses given`,
                      linetype = imputed
                      ))
  }
  out +
    geom_line() + theme_pubclean() +
    scale_linetype(guide = "none")
}
plot_dose_comp <- function(data){
  ggplot(data %>%
           mutate(`First Dose` = INTERNAL_cum_first_in_comp,
                  `Second Dose` = INTERNAL_cum_second_in_comp
           ) %>%
           pivot_longer(c(`First Dose`,
                          `Second Dose`),
                        values_to = "Cumulative Vaccinated,\nadjusted for delay",
                        names_to = "Dose:"),
         aes(x = Date, y = `Cumulative Vaccinated,\nadjusted for delay`, colour = `Dose:`, linetype = imputed)) +
    geom_line() + theme_pubclean() +
    scale_linetype(guide = "none")
}
plot_dose_ratio <- function(data, waning){
  if(waning){
    out <- ggplot(data %>%
                    rename(`Adjusted for` = INTERNAL_dose_ratio_comp,
                           `Raw` = dose_ratio) %>%
                    pivot_longer(c(`Adjusted for`, `Raw`),
                                 names_to = "Delay", values_to = "Dose Ratio"),
                  aes(x = Date, y = `Dose Ratio`, linetype = imputed, colour = Delay))
  } else {
    out <- ggplot(data %>%
             rename(`Dose Ratio` = dose_ratio),
           aes(x = Date, y = `Dose Ratio`, linetype = imputed))
  }
  out +
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
      select(shift_start, shift_end) %>%
      unique() %>%
      pivot_longer(c(shift_start, shift_end), values_to = "x")
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
