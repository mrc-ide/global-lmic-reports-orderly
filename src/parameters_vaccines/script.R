##Set-up VEs for Vaccine Types
ves_by_type <- read_csv("vaccine_efficacy_groups.csv")
master_ves <- tribble(
  ~Variant, ~Dose, ~Protection, ~Efficacy,
  "Wild", "Partial", "Infection", 0.6,
  "Wild", "Full", "Infection", 0.8,
  "Wild", "Partial", "Hospitalisation", 0.8,
  "Wild", "Full", "Hospitalisation", 0.98,
  "Delta", "Partial", "Infection", 0.224,
  "Delta", "Full", "Infection", 0.646,
  "Delta", "Partial", "Hospitalisation", 0.75,
  "Delta", "Full", "Hospitalisation", 0.94,
  "Omicron", "Partial", "Infection", 0,
  "Omicron", "Full", "Infection", 0.1,
  "Omicron", "Partial", "Hospitalisation", 0.4192211,
  "Omicron", "Full", "Hospitalisation", 0.5254237
)
#calculate omicron changes based on the changes from delta in the old generalised ves
changes <- master_ves %>%
  filter(Variant != "Wild") %>%
  pivot_wider(names_from = Variant, values_from = Efficacy) %>%
  transmute(
    dose = Dose,
    endpoint = Protection,
    p_change = Omicron/Delta,
  )
omicron_efficacies <- ves_by_type %>%
  filter(variant == "Delta") %>%
  left_join(
    changes,
    by = c("dose", "endpoint")
  ) %>%
  mutate(
    variant = "Omicron",
    efficacy = efficacy * p_change
  ) %>%
  select(!p_change)

ves_by_type <- ves_by_type %>%
  rbind(
    omicron_efficacies
  )
##Fit Waning Curves
simulate_time <- 2*365
calc_eff_partial_gen <- odin({
  initial(C1) <- 1
  initial(C2) <- 0
  deriv(C1) <- -w_p*C1
  deriv(C2) <- w_p*C1
  output(ve_d) <- C1 * pV_1_d + C2 * pV_2_d
  output(ve_i) <- C1 * pV_1_i + C2 * pV_2_i
  w_p <- user()
  pV_1_d <- user()
  pV_2_d <- user()
  pV_1_i <- user()
  pV_2_i <- user()
})
calc_eff_full_gen <- odin({
  initial(C1) <- 1
  initial(C2) <- 0
  initial(C3) <- 0
  deriv(C1) <- -w_1*C1
  deriv(C2) <- w_1*C1 - w_2*C2
  deriv(C3) <- w_2*C2
  output(ve_d) <- C1 * fV_1_d + C2 * fV_2_d + C3 * fV_3_d
  output(ve_i) <- C1 * fV_1_i + C2 * fV_2_i + C3 * fV_3_i
  w_1 <- user()
  w_2 <- user()
  fV_1_d <- user()
  fV_2_d <- user()
  fV_3_d <- user()
  fV_1_i <- user()
  fV_2_i <- user()
  fV_3_i <- user()
})
df <- ves_by_type %>% filter(vaccine_type == "mRNA", dose == "Partial", variant == "Wild")
fit_curve <- function(df) {
  parameter_infection <- df %>% filter(endpoint == "Infection") %>% pull(efficacy)
  parameter_hospitalisation <- df %>% filter(endpoint == "Hospitalisation") %>% pull(efficacy)
  dose <- unique(df$dose)
  variant <- unique(df$variant)

  message(paste0(dose, "; ", variant, "; ", unique(df$vaccine_type)))
  #median values for now
  k <- 2.9
  h_s <- 69
  h_l <- 431
  t_s <- 95
  t_l <- 365
  #get scaling
  n_50_d <- 0.027
  n_50_i <- 0.113
  #calculate initial AB value
  ab_d <- c(10^((-log((1/parameter_hospitalisation) - 1)/k) + log10(n_50_d)),
            rep(NA, simulate_time))
  ab_i <- c(10^((-log((1/parameter_infection) - 1)/k) + log10(n_50_i)),
            rep(NA, simulate_time))
  #compute AB curve
  for(t in seq_len(simulate_time)){
    #get deacy rate
    if(t < t_s){
      decay <- 1/h_s
    } else if(t < t_l){
      along <- ((t - t_s)/(t_l - t_s))
      decay <- (1/h_s)*(1-along) + (1/h_l) * along
    } else {
      decay <- 1/h_l
    }
    #update ab
    ab_d[t + 1] <- ab_d[t] * exp(-decay)
    ab_i[t + 1] <- ab_i[t] * exp(-decay)
  }
  #convert to efficacy curve
  ve_d <- 1/(1 + exp(-k*(log10(ab_d) - log10(n_50_d))))
  ve_i <- 1/(1 + exp(-k*(log10(ab_i) - log10(n_50_i))))

  if(dose == "Partial") {
    calc_eff_partial <- calc_eff_partial_gen$new(user = list(
      pV_1_d = parameter_hospitalisation,
      pV_2_d = 0,
      pV_1_i = parameter_infection,
      pV_2_i = 0,
      w_p = 0
    ))
    err_func <- function(pars) {
      #detect if ve is 0
      if(length(pars) < 3){
        pars <- c(pars, 0)
      }

      calc_eff_partial$set_user(
        user = list(
          pV_1_d =  parameter_hospitalisation,
          pV_1_i = parameter_infection,
          w_p = pars[1],
          pV_2_d = pars[2],
          pV_2_i = pars[3]
        )
      )
      mod_value <- calc_eff_partial$run(t = c(0, seq_len(simulate_time)))
      (sum((ve_d - mod_value[, "ve_d"])^2) + sum((ve_i - mod_value[, "ve_i"])^2)) %>% sqrt %>% log
    }
    lower = list(w_p = 1/(3*365), pV_2_d = 0.001, pV_2_i = 0)
    upper = list(w_p = 1/30, pV_2_d = parameter_hospitalisation, pV_2_i = parameter_infection)
    par = list(w_p = 1/365, pV_2_d = parameter_hospitalisation/2, pV_2_i = parameter_infection/2)
    if(parameter_infection == 0){
      lower$pV_2_i <- NULL
      upper$pV_2_i <- NULL
      par$pV_2_i <- NULL
    }
    if((variant == "Wild" & unique(df$vaccine_type) == "mRNA")){
      lower$pV_2_i <- 0.0001
    }
  } else {
    calc_eff_full <- calc_eff_full_gen$new(user = list(
      fV_1_d = parameter_hospitalisation,
      fV_2_d = 0,
      fV_3_d = 0,
      fV_1_i = parameter_infection,
      fV_2_i = 0,
      fV_3_i = 0,
      w_1 = 0,
      w_2 = 0
    ))
    err_func <- function(pars) {
      #detect if ve is 0
      if(length(pars) < 6){
        pars <- c(pars, rep(0, 2))
      }
      calc_eff_full$set_user(
        user = list(
          w_1 = pars[1],
          w_2 = pars[2],
          fV_1_d = parameter_hospitalisation,
          fV_1_i = parameter_infection,
          fV_2_d = pars[3],
          fV_3_d = pars[4],
          fV_2_i = pars[5],
          fV_3_i = pars[6]
        )
      )
      mod_value <- calc_eff_full$run(t = c(0, seq_len(simulate_time)))
      (sum((ve_d - mod_value[, "ve_d"])^2) + sum((ve_i - mod_value[, "ve_i"])^2)) %>%
        sqrt %>% log

    }
    lower = list(
      w_1 = 1/(3*365),
      w_2 = 1/(3*365),
      fV_2_d = 0,
      fV_3_d = 0,
      fV_2_i = 0,
      fV_3_i = 0
    )
    upper = list(
      w_1 = 1/30,
      w_2 = 1/30,
      fV_2_d = parameter_hospitalisation,
      fV_3_d = parameter_hospitalisation,
      fV_2_i = parameter_infection,
      fV_3_i = parameter_infection
    )
    par =  list(
      w_1 = 1/365,
      w_2 = 1/365,
      fV_2_d = parameter_hospitalisation/2,
      fV_3_d = parameter_hospitalisation/2,
      fV_2_i = parameter_infection/2,
      fV_3_i = parameter_infection/2
    )
    if(parameter_infection == 0){
      lower$fV_2_i <- NULL
      upper$fV_2_i <- NULL
      par$fV_2_i <- NULL
      lower$fV_3_i <- NULL
      upper$fV_3_i <- NULL
      par$fV_3_i <- NULL
    }
    if((variant == "Omicron" & unique(df$vaccine_type) == "Johnson&Johnson") |
       (variant %in% c("Delta", "Wild") & unique(df$vaccine_type) == "mRNA")){
      lower$fV_3_i <- 0.00001
    }
    # if(dose == "Full" & variant == "Omicron" & unique(df$vaccine_type) == "Johnson&Johnson"){
    #   lower$fV_2_d <- 0.001
    # } else if(dose == "Full" & variant == "Delta" & unique(df$vaccine_type) == "mRNA"){
    #   lower$fV_3_i <- 0.00001
    # } else if(dose == "Full" & variant == "Wild" & unique(df$vaccine_type) == "mRNA"){
    #   lower$fV_3_i <- 0.00001
    # }
  }
  res <- optim(unlist(par), fn = err_func, method = "L-BFGS-B", lower = lower, upper = upper)
  if(res$convergence != 0){
    stop(res$message)
  }
  if(parameter_infection == 0){
    if(dose == "Partial"){
      res$par[3] <- 0
    } else {
      res$par[6] <- 0
    }
  }
  if(dose == "Full"){
    #ensure is decreasing
    if(res$par[4] > res$par[3]){
      res$par[4] <- res$par[3]
    }
    if(res$par[6] > res$par[5]){
      res$par[6] <- res$par[5]
    }
  }
  #make plot
  if(dose == "Partial"){
    calc_eff_partial$set_user(
      user = list(
        w_p = res$par[1],
        pV_1_d =  parameter_hospitalisation,
        pV_1_i = parameter_infection,
        pV_2_d = res$par[2],
        pV_2_i = res$par[3]
      )
    )
    mod_value <- calc_eff_partial$run(t = c(0, seq_len(simulate_time)))
    ve_f_d <- mod_value[, "ve_d"]
    ve_f_i <- mod_value[, "ve_i"]
  } else {
    calc_eff_full$set_user(
      user = list(
        w_1 = res$par[1],
        w_2 = res$par[2],
        fV_1_d = parameter_hospitalisation,
        fV_1_i = parameter_infection,
        fV_2_d = res$par[3],
        fV_3_d = res$par[4],
        fV_2_i = res$par[5],
        fV_3_i = res$par[6]
      )
    )
    mod_value <- calc_eff_full$run(t = c(0, seq_len(simulate_time)))
    ve_f_d <- mod_value[, "ve_d"]
    ve_f_i <- mod_value[, "ve_i"]
  }
  p <- ggplot(
    tibble(
      t = rep(c(0, seq_len(simulate_time)), 4),
      `Protection:` = c(rep("Disease", (simulate_time + 1)*2), rep("Infection", (simulate_time + 1)*2)),
      `Version:` = rep(c(rep("AB Process", simulate_time + 1), rep("Booster Model", simulate_time + 1)), 2),
      `Vaccine Efficacy` = c(ve_d, ve_f_d, ve_i, ve_f_i)
    )
    , aes(x = t, y = `Vaccine Efficacy`, colour = `Protection:`, linetype = `Version:`)
  ) +
    geom_line() +
    labs(y = "Vaccine Efficacy", x = "Days Since Dose", title = paste0("Dose: ", dose, ", Variant: ", variant, ", Type: ", df$vaccine_type[1])) +
    ggpubr::theme_pubclean()

  out <- as.data.frame(res$par)
  out <- mutate(out, parameter = rownames(out),
                parameter = if_else(parameter == "", "pV_2_i", parameter)) %>%
    rename(value = `res$par`) %>%
    rbind(tibble(
      value = c(parameter_hospitalisation, parameter_infection),
      parameter = if_else(rep(dose == "Full", 2), c("fV_1_d", "fV_1_i"), c("pV_1_d", "pV_1_i"))
    ))
  rownames(out) <- NULL
  out$dose <- dose
  out$variant <- variant
  out$platform <- unique(df$vaccine_type)
  #return parameters
  return(list(out, p))
}
set.seed(1000100001)
values <-
  ves_by_type %>%
  group_by(vaccine_type, dose, variant) %>%
  group_split() %>%
  map(fit_curve)

#split into plots and data
plots <- map(values, ~.x[[2]])
efficacies <- map(values, ~.x[[1]])

#add randomness (ideally we should add this at the start but fitting process is too sensitive)
N_samples <- 100
random_efficacies <- map(efficacies, function(x){
  message(paste0(x$dose[1], "; ", x$variant[1], "; ", x$platform[1]))
  map_dfr(seq_len(N_samples), function(it){
    gamma_var <- 10
    beta_var <- 0.005
    percentage_zero <- 0.5
    min_value <- beta_var * 1.1
    zero_inf_beta <- function(n, alpha, beta, p_0){
      zero <- as.numeric(runif(n) < p_0)
      rbeta(n, alpha, beta)*abs(zero - 1)
    }
    calculate_alpha_beta <- function(mean){
      mean <- if_else(mean < min_value, min_value, mean)
      (mean)*((mean)*(1-(mean))/beta_var - 1)
    }
    calculate_beta_beta <- function(mean, alpha){
      mean <- if_else(mean < min_value, min_value, mean)
      (alpha - mean*alpha)/mean
    }
    samples <- x %>%
      mutate(
        zero_inflation_percentage = if_else(value == 0, percentage_zero, 0),
        alpha = case_when(
          parameter %in% c("w_1", "w_2", "w_p") ~ (1/value) * (1/value)/gamma_var,
          value == 0 ~ calculate_alpha_beta(value),
          TRUE ~ calculate_alpha_beta(value)
        ),
        beta = case_when(
          parameter %in% c("w_1", "w_2", "w_p") ~ (1/value)/gamma_var,
          value == 0 ~ calculate_beta_beta(value, alpha),
          TRUE ~ calculate_beta_beta(value, alpha)
        ),
        distribution = case_when(
          parameter %in% c("w_1", "w_2", "w_p") ~ "gamma",
          value == 0 ~ "zero inflated beta",
          TRUE ~ "beta"
        ),
        value = case_when(
          parameter %in% c("w_1", "w_2", "w_p") ~ 1/rgamma(n(), alpha, beta),
          value == 0 ~ zero_inf_beta(n(), alpha, beta, zero_inflation_percentage),
          TRUE ~ rbeta(n(), alpha, beta)
        ),
        sample = it
      )
    #correct so not increasing in full dose protection,ISSUE: should adjust for this
    #earlier
    samples %>%
      mutate(
        temp = str_remove(parameter, "\\d")
      ) %>%
      group_by(temp) %>%
      arrange(temp, parameter) %>%
      mutate(
        value = if_else(
          !str_detect(temp, "w") & value < lead(value, 1, default = 0),
          lead(value, 1, default = 0),
          value
        )
      ) %>%
      ungroup() %>%
      select(!temp)
  })
})

##Calibration plots

pdf("calibration.pdf")
map(seq_along(plots), function(x){
  plot <- plots[[x]]
  #generate curves for random efficacies
  generate_curve <- function(df){
    user <- as.list(pull(df, value, parameter))
    if (df$dose[1] == "Full") {
      mod <- calc_eff_full_gen$new(user = user)
    } else {
      mod <- calc_eff_partial_gen$new(user = user)
    }
    tt <- c(0, seq_len(simulate_time))
    mod$run(t = c(0, seq_len(simulate_time))) %>%
      as_tibble() %>%
      select(t, ve_d, ve_i) %>%
      mutate(sample = df$sample[1],
             t = as.integer(t),
             ve_d = as.numeric(ve_d),
             ve_i = as.numeric(ve_i))
  }
  random_bounds <- random_efficacies[[x]] %>%
    group_by(sample) %>%
    group_split %>%
    map_dfr(generate_curve) %>%
    group_by(t) %>%
    summarise(
      ve_d_025 = quantile(ve_d, 0.025),
      ve_d_975 = quantile(ve_d, 0.975),
      ve_i_025 = quantile(ve_i, 0.025),
      ve_i_975 = quantile(ve_i, 0.975),
      .groups = "drop"
    ) %>%
    pivot_longer(!t, names_to = "temp", values_to = "y") %>%
    mutate(
      `Protection:` = if_else(str_detect(temp, "_i_"), "Infection", "Disease"),
      bound = if_else(str_detect(temp, "_025"), "lower", "upper")
    ) %>%
    select(!temp) %>%
    pivot_wider(names_from = bound, values_from = y)
  plot + geom_ribbon(data = random_bounds, aes(
    x = t, ymin = lower , ymax = upper, fill = `Protection:`
  ), inherit.aes = FALSE, alpha = 0.1)
})
dev.off()

##Create Sampling Function

#add omicon sub unit values
random_efficacies <- map_dfr(random_efficacies, ~.x)
random_efficacies <- random_efficacies %>%
  rbind(
    random_efficacies %>% filter(platform == "Johnson&Johnson") %>% mutate(dose = "Partial", value = 0)
  ) %>%
  rbind(
    random_efficacies %>% filter(variant == "Omicron") %>% mutate(variant = "Omicron Sub Variant")
  ) %>%
  arrange(
    sample, platform, variant, dose, parameter
  )

#empty environment of everything not relevant
rm(list = setdiff(ls(), c("N_samples", "random_efficacies")))

sample_vaccine_efficacies <- function(n, platforms){
  #uniformly samples
  platforms <- sample(platforms, n, replace = TRUE)
  #now for each plat form draw a random sample from the df
  output <- map_dfr(seq_along(platforms),
      ~random_efficacies %>% filter(platform == platforms[.x], sample == sample(N_samples, 1)) %>% select(dose, variant, parameter, value),
      .id = "sample")
  #make per variant and format into model compatible formats
  variants <- unique(output$variant)
  output <- map(variants, function(this_variant){
    df <- filter(output, variant == this_variant)
    samples <- map(unique(df$sample), function(this_sample){
      values <- as.list(filter(df, sample == this_sample) %>%
        pull(value, parameter))
      #convert to format
      list(
        dur_V = 1/c(values$w_1, values$w_2, values$w_p),
        vaccine_efficacy_infection = c(values$pV_1_i, values$pV_2_i, values$fV_1_i, values$fV_2_i, values$fV_3_i),
        vaccine_efficacy_disease = c(values$pV_1_d, values$pV_2_d, values$fV_1_d, values$fV_2_d, values$fV_3_d)
      )
    })
    #reformat so easier to use
    output <- map(names(samples[[1]]), function(x){map(samples, ~.x[[x]])})
    names(output) <- names(samples[[1]])
    output
  })
  names(output) <- variants
  output
}

saveRDS(list(sample_vaccine_efficacies = sample_vaccine_efficacies), "vaccine_params.Rds")
