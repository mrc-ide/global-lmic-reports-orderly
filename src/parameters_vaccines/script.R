setwd("/home/gregbarnsley/Documents/imperial/cepi/new_fitting/global-lmic-reports-orderly/")
orderly::orderly_develop_start("parameters_vaccines")
setwd("/home/gregbarnsley/Documents/imperial/cepi/new_fitting/global-lmic-reports-orderly/src/parameters_vaccines")
##Set-up VEs for Vaccine Types
ves_by_type <- read_csv("vaccine_efficacy_groups.csv") %>%
  mutate(dose = if_else(dose == "Partial", "First", "Second"))
#assumed numbers, booster restores efficacy NOTE Ideally we'd find better numbers, might be too many categories
master_ves <- tribble(
  ~Variant, ~Dose, ~Endpoint, ~Efficacy,
  "Wild", "First", "Infection", 0.6,
  "Wild", "Second", "Infection", 0.8,
  "Wild", "Booster", "Infection", 0.85,
  "Wild", "First", "Hospitalisation", 0.8,
  "Wild", "Second", "Hospitalisation", 0.98,
  "Wild", "Booster", "Hospitalisation", 0.98,
  "Delta", "First", "Infection", 0.224,
  "Delta", "Second", "Infection", 0.646,
  "Delta", "Booster", "Infection", 0.85,
  "Delta", "First", "Hospitalisation", 0.75,
  "Delta", "Second", "Hospitalisation", 0.94,
  "Delta", "Booster", "Hospitalisation", 0.98,
  "Omicron", "First", "Infection", 0,
  "Omicron", "Second", "Infection", 0.1,
  "Omicron", "Booster", "Infection", 0.55,
  "Omicron", "First", "Hospitalisation", 0.1557433,
  "Omicron", "Second", "Hospitalisation", 0.449152542372881,
  "Omicron", "Booster", "Hospitalisation", 0.90
)

#calculate omicron changes based on the changes from delta in the old generalised ves
logit <- function(x){
  log(x/(1-x))
}
inv_logit <- function(x){
  1/(1+exp(-x))
}
#assume that boosters are either mRNA or whol virus
changes_booster <- master_ves %>%
  pivot_wider(names_from = Dose, values_from = Efficacy) %>%
  transmute(
    variant = Variant,
    endpoint = Endpoint,
    p_change = logit(Booster) - logit(Second)
  )
booster_efficacies <- ves_by_type %>%
  filter(dose == "Second" & vaccine_type %in% c("mRNA", "Whole Virus")) %>%
  left_join(
    changes_booster,
    by = c("variant", "endpoint")
  ) %>%
  mutate(
    dose = "Booster",
    efficacy = inv_logit(logit(efficacy) + p_change)
  ) %>%
  select(!p_change)
ves_by_type <- ves_by_type %>%
  rbind(
    booster_efficacies
  )
changes_omicron <- master_ves %>%
  filter(Variant != "Wild") %>%
  pivot_wider(names_from = Variant, values_from = Efficacy) %>%
  transmute(
    dose = Dose,
    endpoint = Endpoint,
    p_change = Omicron/Delta,
  )
omicron_efficacies <- ves_by_type %>%
  filter(variant == "Delta") %>%
  left_join(
    changes_omicron,
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
  ) %>%
  arrange(vaccine_type, variant, endpoint, dose)

#overwrite with omicron data for boosters where possible
ves_by_type <- ves_by_type %>%
  mutate(
    start_period = 7, #what period of time does this VE cover (assume 1-5 weeks)
    end_period = 7*6,
    efficacy = case_when(
      vaccine_type == "mRNA" & dose == "booster" & endpoint == "Infection" & variant == "Omicron" ~
        0.65, #https://www.nejm.org/doi/full/10.1056/NEJMoa2119451
      vaccine_type == "mRNA" & dose == "booster" & endpoint == "Hospitalisation" & variant == "Omicron" ~
        0.82, #https://www.sciencedirect.com/science/article/pii/S0264410X22005230
      vaccine_type == "Whole Virus" & dose == "booster" & endpoint == "Hospitalisation" & variant == "Omicron" ~
        0.80, #evidence suggest that its similar to mRNA so we'll keep them the same
      #https://doi.org/10.1016/S2214-109X(22)00112-7 #Not the right variant but presumably will be lower?
      #some evidence that it's much higher https://www.medrxiv.org/content/10.1101/2022.03.22.22272769v1#:~:text=Two%20doses%20of%20either%20vaccine,%3A%2067.8%25%2C%2079.2%25).
      TRUE ~ efficacy
    )
  )
#subvariant data where possible

##Fit Waning Curves
simulate_time <- 2*365
calc_eff_gen <- odin({
  initial(C1) <- 1
  initial(C2) <- 0
  initial(C3) <- 0
  deriv(C1) <- -w_1*C1
  deriv(C2) <- w_1*C1 - w_2*C2
  deriv(C3) <- w_2*C2
  output(ve_d) <- (C1 * ved + C2 * ved_2 + C3 * ved_3)
  output(ve_i) <- C1 * vei
  w_1 <- user()
  w_2 <- user()
  ved <- user()
  vei <- user()
  ved_2 <- user()
  ved_3 <- user()
})
fit_curve <- function(df) {
  parameter_infection <- df %>% filter(endpoint == "Infection") %>% pull(efficacy)
  parameter_hospitalisation <- df %>% filter(endpoint == "Hospitalisation") %>% pull(efficacy)
  dose <- unique(df$dose)
  variant <- unique(df$variant)
  platform <- unique(df$vaccine_type)

  message(paste0(dose, "; ", variant, "; ", platform))

  if(dose == "First"){
    out <- tibble(
      parameter = c("pV_i", "pV_d"),
      value = c(parameter_infection, parameter_hospitalisation),
      dose = dose,
      variant = variant,
      platform = platform
    )
    p <- NULL
  } else {
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

    calc_eff <- calc_eff_gen$new(user = list(
      ved = parameter_hospitalisation,
      vei = parameter_infection,
      ved_2 = 1,
      ved_3 = 1,
      w_1 = 0,
      w_2 = 0
    ))
    err_lines <- function(l1, l2){
      sum((l1 - l2)^2/l1)
      #scale it so the lower values have more weight
    }
    err_func <- function(pars) {
      calc_eff$set_user(
        user = list(
          w_1 = pars[1],
          w_2 = pars[2],
          ved_2 = pars[3] * parameter_hospitalisation,
          ved_3 = pars[3] * parameter_hospitalisation * pars[4],
          ved = parameter_hospitalisation,
          vei = parameter_infection
        )
      )
      mod_value <- calc_eff$run(t = c(0, seq_len(simulate_time)))

      plot((ve_d - mod_value[, "ve_d"])^2/ve_d)

      log(sqrt(err_lines(ve_d, mod_value[, "ve_d"]) + err_lines(ve_i, mod_value[, "ve_i"])))

    }
    lower = list(
      w_1 = 1/(3*365),
      w_2 = 1/(3*365),
      ved_2 = 0,
      ved_3 = 0
    )
    upper = list(
      w_1 = 1/30,
      w_2 = 1/30,
      ved_2 = 1,
      ved_3 = 1
    )
    par =  list(
      w_1 = 1/365,
      w_2 = 1/365,
      ved_2 = 0.5,
      ved_3 = 0.5
    )
    res <- dfoptim::nmkb(unlist(par), fn = err_func, lower = unlist(lower), upper = unlist(upper))
    if(res$convergence != 0){
      stop(res$message)
    }
    
    V_d_2 <- res$par[3] * parameter_hospitalisation
    V_d_3 <- res$par[3] * parameter_hospitalisation * res$par[4]

    #make plot
    calc_eff$set_user(
      user = list(
        w_1 = res$par[1],
        w_2 = res$par[2],
        ved_2 = V_d_2,
        ved_3 = V_d_3,
        ved = parameter_hospitalisation,
        vei = parameter_infection
      )
    )
    mod_value <- calc_eff$run(t = c(0, seq_len(simulate_time)))
    ve_f_d <- mod_value[, "ve_d"]
    ve_f_i <- mod_value[, "ve_i"]
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
    out <- mutate(out, parameter = c("w_1", "w_2", "V_2_d", "V_3_d")) %>%
      rename(value = `res$par`) %>%
      mutate(
        value = case_when(
          parameter == "V_2_d" ~ V_d_2,
          parameter == "V_3_d" ~ V_d_3,
          TRUE ~ value
        )
      ) %>% 
      rbind(tibble(
        value = c(parameter_hospitalisation, parameter_infection),
        parameter = c("V_1_d", "V_i")
      )) %>% 
      mutate(
        parameter = if_else(rep(dose == "Booster", 6), paste0("b", parameter), paste0("f", parameter))
      )
    rownames(out) <- NULL
    out$dose <- dose
    out$variant <- variant
    out$platform <- platform
  }
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
    min_value <- beta_var * 1.2
    max_value <- 0.99
    zero_inf_beta <- function(n, alpha, beta, p_0){
      zero <- as.numeric(runif(n) < p_0)
      rbeta(n, alpha, beta)*abs(zero - 1)
    }
    calculate_alpha_beta <- function(mean){
      mean <- case_when(
        mean < min_value ~ min_value,
        mean > max_value ~ max_value,
        TRUE ~ mean)
      (mean)*((mean)*(1-(mean))/beta_var - 1)
    }
    calculate_beta_beta <- function(mean, alpha){
      mean <- case_when(
        mean < min_value ~ min_value,
        mean > max_value ~ max_value,
        TRUE ~ mean)
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
        ),
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
  message(x)
  if(random_efficacies[[x]]$dose[1] == "First"){
    NULL
  } else {

    plot <- plots[[x]]
    #generate curves for random efficacies
    generate_curve <- function(df){
      user <- as.list(pull(df, value, parameter))
      if (df$dose[1] == "Second") {
        mod <- calc_eff_primary_gen$new(user = user)
      } else {
        mod <- calc_eff_booster_gen$new(user = user)
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
    print(plot + geom_ribbon(data = random_bounds, aes(
      x = t, ymin = lower , ymax = upper, fill = `Protection:`
    ), inherit.aes = FALSE, alpha = 0.1))
  }
})
dev.off()

##Create Sampling Function
random_efficacies <- map_dfr(random_efficacies, ~.x)
#add first dose eff for J&J
random_efficacies <- random_efficacies %>%
  rbind(
    random_efficacies %>% filter(platform == "mRNA" & dose == "First") %>%
      mutate(platform = "Johnson&Johnson", value = 0)
  )
#correct so not decreasing with more vaccinations
random_efficacies <- random_efficacies %>%
  group_by(variant, platform, sample) %>%
  select(variant, platform, sample, parameter, value) %>%
  pivot_wider(names_from = parameter, values_from = value) %>%
  mutate(
    fV_1_d = if_else(
      fV_1_d < pV_1_d, pV_1_d, fV_1_d
    ),
    bV_1_d = if_else(
      bV_1_d < fV_1_d, fV_1_d, bV_1_d
    ),
    fV_1_i = if_else(
      fV_1_i < pV_1_i, pV_1_i, fV_1_i
    ),
    bV_1_i = if_else(
      bV_1_i < fV_1_i, fV_1_i, bV_1_i
    )
  ) %>%
  pivot_longer(!c(variant, platform, sample), names_to = "parameter", values_to = "value") %>%
  left_join(
    random_efficacies %>% select(!value),
    by = c("variant", "platform", "sample", "parameter")
  ) %>%
  ungroup() %>%
  filter(!is.na(dose))


#add Omicron sub unit efficacies
random_efficacies <- random_efficacies %>%
  rbind(
    random_efficacies %>% filter(variant == "Omicron") %>% mutate(variant = "Omicron Sub-Variant")
  ) %>%
  arrange(
    sample, platform, variant, dose, parameter
  ) %>%
  mutate(
    platform = if_else(platform == "Johnson&Johnson", "Single-Dose", platform)
  )

#Scale for breakthrough infections
random_efficacies <- random_efficacies %>%
  filter(str_detect(parameter, "V")) %>%
  mutate(temp_1 = str_remove(parameter, "[di]"),
         temp_2 = str_sub(parameter, -1, -1)) %>%
  select(dose, variant, platform, sample, temp_1, temp_2, value) %>%
  pivot_wider(names_from = temp_2, values_from = value) %>%
  mutate(
    d = (d - i)/(1 - i),
    d = if_else(d < 0, 0, d)
  ) %>%
  pivot_longer(c(d, i), names_to = "temp_2", values_to = "value_new") %>%
  mutate(parameter = paste0(temp_1, temp_2)) %>%
  select(!c(temp_1, temp_2)) %>%
  right_join(
    random_efficacies,
    by = c("dose", "variant", "platform", "sample", "parameter")
  ) %>%
  mutate(
    value = if_else(is.na(value_new), value, value_new)
  ) %>%
  select(!value_new)

#empty environment of everything not relevant
rm(list = setdiff(ls(), c("N_samples", "random_efficacies")))

sample_vaccine_efficacies <- function(n, platforms){
  if("mRNA" %in% platforms){
    booster_platform <- "mRNA"
  } else {
    booster_platform <- "Whole Virus"
  }
  #uniformly sample
  platforms <- sample(platforms, n, replace = TRUE)
  #now for each platform draw a random sample from the df
  output <- map_dfr(seq_along(platforms),
      ~random_efficacies %>%
        filter(platform == platforms[.x], sample == sample(N_samples, 1), dose != "Booster") %>%
        select(dose, variant, parameter, value) %>% #add booster sample
        rbind(
          random_efficacies %>%
            filter(platform == booster_platform, sample == sample(N_samples, 1), dose == "Booster") %>%
            select(dose, variant, parameter, value)
        ) %>%
        mutate(sample = .x)
  )
  #make per variant and format into model compatible formats
  variants <- unique(output$variant)
  output <- map(variants, function(this_variant){
    df <- filter(output, variant == this_variant)
    samples <- map(unique(df$sample), function(this_sample){
      values <- as.list(filter(df, sample == this_sample) %>%
        pull(value, parameter))
      #convert to format
      list(
        dur_V = 1/c(values$w_p, values$w_1, values$w_2),
        vaccine_efficacy_infection = c(values$pV_1_i, values$fV_1_i, values$fV_2_i, values$bV_1_i, values$bV_2_i, values$bV_3_i),
        vaccine_efficacy_disease = c(values$pV_1_d, values$fV_1_d, values$fV_2_d, values$bV_1_d, values$bV_2_d, values$bV_3_d)
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
