set.seed(1000101)

#setup colours
library(scales)
colours <- hue_pal()(3)
names(colours) <- c("Hospitalisation", "Hospitalisation\n(Scaled for Breakthrough)", "Infection")

n_samples <- 50
##Set-up VEs for Vaccine Types
ves_by_type <- read_csv("vaccine_efficacy_groups.csv") %>%
  mutate(dose = if_else(dose == "Partial", "First", "Second")) %>%
  #ensure ve does not increase from Wild to Delta
  pivot_wider(names_from = variant, values_from = efficacy) %>%
  mutate(
    Delta = pmin(Delta, Wild)
  ) %>%
  pivot_longer(
    cols = c("Wild", "Delta"),
    names_to = "variant",
    values_to = "efficacy"
  ) %>%
  mutate(
    vaccine_type = if_else(vaccine_type == "Johnson&Johnson", "Single-Dose", vaccine_type)
  )


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
  "Delta", "Booster", "Hospitalisation", 0.98
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
  ) %>% #average across variants for simplicities sake
  group_by(endpoint)  %>%
  summarise(p_change = mean(p_change))

booster_efficacies <- ves_by_type %>%
  filter(dose == "Second" & vaccine_type %in% c("mRNA", "Whole Virus")) %>%
  left_join(
    changes_booster,
    by = c("endpoint")
  ) %>%
  mutate(
    dose = "Booster",
    efficacy = inv_logit(logit(efficacy) + p_change)
  ) %>%
  select(!p_change)
ves_by_type <- ves_by_type %>%
  rbind(
    booster_efficacies
  ) %>%
  arrange(vaccine_type, variant, endpoint, dose)

#ves_by_type %>%
#  pivot_wider(names_from = variant, values_from = efficacy) %>%
#  select(vaccine_type, endpoint, dose, Wild, Delta, Omicron) %>%
#  filter(vaccine_type == "mRNA")
#load parameter chains
ab_params <- list(
  ni50 = -1.040957,
  ns50 = -1.675559,
  k = 3.182298,
  hl_s = 34.53602,
  hl_l = 573.5029,
  period_s = 75.1927
)

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
  output(ve_i) <- (C1 * vei + C2 * vei_2 + C3 * vei_3)
  w_1 <- user()
  w_2 <- user()
  ved <- user()
  ved_2 <- user()
  ved_3 <- user()
  vei <- user()
  vei_2 <- user()
  vei_3 <- user()
})

simulate_ab <- function(t, initial_ab, h_s, h_l, t_s) {
  pi1 <- -log(2)/h_s
  pi2 <- -log(2)/h_l
  initial_ab * (
    (exp(pi1 * t + pi2 * t_s) + exp(pi1 * t_s + pi2 * t))/
    (exp(pi1 * t_s) + exp(pi2 * t_s))
  )
}

ab_to_ve <- function(ab, n50, k){
  1/(1 + exp(-k * (log10(ab) - n50)))
}

err_lines <- function(l1, l2){
  #sum((l1 - l2)^2/l1)
  sum(((l1 - l2)/l1)^2)
  #scale it so the lower values have more weight
}

calculate_ve <- function(p1, p2, p3) {
  ved <- p1
  ved_2 <- ved * p2
  ved_3 <- ved_2 * p3
  c(ved, ved_2, ved_3)
}

first_doses <- ves_by_type %>%
  filter(dose == "First") %>%
  group_by(vaccine_type, dose, variant) %>%
  group_split() %>%
  map(function(df){
    parameter_infection <- df %>% filter(endpoint == "Infection") %>% pull(efficacy)
    parameter_hospitalisation <- df %>% filter(endpoint == "Hospitalisation") %>% pull(efficacy)
    dose <- unique(df$dose)
    variant <- unique(df$variant)
    platform <- unique(df$vaccine_type)
    tibble(
      parameter = c("pV_i", "pV_d"),
      value = c(parameter_infection, (parameter_hospitalisation - parameter_infection)/(1 - parameter_infection)),
      dose = dose,
      variant = variant,
      platform = platform
    )
  })
plot_legend <- NULL
other_doses <- ves_by_type %>%
  filter(dose != "First") %>%
  group_by(vaccine_type, dose, variant) %>%
  group_split() %>%
  map(function(df){
    parameter_infection <- df %>% filter(endpoint == "Infection") %>% pull(efficacy)
    parameter_hospitalisation <- df %>% filter(endpoint == "Hospitalisation") %>% pull(efficacy)

    dose <- unique(df$dose)
    variant <- unique(df$variant)
    platform <- unique(df$vaccine_type)

    message(paste0(dose, "; ", variant, "; ", platform))

    t_plot <- seq(0, simulate_time, length.out = 100)

    calc_eff <- calc_eff_gen$new(user = list(
      ved = parameter_hospitalisation,
      vei = parameter_infection,
      ved_2 = 1,
      ved_3 = 1,
      vei_2 = 1,
      vei_3 = 1,
      w_1 = 0,
      w_2 = 0
    ))

    #need to calculate initial dose level
    #assume ve is 30 days after dose
    t_measure <- 30
    #just assume the initia AB is
    err_func <- function(initial_ab) {
      ab_t_measure <- simulate_ab(t_measure, initial_ab, ab_params$hl_s, ab_params$hl_l, ab_params$period_s)
      ve_i <- ab_to_ve(ab_t_measure, ab_params$ni50, ab_params$k)
      ve_d <- ab_to_ve(ab_t_measure, ab_params$ns50, ab_params$k)
      err_lines(parameter_infection, ve_i) +
        err_lines(parameter_hospitalisation, ve_d)
    }
    res <- optimize(err_func, interval = c(0, 10), maximum = FALSE)
    initial_ab <- res$minimum
    abs <- simulate_ab(0:simulate_time, initial_ab,  ab_params$hl_s,  ab_params$hl_l,  ab_params$period_s)
    ve_i <- ab_to_ve(abs, ab_params$ni50, ab_params$k)
    ve_d <- ab_to_ve(abs, ab_params$ns50, ab_params$k)
    #plot for diagnostics
    initial_ab_plot <- tibble(
      t = 0:(2*t_measure),
      Hospitalisation = ve_d[1:(2*t_measure + 1)],
      Infection = ve_i[1:(2*t_measure + 1)]
    ) %>% 
      pivot_longer(cols = -t, names_to = "Endpoint", values_to = "Efficacy") %>%
      mutate(
        Endpoint = factor(Endpoint, levels = c("Infection", "Hospitalisation", "Hospitalisation\n(Scaled for Breakthrough)"))
      ) %>%
      ggplot(aes(x = t, y = Efficacy, colour = Endpoint)) +
      geom_line(show.legend = FALSE) +
      geom_point(data = tibble(
        t = t_measure,
        Efficacy = c(parameter_infection, parameter_hospitalisation),
        Endpoint = c("Infection", "Hospitalisation")
      ) %>%
        mutate(
          Endpoint = factor(Endpoint, levels = c("Infection", "Hospitalisation", "Hospitalisation\n(Scaled for Breakthrough)"))
        ), shape = "x", size = 5, show.legend = FALSE) +
      ggpubr::theme_pubclean() +
      scale_colour_manual(values = colours) +
      labs(
        title = "Fit of initial Immune Response level, X shows the input efficacies",
        y = "Vaccine Efficacy", x = "Days Since Dose"
      )

    #scale for break through infection
    ve_d <- (ve_d - ve_i)/(1 - ve_i)

    err_func <- function(pars) {
      veds <- calculate_ve(pars[3], pars[4], pars[5])
      veis <- calculate_ve(pars[6], pars[7], pars[8])
      calc_eff$set_user(
        user = list(
          w_1 = pars[1],
          w_2 = pars[2],
          ved = veds[1],
          ved_2 = veds[2],
          ved_3 = veds[3],
          vei = veis[1],
          vei_2 = veis[2],
          vei_3 = veis[3]
        )
      )
      mod_value <- calc_eff$run(t = seq(0, simulate_time))
      log(sqrt(err_lines(ve_d, mod_value[, "ve_d"]) + err_lines(ve_i, mod_value[, "ve_i"])))
    }
    lower = list(
      w_1 = 1/(3*365),
      w_2 = 1/(3*365),
      ved_1 = 0,
      ved_2 = 0,
      ved_3 = 0,
      vei_1 = 0,
      vei_2 = 0,
      vei_3 = 0
    )
    upper = list(
      w_1 = 1/30,
      w_2 = 1/30,
      ved_1 = 1,
      ved_2 = 1,
      ved_3 = 1,
      vei_1 = 1,
      vei_2 = 1,
      vei_3 = 1
    )
    par =  list(
      w_1 = 1/365,
      w_2 = 1/365,
      ved_1 = 0.5,
      ved_2 = 0.5,
      ved_3 = 0.5,
      vei_1 = 0.5,
      vei_2 = 0.5,
      vei_3 = 0.5
    )
    res <- dfoptim::nmkb(unlist(par), fn = err_func, lower = unlist(lower), upper = unlist(upper), control = list(maxfeval = 5000))
    if(res$convergence != 0){
      stop(res$message)
    }

    #add randomness (should do it in fitting really)
    out <- map(seq_len(n_samples), function(i){
      pars <- res$par
      if(i != 1){
        pars[1:2] <- 1/pars[1:2]
        pars <- pars + runif(length(pars), -0.05, 0.05)
        pars[1:2] <- 1/pars[1:2]
        pars[pars < 0] <- 0
        pars[pars > 1] <- 1
      }
      veds <- calculate_ve(pars[3], pars[4], pars[5])
      veis <- calculate_ve(pars[6], pars[7], pars[8])
      tibble(
        value = c(pars[1:2], veds, veis),
        parameter = names(par),
        iteration = i
      )
    })

    #plot
    t_plot <- 0:simulate_time
    p <- map_dfr(out, function(pars){
      calc_eff$set_user(
        user = list(
          w_1 = pars$value[1],
          w_2 = pars$value[2],
          ved = pars$value[3],
          ved_2 = pars$value[4],
          ved_3 = pars$value[5],
          vei = pars$value[6],
          vei_2 = pars$value[7],
          vei_3 = pars$value[8]
        )
      )
      mod_value <- calc_eff$run(t = t_plot)
      tibble(
        t = rep(t_plot, 2),
        value = c(mod_value[, "ve_d"], mod_value[, "ve_i"]),
        Endpoint = c(rep("Hospitalisation\n(Scaled for Breakthrough)", length(t_plot)), rep("Infection", length(t_plot))),
        iteration = pars$iteration[1]
      )
    }) %>%
      mutate(
        model = "Booster Model"
      ) %>%
      rbind(
        tibble(
          t = t_plot,
          value = ve_d,
          Endpoint = "Hospitalisation\n(Scaled for Breakthrough)",
          iteration = 0,
          model = "AB Process"
        )
      ) %>%
      rbind(
        tibble(
          t = t_plot,
          value = ve_i,
          Endpoint = "Infection",
          iteration = 0,
          model = "AB Process"
        )
      ) %>%
      mutate(
        group_par = paste0(Endpoint, "_", iteration, "_", model),
        alpha = ifelse(model == "AB Process", 1, 0.25)
      ) %>%
      mutate(
        Endpoint = factor(Endpoint, levels = c("Infection", "Hospitalisation", "Hospitalisation\n(Scaled for Breakthrough)"))
      ) %>%
      ggplot(aes(x = t, y = value, color = Endpoint, linetype = model, group = group_par, alpha = alpha)) +
        geom_line() +
      labs(y = "Vaccine Efficacy", x = "Days Since Dose", title = paste0("Dose: ", dose, ", Variant: ", variant, ", Type: ", df$vaccine_type[1]), linetype = "Model", colour = "Endpoint") +
      ggpubr::theme_pubclean() + scale_alpha(guide = 'none') + 
      scale_colour_manual(values = colours, drop = FALSE) +
      ylim(c(0, 1))

    if(is.null(plot_legend)){
      plot_legend <<- get_legend(p, position = "top") #only need one, all identical
      dev.off() #don't know why it does this
    }
    p <- list(
      p = p, ab = initial_ab_plot
    )

    out <- map_dfr(out, ~.x)

    out$dose <- dose
    out$variant <- variant
    out$platform <- platform

    list(
      out = out,
      plot = p
    )
})
#split into plots and data
plots <- map(other_doses, ~.x$plot)
other_doses <- map(other_doses, ~.x$out)

##Calibration plots

platform <- map_chr(other_doses, ~.x$platform[1])
platforms <- unique(platform)
variant <- map_chr(other_doses, ~.x$variant[1])
variants <- c("Wild", "Delta")
dose <- map_chr(other_doses, ~.x$dose[1])
doses <- unique(dose)

p <- map(platforms, function(plat){
  dose_list <- map(doses, function(dos){
    var_list <- map(variants, function(vari){
      index <- detect_index(other_doses, ~.x$platform[1] == plat & .x$variant[1] == vari & .x$dose[1] == dos)
      if(index > 0){
        ggarrange(plots[[index]]$p + theme(legend.position = "none"), plots[[index]]$ab, ncol = 1, heights = c(0.6, 0.2))
      }
    })
    if(every(var_list, is.null)){
      return(NULL)
    } else {
      return(ggarrange(plotlist = var_list, nrow = 1))
    }
  }) %>%
    compact()
  list(
    height = (length(dose_list) + 0.1)/(0.1 + 2),
    plot = ggarrange(plotlist = c(list(legend), dose_list), heights = c(0.1, rep(1, length(dose_list))), ncol = 1, legend.grob = plot_legend)
  )
})

dir.create("plots", showWarnings = FALSE)
iwalk(p, \(x, idx) {
  ggsave(paste0("plots/", idx, ".png"), x$plot, width = 15, height = 20 * x$height, bg = "white")
})

zip("plots.zip", file.path("plots", list.files("plots")))

unlink("plots", recursive = TRUE)

rm(p, plots)

first_doses <- map_dfr(first_doses, ~.x)
other_doses <- map_dfr(other_doses, ~.x)

#add first dose eff for J&J
first_doses <- first_doses %>%
  rbind(
    first_doses %>% filter(platform == "mRNA") %>%
      mutate(platform = "Single-Dose", value = 0)
  )
#sense check so that first < second < booster
silent <- other_doses %>%
  rbind(first_doses %>%
  mutate(iteration = 0)) %>%
  group_by(platform, variant) %>%
  group_split() %>%
  map(function(df){
    map(c("i", "d"), function(type){
      first <- df %>% filter(iteration == 0, parameter == paste0("pV_", type)) %>% pull(value)
      second <- df %>% filter(iteration != 0, parameter == paste0("ve", type, "_", 1), dose == "Second") %>% pull(value)
      if(any(second < first)){
        warning(paste0("Second dose is less than first dose in ", type, " for ", df$platform[1], " ", df$variant[1]))
      }
      if("Booster" %in% df$dose){
        booster <- df %>% filter(iteration != 0, parameter == paste0("ve", type, "_", 1), dose == "Booster") %>% pull(value)
        if(any(unlist(map(booster, ~.x <= second)))){
          warning(paste0("Booster dose is less than second dose in ", type, " for ", df$platform[1], " ", df$variant[1]))
        }
      }
    })
  })

other_doses <- other_doses %>%
  mutate(
    parameter = case_when(
      parameter == "w_1" & dose == "Second" ~ "fw_1",
      parameter == "w_2" & dose == "Second"~ "fw_2",
      parameter == "ved_1" & dose == "Second" ~ "fV_d_1",
      parameter == "ved_2" & dose == "Second" ~ "fV_d_2",
      parameter == "ved_3" & dose == "Second" ~ "fV_d_3",
      parameter == "vei_1" & dose == "Second" ~ "fV_i_1",
      parameter == "vei_2" & dose == "Second" ~ "fV_i_2",
      parameter == "vei_3" & dose == "Second" ~ "fV_i_3",
      parameter == "w_1" & dose == "Booster" ~ "bw_1",
      parameter == "w_2" & dose == "Booster"~ "bw_2",
      parameter == "ved_1" & dose == "Booster" ~ "bV_d_1",
      parameter == "ved_2" & dose == "Booster" ~ "bV_d_2",
      parameter == "ved_3" & dose == "Booster" ~ "bV_d_3",
      parameter == "vei_1" & dose == "Booster" ~ "bV_i_1",
      parameter == "vei_2" & dose == "Booster" ~ "bV_i_2",
      parameter == "vei_3" & dose == "Booster" ~ "bV_i_3",
      TRUE ~ parameter
    )
  )
#empty environment of everything not relevant
rm(list = setdiff(ls(), c("first_doses", "other_doses", "n_samples")))

sample_vaccine_efficacies <- function(n, platforms){
  if("mRNA" %in% platforms){
    booster_platform <- "mRNA"
  } else {
    booster_platform <- "Whole Virus"
  }
  #uniformly sample
  platforms <- sample(platforms, n, replace = TRUE)
  its_primary <- sample(n_samples, n, replace = TRUE)
  its_booster <- sample(n_samples, n, replace = TRUE)
  #now for each platform draw a random sample from the df
  output <- map_dfr(seq_len(n),
      ~first_doses %>%
        filter(platform == platforms[.x]) %>%
        rbind(
          other_doses %>%
            filter(platform == platforms[.x], iteration == its_primary[.x], dose == "Second") %>%
            rbind(
              other_doses %>%
                filter(platform == booster_platform, iteration == its_booster[.x], dose == "Booster")
            ) %>%
            select(!iteration)
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
        dur_V = 1/c(values$fw_1, values$fw_2, values$bw_1, values$bw_2),
        vaccine_efficacy_infection = c(values$pV_i, values$fV_i_1, values$fV_i_2, values$fV_i_3, values$bV_i_1, values$bV_i_2, values$bV_i_3),
        vaccine_efficacy_disease = c(values$pV_d, values$fV_d_1, values$fV_d_2, values$fV_d_3, values$bV_d_1, values$bV_d_2, values$bV_d_3)
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
