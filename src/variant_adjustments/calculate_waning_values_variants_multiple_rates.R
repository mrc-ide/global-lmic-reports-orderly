library(tidyverse)
library(odin)
library(ggpubr)
set.seed(1001)
initial_efficacies <- tribble(
  ~Dose, ~Variant, ~Disease, ~Infection,
  "Partial", "Wild", 0.8, 0.6,
  "Full", "Wild", 0.98, 0.80,
  "Partial", "Delta", 0.75, 0.224,
  "Full", "Delta", 0.94, 0.646,
  "Partial", "Omicron", 0.4192211, 0,
  "Full", "Omicron", 0.5254237, 0.1
)
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

values <- pmap(initial_efficacies, function(...){
  parameters <- tibble(...)
  simulate_time <- 2*365
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
  ab_d <- map(parameters$Disease, ~c(10^((-log((1/.x) - 1)/k) + log10(n_50_d)),
          rep(NA, simulate_time)))
  ab_i <- map(parameters$Infection, ~c(10^((-log((1/.x) - 1)/k) + log10(n_50_i)),
            rep(NA, simulate_time)))
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
    ab_d <- map(ab_d, function(x){
      x[t+1] <- x[t]*exp(-decay)
      return(x)
    })
    ab_i <- map(ab_i, function(x){
      x[t+1] <- x[t]*exp(-decay)
      return(x)
    })
  }
  #convert to efficacy curve
  ve_d <- map(ab_d, ~1/(1 + exp(-k*(log10(.x) - log10(n_50_d)))))
  ve_i <- map(ab_i, ~1/(1 + exp(-k*(log10(.x) - log10(n_50_i)))))
  #fit parameters
  if(parameters$Dose == "Partial"){
    calc_eff_partial <- calc_eff_partial_gen$new(user = list(
      pV_1_d = parameters$Disease,
      pV_2_d = 0,
      pV_1_i = parameters$Infection,
      pV_2_i = 0,
      w_p = 0
    ))

    err_func <- function(pars) {
      pars <- c(pars, 0) #add an empty term of the omicron I waned (=0)
      map_dbl(seq_along(parameters$Disease), function(x){
        calc_eff_partial$set_user(
          user = list(
            pV_1_d =  parameters$Disease[x],
            pV_1_i = parameters$Infection[x],
            w_p = pars[1 + 3*(x-1)],
            pV_2_d = pars[2 + 3*(x-1)],
            pV_2_i = pars[3 + 3*(x-1)]
          )
        )
        mod_value <- calc_eff_partial$run(t = c(0, seq_len(simulate_time)))
        sum((ve_d[[x]] - mod_value[, "ve_d"])^2) + sum((ve_i[[x]] - mod_value[, "ve_i"])^2)
      }) %>% sum %>% sqrt %>% log

    }
    if(parameters$Variant == "Delta"){
      lower = list(w_p = 1/(3*365), pV_2_d = 0.001, pV_2_i = 0.001)
    } else {
      lower = list(w_p = 1/(3*365), pV_2_d = 0.001, pV_2_i = 0)
    }
    upper = list(w_p = 1/30, pV_2_d = parameters$Disease, pV_2_i = parameters$Infection)
    par = list(w_p = 1/365, pV_2_d = parameters$Disease/2, pV_2_i = parameters$Infection/2)
    if(parameters$Variant == "Omicron"){
      lower$pV_2_i <- NULL
      upper$pV_2_i <- NULL
      par$pV_2_i <- NULL
    }
  } else {
    calc_eff_full <- calc_eff_full_gen$new(user = list(
      fV_1_d = parameters$Disease[1],
      fV_2_d = 0,
      fV_3_d = 0,
      fV_1_i = parameters$Infection[1],
      fV_2_i = 0,
      fV_3_i = 0,
      w_1 = 0,
      w_2 = 0
    ))
    err_func <- function(pars) {
      map_dbl(seq_along(parameters$Disease), function(x){
        calc_eff_full$set_user(
          user = list(
            w_1 = pars[1],
            w_2 = pars[2],
            fV_1_d =  parameters$Disease[x],
            fV_1_i = parameters$Infection[x],
            fV_2_d = pars[1 + 2 + 4*(x-1)],
            fV_3_d = pars[2 + 2 + 4*(x-1)],
            fV_2_i = pars[3 + 2 + 4*(x-1)],
            fV_3_i = pars[4 + 2 + 4*(x-1)]
          )
        )
        mod_value <- calc_eff_full$run(t = c(0, seq_len(simulate_time)))
        sum((ve_d[[x]] - mod_value[, "ve_d"])^2) + sum((ve_i[[x]] - mod_value[, "ve_i"])^2)
      }) %>% sum %>% sqrt %>% log
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
      fV_2_d = parameters$Disease[1],
      fV_3_d = parameters$Disease[1],
      fV_2_i = parameters$Infection[1],
      fV_3_i = parameters$Infection[1]
    )
    par =  list(
      w_1 = 1/365,
      w_2 = 1/365,
      fV_2_d = parameters$Disease[1]/2,
      fV_3_d = parameters$Disease[1]/2,
      fV_2_i = parameters$Infection[1]/2,
      fV_3_i = parameters$Infection[1]/2
    )
  }
  res <- optim(par, fn = err_func, method = "L-BFGS-B", lower = lower, upper = upper)
  if(res$convergence != 0){
    print(paste0("Dose: ", parameters$Dose, ", Variant: ", parameters$Variant))
    stop(res$message)
  }
  if(parameters$Variant == "Delta" & parameters$Dose == "Partial" & res$par[3] == 0.001){
    res$par[3] <- 0
  }
  #make plot
  p <- map(seq_along(parameters$Disease), function(x){
    if(parameters$Dose == "Partial"){
      calc_eff_partial$set_user(
        user = list(
          w_p = res$par[1],
          pV_1_d =  parameters$Disease[x],
          pV_1_i = parameters$Infection[x],
          pV_2_d = res$par[1 + x*2 - 1],
          pV_2_i = res$par[1 + x*2 - 1]
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
          fV_1_d =  parameters$Disease[x],
          fV_1_i = parameters$Infection[x],
          fV_2_d = res$par[1 + 2 + 4*(x-1)],
          fV_3_d = res$par[2 + 2 + 4*(x-1)],
          fV_2_i = res$par[3 + 2 + 4*(x-1)],
          fV_3_i = res$par[4 + 2 + 4*(x-1)]
        )
      )
      mod_value <- calc_eff_full$run(t = c(0, seq_len(simulate_time)))
      ve_f_d <- mod_value[, "ve_d"]
      ve_f_i <- mod_value[, "ve_i"]
    }
    ggplot(
      tibble(
        t = rep(c(0, seq_len(simulate_time)), 4),
        `Protection:` = c(rep("Disease", (simulate_time + 1)*2), rep("Infection", (simulate_time + 1)*2)),
        `Version:` = rep(c(rep("AB Process", simulate_time + 1), rep("Weibull", simulate_time + 1)), 2),
        `Vaccine Efficacy` = c(ve_d[[x]], ve_f_d, ve_i[[x]], ve_f_i)
      )
      , aes(x = t, y = `Vaccine Efficacy`, colour = `Protection:`, linetype = `Version:`)
    ) +
      geom_line() +
      labs(y = "Vaccine Efficacy", x = "Days Since Dose", title = paste0("Dose: ", parameters$Dose, ", Variant: ", parameters$Variant)) +
      ggpubr::theme_pubclean()
  })
  out <- as.data.frame(res$par)
  out$Dose <- parameters$Dose
  out$Variant <- parameters$Variant
  #return parameters
  return(list(out, p))
})

ggsave("meta/variant_fitted_value_plot.pdf",
       ggarrange(plotlist = map(values, ~.x[[2]][[1]]), common.legend = TRUE, ncol = 2, nrow = 3)
)
map(values, ~.x[[1]])
