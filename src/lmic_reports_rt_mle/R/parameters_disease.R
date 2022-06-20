# Produces functions to sample from the distributions on COVID disease related
# parameters, non-variant dependant.


#method: sample overall IFR then adjust constituent values to ensure its consistent
#also sample no treatment parameters ensuring they are larger than the treatment parametesr
sample_ifr <- function(n, iso3c){
  scale_func <- function(scale, central, min, max){
    if(scale < 0.5){
      min + (central - min)*scale/0.5
    } else {
      central + (max - central)*(scale-0.5)/0.5
    }
  }
  ifr_report <- tribble(
    ~central, ~low, ~high,
    0.00, 0.00, 0.03,
    0.01, 0.00, 0.06,
    0.01, 0.00, 0.11,
    0.02, 0.00, 0.18,
    0.03, 0.00, 0.30,
    0.04, 0.00, 0.46,
    0.06, 0.01, 0.71,
    0.10, 0.01, 1.03,
    0.16, 0.02, 1.47,
    0.24, 0.03, 2.03,
    0.38, 0.05, 2.74,
    0.60, 0.10, 3.64,
    0.94, 0.18, 4.79,
    1.47, 0.35, 6.27,
    2.31, 0.65, 8.21,
    3.61, 1.21, 10.8,
    #caculate these based on the countries population
    weighted.mean(c(5.66, 8.86, 17.37), squire::get_elderly_population(iso3c = iso3c)$n),
    weighted.mean(c(2.23, 4.06, 9.7), squire::get_elderly_population(iso3c = iso3c)$n),
    weighted.mean(c(14.37, 19.36, 31.12), squire::get_elderly_population(iso3c = iso3c)$n)
  )
  samples <- map(seq_len(n), function(it){
    #get the IFR to recreate
    scale <- rbeta(1, 2, 2)

    target_ifr <-
      scale_func(scale, ifr_report$central, ifr_report$low, ifr_report$high)/100

    pars <- map(seq_along(target_ifr), function(a){
      #try prioritise keeping prob hosp, then prob severe and then the ratio of deaths the same
      N <- 50
      parameters <- expand_grid(
        prob_hosp = seq(nimue:::probs$prob_hosp[a], 1, length.out = N),
        prob_severe = seq(nimue:::probs$prob_severe[a], 1, length.out = N),
        prop_icu = seq(0.8, 1, length.out = N)
      )
      fitting <- TRUE
      index <- 0
      while(fitting){
        index <- index + 1
        p_s_d <- target_ifr[a]/parameters$prob_hosp[index] * 1/(parameters$prob_severe[index] + (1 - parameters$prob_severe[index]) * (1 - parameters$prop_icu[index])/parameters$prop_icu[index])
        p_ns_d <- p_s_d*(1-parameters$prop_icu[index])/parameters$prop_icu[index]
        if(p_s_d >= 0 & p_s_d <= 1 & p_ns_d >= 0 & p_ns_d <= 1){
          fitting <- FALSE
        }
      }
      #add non treatment parameters
      var <- 0.001
      means <- c(nimue:::probs$prob_severe_death_no_treatment,
                 nimue:::probs$prob_non_severe_death_no_treatment)
      alphas <- means*(means*(1-means)/var - 1)
      betas <- (alphas - means*alphas)/means
      prob_severe_death_no_treatment <- rbeta(1, alphas[1], betas[1])
      prob_non_severe_death_no_treatment <- rbeta(1, alphas[2], betas[2])

      list(prob_hosp = parameters$prob_hosp[index],
           prob_severe = parameters$prob_severe[index],
           prob_non_severe_death_treatment = p_ns_d,
           prob_severe_death_treatment = p_s_d,
           prob_severe_death_no_treatment = prob_severe_death_no_treatment,
           prob_non_severe_death_no_treatment = prob_non_severe_death_no_treatment
      )
    })
    out <- map(
      names(pars[[1]]),
      function(x){map_dbl(pars, ~.x[[x]])}
    )
    names(out) <- names(pars[[1]])
    out
  })
  #rearrange for easier usage
  out <- map(
    names(samples[[1]]),
    function(x){map(samples, ~.x[[x]])}
  )
  names(out) <- names(samples[[1]])
  out
}

sample_duration_natural_immunity <- function(n){
  rgamma(n, 20, 4/73)
}
