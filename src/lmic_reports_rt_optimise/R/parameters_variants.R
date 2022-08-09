# Produces functions to sample from the distributions on variant related
# parameters, excluding vaccine efficacies.


#immune escape against Wild type
sample_variant_immune_escape <- function(n, variants){
  parameters <- tribble(
    ~variant, ~p1, ~p2,
    "Delta", 1.014, 2,
    "Omicron", 2.54, 2,
    "Omicron Sub-Variant", 1.014, 2 #ISSUE:: GET A REAL NUMBER
  ) %>%
    filter(variant %in% variants)
  out <- map(seq_len(n), function(it){
    #generate values
    immune_escapes <- map(transpose(parameters), ~rbeta(1, .x$p1, .x$p2))
    names(immune_escapes) <- parameters$variant
    #convert to escape against wild type, assume the order in parameters is the correct order
    for(i in seq_along(immune_escapes)[-1]){
      immune_escapes[[i]] <- (1-immune_escapes[[i-1]]) * immune_escapes[[i]] +
        immune_escapes[[i-1]]
    }
    immune_escapes
  })
  output <- map(names(out[[1]]), function(x){map_dbl(out, ~.x[[x]])})
  names(output) <- names(out[[1]])
  output
}

#against Wild type
sample_variant_prob_hosp <- function(n, variants){
  parameters <- tribble(
    ~variant, ~p1, ~p2,
    "Delta", log(1.45), 0.15,
    "Omicron", log(0.59), 0.08,
    "Omicron Sub-Variant", log(1), 0.08 #ISSUE:: GET A REAL NUMBER
  ) %>%
    filter(variant %in% variants)
  out <- map(seq_len(n), function(it){
    #generate values
    prob_multipliers <- map(transpose(parameters), ~rlnorm(1, .x$p1, .x$p2))
    names(prob_multipliers) <- parameters$variant
    #convert against wild type, assume the order in parameters is the correct order
    for(i in seq_along(prob_multipliers)[-1]){
      prob_multipliers[[i]] <- prob_multipliers[[i-1]] * prob_multipliers[[i]]
    }
    prob_multipliers
  })
  output <- map(names(out[[1]]), function(x){map_dbl(out, ~.x[[x]])})
  names(output) <- names(out[[1]])
  output
}

sample_variant_prob_severe <- function(n, variants){
  parameters <- tribble(
    ~variant, ~p1, ~p2,
    "Delta", log(1), 0.08,
    "Omicron", log(0.34), 0.45,
    "Omicron Sub-Variant", log(1), 0.08 #ISSUE:: GET A REAL NUMBER
  ) %>%
    filter(variant %in% variants)
  out <- map(seq_len(n), function(it){
    #generate values
    prob_multipliers <- map(transpose(parameters), ~rlnorm(1, .x$p1, .x$p2))
    names(prob_multipliers) <- parameters$variant
    #convert against wild type, assume the order in parameters is the correct order
    for(i in seq_along(prob_multipliers)[-1]){
      prob_multipliers[[i]] <- prob_multipliers[[i-1]] * prob_multipliers[[i]]
    }
    prob_multipliers
  })
  output <- map(names(out[[1]]), function(x){map_dbl(out, ~.x[[x]])})
  names(output) <- names(out[[1]])
  output
}

estimate_generation_time <- function(variant_timings, start_date){
  #update generation times
  #wild
  latent <- nimue:::default_durations()$dur_E
  case <- nimue:::default_durations()$dur_ICase
  gen_time <- 6.3 #https://doi.org/10.1016/S1473-3099(20)30287-5
  p_mild <- 0.26 + 0.65
  mild <- (gen_time - latent - ((1 - p_mild) * case / 2)) * 2 / p_mild
  #delta
  gen_time_delta <- 4.7 #https://doi.org/10.1016/S1473-3099(22)00001-9
  p_mild_delta <- 1-(369 + 129)/(6681 + 2001) #not aviabile so taken from https://doi.org/10.1016/S1473-3099(21)00475-8
  #just scale down all values
  #this is easier than extracting their estimates of latent etc
  p_gen_change_delta <- gen_time_delta/gen_time
  latent_delta <- p_gen_change_delta * latent
  mild_delta <- p_gen_change_delta * mild * p_mild / p_mild_delta
  case_delta <- p_gen_change_delta * case * (1 - p_mild) / (1 - p_mild_delta)
  #Omicron
  gen_time_omicron <- 3.6 #https://www.gov.uk/government/publications/spi-m-o-consensus-statement-on-covid-19-6-january-2022/spi-m-o-consensus-statement-on-covid-19-6-january-2022#fn:1
  p_mild_omicron <- p_mild_delta#just keep same as delta for now
  p_gen_change_omicron <- gen_time_omicron/gen_time_delta
  latent_omicron <- p_gen_change_omicron * latent_delta
  mild_omicron <- p_gen_change_omicron * mild_delta * p_mild_delta / p_mild_omicron
  case_omicron <- p_gen_change_omicron * case_delta * (1 - p_mild_delta) / (1 - p_mild_omicron)
  #correct format
  latent <- list(Wild = latent, Delta = latent_delta, Omicron = latent_omicron)
  mild <- list(Wild = mild, Delta = mild_delta, Omicron = mild_omicron)
  case <- list(Wild = case, Delta = case_delta, Omicron = case_omicron)
  #ensure the same variants
  variant_timings <- variant_timings %>%
    filter(variant %in% names(latent))
  latent <- latent[c("Wild", variant_timings$variant)]
  mild <- mild[c("Wild", variant_timings$variant)]
  case <- case[c("Wild", variant_timings$variant)]
  #adjust over time
  latent <- variant_par_changes_over_time(variant_timings, latent, start_date)
  mild <- variant_par_changes_over_time(variant_timings, mild, start_date)
  case <- variant_par_changes_over_time(variant_timings, case, start_date)
  #name correctly and return
  list(
    dur_E = latent$var,
    tt_dur_E = latent$tt,
    dur_IMild = mild$var,
    tt_dur_IMild = mild$tt,
    dur_ICase = case$var,
    tt_dur_ICase = case$tt
  )
}

estimate_healthcare_durations <- function(variant_timings, start_date){
  #wild
  mv_survive <- nimue:::default_durations()$dur_get_mv_survive
  mv_die <- nimue:::default_durations()$dur_get_mv_die
  ox_survive <- nimue:::default_durations()$dur_get_ox_survive
  ox_die <- nimue:::default_durations()$dur_get_ox_survive
  #delta
  mv_survive_delta <- mv_survive
  mv_die_delta <- mv_die
  ox_survive_delta <- 5 #https://www.medrxiv.org/content/10.1101/2022.01.11.22269045v1
  ox_die_delta <- ox_die
  #omicron
  mv_survive_omicron <- 8 #INARC report Juky 2022
  mv_die_omicron <- 5
  ox_survive_omicron <- 5
  ox_die_omicron <- 6
  #correct format
  mv_survive <- list(Wild = mv_survive, Delta = mv_survive_delta, Omicron = mv_survive_omicron)
  mv_die <- list(Wild = mv_die, Delta = mv_die_delta, Omicron = mv_die_omicron)
  ox_survive <- list(Wild = ox_survive, Delta = ox_survive_delta, Omicron = ox_survive_omicron)
  ox_die <- list(Wild = ox_die, Delta = ox_die_delta, Omicron = ox_die_omicron)
  #ensure the same variants
  variant_timings <- variant_timings %>%
    filter(variant %in% names(mv_survive))
  mv_survive <- mv_survive[c("Wild", variant_timings$variant)]
  mv_die <- mv_die[c("Wild", variant_timings$variant)]
  ox_survive <- ox_survive[c("Wild", variant_timings$variant)]
  ox_die <- ox_die[c("Wild", variant_timings$variant)]
  #adjust over time
  mv_survive <- variant_par_changes_over_time(variant_timings, mv_survive, start_date)
  mv_die <- variant_par_changes_over_time(variant_timings, mv_die, start_date)
  ox_survive <- variant_par_changes_over_time(variant_timings, ox_survive, start_date)
  ox_die <- variant_par_changes_over_time(variant_timings, ox_die, start_date)
  #name correctly and return
  list(
    dur_get_mv_survive = mv_survive$var,
    tt_dur_get_mv_survive = mv_survive$tt,
    dur_get_mv_die = mv_die$var,
    tt_dur_get_mv_die = mv_die$tt,
    dur_get_ox_survive = ox_survive$var,
    tt_dur_get_ox_survive = ox_survive$tt,
    dur_get_ox_die = ox_die$var,
    tt_dur_get_ox_die = ox_die$tt
  )
}
