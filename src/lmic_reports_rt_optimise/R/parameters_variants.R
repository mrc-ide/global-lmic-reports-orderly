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
