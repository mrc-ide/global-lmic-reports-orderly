curve <- function(t, centre, growth){
  #logistic bound to 5% and 95% so we don't get infinite values
  1/(1 + exp(-growth*(t - centre)))
}

likelihood <- function(df, centres, growths, background_chance){
  #only keep dates where there are some sequencing
  tts <- unique(df$t)
  variants <- sort(unique(df$variant))[-1]
  map_dbl(tts, function(tt){
    this_data <- filter(df, t == tt)
    probs <- curve(tt, centres, growths)
    #assume each subsequent variant takes priority of the previous
    #order based on when each variant came into effect
    probs <- probs[order(centres)]

    # ##CARRY THIS ON AT SOME POINT, JUST KEEP GROWTH THE SAME FOR NOW
    # adjusted_probs <- rep(NA, length(probs))
    # for(i in rev(seq_along(probs))){
    #   #if check if there are any parameters ahead of this
    #   if(i == length(adjusted_probs)){
    #     #then the prob is not adjusted
    #     adjusted_probs[i] <- probs[i]
    #   } else {
    #     #subtract the previous probability (set to 0 if it would be negative)
    #     max_prob <- max(probs[(i + 1):length(probs)])
    #     adjusted_probs[i] <- if_else(
    #       probs[i] < max_prob,
    #       0,
    #       probs[i] - max_prob
    #     )
    #   }
    # }

    #TEMP CODE, assumes  that growths are all equal:
    adjusted_probs <- probs - lead(probs, default = 0)


    p_wild <- 1 - sum(adjusted_probs)
    names(p_wild) <- "Wild"
    adjusted_probs <- c(adjusted_probs, p_wild)
    adjusted_probs <- adjusted_probs[map_int(this_data$variant, ~which(names(adjusted_probs) == .x))]
    #assume a minimal background chance of seeing a variant independent of its
    #status i.e. uniformly selected from variants
    adjusted_probs <- background_chance * 1/length(adjusted_probs) +
      (1 - background_chance) * adjusted_probs

    dmultinom(this_data$count, size = sum(this_data$count), prob = adjusted_probs,
              log = TRUE)
  }) %>%
    sum()
}

#fit world wide with drjacoby, will do for now
worldwide <- variants_df %>%
  ungroup() %>%
  mutate(t = as.numeric(date - min(date))) %>%
  group_by(t, variant) %>%
  summarise(count = sum(count), .groups = "drop")%>%
  group_by(t) %>%
  filter(sum(count) > 0) %>%
  ungroup() %>%
  arrange(t, variant)
max_centre <- max(worldwide$t)
variants <- unique(worldwide$variant)[-1]
df_params <- tibble(
  name = paste0(variants, "_centre"),
  min = 0,
  max = max_centre
) %>%
  # rbind(
  #   tibble(
  #     name = paste0(variants, "_growth"),
  #     min = 0.01,
  #     max = 1
  #   )
  # ) %>%
  add_row(
    name = "growth",
    min = 0.01,
    max = 1
  ) %>%
  add_row(
    name = "p_background",
    min = 0,
    max = 0.5
  )
dr_j_likelihood <- function(params, data, misc){
  centres <- params[stringr::str_subset(names(params), "_centre")]
  names(centres) <- stringr::str_remove(names(centres), "_centre")
  # growths <- params[stringr::str_subset(names(params), "_growth")]
  # names(growths) <- stringr::str_remove(names(growths), "_growth")
  growths <- rep(params["growth"], length(centres))
  names(growths) <- names(centres)
  likelihood(data$df, centres, growths, params["p_background"])
}
dr_j_prior <- function(params, misc){
  #just unifom for now
  -10
}

mcmc <- drjacoby::run_mcmc(data = list(df = worldwide),
                           df_params = df_params,
                           loglike = dr_j_likelihood,
                           logprior = dr_j_prior,
                           burnin = 100,
                           samples = 100)
