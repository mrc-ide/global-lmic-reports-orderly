format_covariants <- function(raw_data, variants){
  #give the variants names so we don't get hundreds of messages about renaming
  #also add total sequences so it doesn't error when no sequences are found
  variants <- map(variants, function(x){
    x <- c(x, "total_sequences")
    names(x) <- seq_along(x)
    return(x)
  })
  #create a data frame for each country
  map_dfr(
    raw_data$countries,
    function(lists){
      df <- tibble(
        week = lists$week,
        sequences = lists$total_sequences
      ) %>%
        #get variants
        cbind(
          map_dfc(variants, ~rowSums(map_dfc(.x, ~lists[[.x]])))
        ) %>%
        mutate(
          across(all_of(names(variants)),
                 ~(.x - sequences)/sequences)
        )
    },
    .id = "iso3c"
  ) %>%
    filter(!iso3c %in% c("Bonaire", "Kosovo")) %>%
    mutate(iso3c = countrycode::countrycode(iso3c, origin = "country.name",
                                            destination = "iso3c"))
}

#This can be improve, once fitting issues resolved
get_start_date <- function(df, variant){
  # df %>%
  #   group_by(iso3c) %>%
  #   filter(sequences > 10) %>%
  #   filter(any(.data[[variant]] > 0.9)) %>%
  #   #isolate the period where variant first became dominant (greater than 0.95)
  #   arrange(iso3c, week) %>%
  #   filter(week <= min(week[.data[[variant]] > 0.9])) %>%
  #   #now take the last week where sequences where less than 5%
  #   filter(any(.data[[variant]] <= 0.05)) %>%
  #   filter(week == max(week[.data[[variant]] <= 0.05])) %>%
  #   select(iso3c, week) %>%
  #   rename(!! paste0(variant, "_start_date") := week)
  df %>%
    group_by(iso3c) %>%
    filter(.data[[variant]] > 0.05) %>%
    #isolate the period where variant first became dominant (greater than 0.95)
    arrange(iso3c, week) %>%
    mutate(
      problematic = min(week[.data[[variant]] > 0.9]) == min(week),
      problematic = if_else(is.na(problematic),
                            FALSE,
                            problematic),
      week = if_else(
        problematic,
        as.character(as.Date(min(week)) - 28), #keeps it the same as before
        min(week)
      )
    ) %>%
    select(iso3c, week) %>%
    unique() %>%
    rename(!! paste0(variant, "_start_date") := week)
}

#function to add the median of the parameters by some measure
add_median <- function(data, join_on){
  for(group in join_on){
    data <- data %>%
      ungroup() %>%
      group_by(across(any_of(group))) %>%
      mutate(
        across(contains("_start_date"),
               ~if_else(
                 is.na(.x),
                 median(.x, na.rm = TRUE),
                 .x
               ))
      )
  }
  return(data)
}

#function to print variant characetistcs
format_variant_characteristics <- function(characteristics){
  make_percentage <- function(x){
    paste0(signif(x*100, 3), "%")
  }
  #make data frame
  tibble(
    `Multiplier on Probability of Hospitalisation` = characteristics$prob_hosp_multiplier,
    `Multiplier on Probability of Severe Symptoms` = characteristics$prob_severe_multiplier,
    `Vaccine Efficacy against Infection` = paste0(c("First: ", "Second: "), make_percentage(characteristics$ve_infection), collapse = ", "),
    `Vaccine Efficacy against Disease` = paste0(c("First: ", "Second: "), make_percentage(characteristics$ve_disease), collapse = ", "),
    `Vaccine Reduction in Transmission` = make_percentage(characteristics$ve_transmission),
    `Immune Escape Percentage` = make_percentage(characteristics$immune_escape),
    `Duration of Treatment if non-severe, given survival` = characteristics$dur_hosp,
    `Duration of Treatment if non-severe, given death` = characteristics$dur_hosp_death,
    `Duration of Treatment if Severe, given survival` = characteristics$dur_ICU,
    `Duration of Treatment if Severe, given death` = characteristics$dur_ICU_death,
  )
}
