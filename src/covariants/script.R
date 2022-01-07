date_0 <- as.Date(date, "%Y-%m-%d")

#read in data
covariants_raw <- fromJSON(
  paste0("https://raw.githubusercontent.com/hodcroftlab/covariants",
         "/master/cluster_tables/EUClusters_data.json")
)

#bonaire is an overseas territory, so we'll remove it
covariants_raw$countries$Bonaire <- NULL
covariants_raw$countries$Kosovo <- NULL

#function to extract variant start date
get_variant_start_date <- function(var_name, variant_names, raw_data){
  start = TRUE
  #only interested in one variant so we need dates, total sequences and delta sequences
  for (country in names(raw_data$countries)) {
    country_df <- raw_data$countries[[country]]
    if (any(variant_names %in% names(country_df))) {
      variant_names_in_df <- intersect(variant_names, names(raw_data$countries[[country]]))
      country_df <- tibble(
        week = raw_data$countries[[country]][["week"]],
        var_seq = colSums(do.call(rbind, raw_data$countries[[country]][variant_names])),
        total_seq = raw_data$countries[[country]][["total_sequences"]]
      ) %>%
        mutate(
          var_prop = var_seq/total_seq
        )
    } else {
      country_df <- as.data.frame(raw_data$countries[[country]][c("week")])
      country_df$var_prop <- 0
    }
    country_df <- select(country_df, week, var_prop) %>%
      mutate(iso3c = countrycode(country, "country.name", "iso3c"),
             iso3c = if_else(is.na(iso3c), country, iso3c))
    if (start) {
      covariants_df <- country_df
      start <- FALSE
    } else {
      covariants_df <- rbind(covariants_df, country_df)
    }
  }
  covariants_df <- covariants_df %>%
    mutate(week = as.Date(week))
  #extract start date
  var_characteristics <- data.frame(iso3c = unique(covariants_df$iso3c))
  #calculate shift start (when delta >10% of cases) and
  var_characteristics <-
    covariants_df %>%
    group_by(iso3c) %>%
    arrange(week) %>%
    #only keep values where delta > 10% of sequences
    filter(var_prop > 0.05) %>%
    filter(week == min(week)) %>%
    rename(start_date = week) %>%
    select(iso3c, start_date) %>%
    right_join(
      var_characteristics,
      by = "iso3c"
    )
  #extract end date
  var_characteristics <-
    covariants_df %>%
    group_by(iso3c) %>%
    arrange(week) %>%
    #only keep values where delta > 10% of sequences
    filter(var_prop > 0.90) %>%
    filter(week == min(week)) %>%
    rename(end_date = week) %>%
    select(iso3c, end_date) %>%
    right_join(
      var_characteristics,
      by = "iso3c"
    )
  #check for consistency
  #only keep countries with start dates
  var_characteristics <- var_characteristics %>%
    filter(!is.na(start_date))
  #check if any end dates are before the start dates
  if (length(filter(var_characteristics, start_date > end_date) %>%
             pull(iso3c)) > 0) {
    #we just swap these around
    var_characteristics <- var_characteristics %>%
      mutate(
        start_date_old = start_date,
        end_date_old = end_date,
        start_date = min(start_date_old, end_date_old, na.rm = T),
        end_date = max(start_date_old, end_date_old, na.rm = T)
      ) %>%
      select(!c(start_date_old, end_date_old))
  }
  #check if the dates are equal
  if (length(filter(var_characteristics, start_date == end_date) %>%
             pull(iso3c)) > 0) {
    #for these countries we'll assume the same average shift time as across
    #the globe we'll assume the start/end date is the middle
    median_days <- var_characteristics %>%
      filter(start_date < end_date) %>%
      mutate(days = as.numeric(end_date - start_date)) %>%
      pull(days) %>%
      median()

    var_characteristics <- var_characteristics %>%
      mutate(
        old_start = start_date,
        old_end = end_date,
        start_date = if_else(
          identical(old_start, old_end),
          old_start - round(median_days / 2),
          old_start
        ),
        end_date = if_else(
          identical(old_start, old_end),
          old_end + round(median_days / 2),
          old_end
        )
      ) %>%
      select(!c(old_start, old_end))
  }
  #calculate duration of shift
  var_characteristics <- var_characteristics %>%
    mutate(shift_duration = as.numeric(end_date - start_date),
           #just use 60 for now to keep things consistent
           shift_duration = 60) %>%
    select(!end_date)
  #name variables with variant
  new_names <- paste0(var_name, c("_start_date", "_shift_duration"))
  var_characteristics[, new_names] <- var_characteristics[, c("start_date", "shift_duration")]
  var_characteristics[, c("start_date", "shift_duration")] <- NULL
  #add an indicator so we know what has been imputed later
  var_characteristics[[paste0(var_name, "_imputed")]] <- FALSE
  #return
  return(var_characteristics)
}

#get the timings for delta and omicron
variant_characteristics <- full_join(
  get_variant_start_date("delta", c("21A (Delta)", "21I (Delta)", "21J (Delta)"), covariants_raw),
  get_variant_start_date("omicron", c("21K (Omicron)", "21L (Omicron)"), covariants_raw)
)

#impute start dates for all countries based on means across sub region, continent
#and world if not there
all_iso3cs <- data.frame(
  iso3c = unique(squire::population$iso3c)
) %>%
  mutate(
    sub_region = countrycode(
      case_when(iso3c == "TWN" ~ "CHN",
                iso3c == "CHI" ~ "GBR",
                TRUE ~ iso3c),
      origin = "iso3c",
      destination = "un.regionsub.name"),
    continent = countrycode(
      case_when(iso3c == "TWN" ~ "CHN",
                iso3c == "CHI" ~ "GBR",
                TRUE ~ iso3c),
      origin = "iso3c",
      destination = "continent")
  )


parameters <- setdiff(names(variant_characteristics), "iso3c")
#function to add the median of the parameters by some measure
add_median <- function(data, join_on){
  if(join_on=="" & nrow(data %>% filter(
    if_any(
      all_of(parameters),
      ~is.na(.x)
    )
  ) ) > 0){
    data %>% filter(
      if_all(
        all_of(parameters),
        ~!is.na(.x)
      )
    ) %>%
      rbind(
        data %>% filter(
          if_any(
            all_of(parameters),
            ~is.na(.x)
          )
        ) %>%
          select(!all_of(parameters)) %>%
          cbind(
            data %>%
              summarise(across(
                all_of(parameters),
                ~median(.x, na.rm = T)
              ))
          )
      )
  } else if(nrow(data %>% filter(
    if_any(
      all_of(parameters),
      ~is.na(.x)
    )
  ) ) > 0){
    data %>% filter(
      if_all(
        all_of(parameters),
        ~!is.na(.x)
      )
    ) %>%
      rbind(
        data %>% filter(
          if_any(
            all_of(parameters),
            ~is.na(.x)
          )
        ) %>%
          select(!all_of(parameters)) %>%
          left_join(
            data %>%
              group_by(across(all_of(join_on))) %>%
              summarise(across(
                all_of(parameters),
                ~median(.x, na.rm = T)
              )),
            by = join_on
          )
      )
  } else{
    data
  }
}

#for data that is missing add sub_region medians
variant_characteristics <- all_iso3cs %>% #add countries with data
  left_join(variant_characteristics,
            by = "iso3c") %>% #add sub_region median
  add_median("sub_region") %>% #add contient median
  add_median("continent") %>% #add world median
  add_median("") %>%
  select(!c(continent, sub_region)) %>% #update imputed indicator
  mutate(
    across(
      ends_with("_imputed"),
      ~if_else(is.na(.x), TRUE, as.logical(.x))
    )
  )


#delta characteristics
variant_characteristics <- variant_characteristics %>%
  mutate(
    #assume 27% immune escape
    delta_immune_escape = 0.27,
    #calculate required changes to the recovery duration
    delta_required_dur_R = 1 / (
      (delta_shift_duration / 365 - log(1 - delta_immune_escape)) / delta_shift_duration
    ),
    #add increased hospitalization
    #(https://www.thelancet.com/journals/laninf/article/PIIS1473-3099(21)00475-8/
    #fulltext#seccestitle150)
    delta_prob_hosp_multiplier = 1.45,
    #decreased vaccine efficacy
    delta_ve_i_low = 0.224,
    delta_ve_i_high = 0.646,
    delta_ve_d_low = 0.75,
    delta_ve_d_high = 0.94,
  )

#calculate vaccine efficacies based on report 49 and 50 and https://www.medrxiv.org/content/10.1101/2021.12.14.21267615v1
#to do

#add omicron adjustments
variant_characteristics <- variant_characteristics %>%
  mutate(
    #assume 50% immune escape same for now
    omicron_immune_escape = 0.50,
    #calculate required changes to the recovery duration
    omicron_required_dur_R = 1 / (
      (omicron_shift_duration / 365 - log(1 - omicron_immune_escape)) / omicron_shift_duration
    ),
    #add decreased hospitalization
    omicron_prob_hosp_multiplier = 0.6,
    #decreased vaccine efficacy, won't touch disease efficacy for now
    omicron_ve_i_low = 0,
    omicron_ve_i_high = 0
  )

#print a plot to check our results
dir.create("calibration")
pdf("calibration/plot.pdf")
print(
  ggplot(variant_characteristics %>%
    arrange(delta_start_date)) +
    geom_segment(aes(
      y = fct_reorder(iso3c, delta_start_date),
      yend = fct_reorder(iso3c, delta_start_date),
      x = delta_start_date,
      xend = delta_start_date + delta_shift_duration
    ),
    alpha = 0.5,
    size = 2) +
    labs(
      x = "Date",
      y = NULL,
      title = "Duration of Delta Shift"
    )
)
print(
  ggplot(variant_characteristics %>%
           arrange(omicron_start_date)) +
    geom_segment(aes(
      y = fct_reorder(iso3c, omicron_start_date),
      yend = fct_reorder(iso3c, omicron_start_date),
      x = omicron_start_date,
      xend = omicron_start_date + omicron_shift_duration
    ),
    alpha = 0.5,
    size = 2) +
    labs(
      x = "Date",
      y = NULL,
      title = "Duration of Delta Shift"
    )
)
dev.off()

#save data
saveRDS(variant_characteristics, "variant_characteristics.rds")
