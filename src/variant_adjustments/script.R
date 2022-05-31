date_0 <- as.Date(date, "%Y-%m-%d")

## Variant Timings

#read in data from CoVariants:
covariants_raw <- jsonlite::fromJSON(
  paste0("https://raw.githubusercontent.com/hodcroftlab/covariants",
         "/master/cluster_tables/EUClusters_data.json")
)

#code sequencing names to VoC, this should be checked semi-regularly
variants <- list(
  "Delta" = c("21A (Delta)", "21I (Delta)", "21J (Delta)"),
  "Omicron" = c("21K (Omicron)", "21L (Omicron)")
)

#format into weekly time-series data
covariants_df <- format_covariants(covariants_raw, variants)

#get start dates for the two sequences
variant_start_dates <- full_join(
  get_start_date(covariants_df, "Delta"),
  get_start_date(covariants_df, "Omicron"),
  by = "iso3c"
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

#for data that is missing add sub_region medians
variant_start_dates <- all_iso3cs %>% #add countries with data
  left_join(variant_start_dates,
            by = "iso3c") %>%
  mutate(Delta_start_date = as.Date(Delta_start_date),
         Omicron_start_date = as.Date(Omicron_start_date)) %>%
  add_median(c("sub_region", "continent", "")) %>%
  select(!c(continent, sub_region))

#set the characteristics for each variant
variant_characteristics <- list(
  Wild = list(
    prob_hosp_multiplier = 1,
    prob_severe_multiplier = 1,

    ve_infection = c(0.6, 0, 0.8, 0.004340038, 0.004340038),
    ve_disease = c(0.8, 0.022339931, 0.98, 0.98, 0.141764382),
    ve_transmission = 0.5,

    immune_escape = 0,
    dur_ICU = 14.8,
    dur_ICU_death = 11.1,
    dur_hosp = 9,
    dur_hosp_death = 9,

    dur_V = c(1/0.009339697, 1/0.008123934, 1/0.006766174)
  ),
  Delta = list(
    prob_hosp_multiplier =  1.45,
    prob_severe_multiplier = 1,

    ve_infection = c(0.224, 0, 0.646, 0.001250519, 0.004340204),
    ve_disease = c(0.75, 0.010415902, 0.94, 0.94, 0.061241989),
    ve_transmission = 0.5,

    immune_escape = 0.27,
    dur_ICU = 14.8,
    dur_ICU_death = 11.1,
    dur_hosp = 5,
    dur_hosp_death = 9,

    shift_duration = 60,

    dur_V = c(1/0.009561369, 1/0.010645451, 1/0.011359099)
  )
)
#calculate omicron values, slightly more complex
hr_hosp <- 0.59
hr_severe_hosp <- 0.34/hr_hosp #adjust severity hazard ratio by reduction in hospitalisation
ve_d_2 <-  mean(0.31, 0.22)/hr_hosp #as we have no value for first dose scale by change in delta values
change_delta <- variant_characteristics$Delta$ve_disease[3]/
  variant_characteristics$Delta$ve_disease[1]
dur_R <- 365
dur_s <- 45

new_dur <- -dur_s/(log(1 - 5.41*(1-exp(-dur_s/dur_R))))
immune_escape <- 1 - exp(dur_s * (1/dur_R - 1/new_dur))

variant_characteristics$Omicron <- list(
    prob_hosp_multiplier = hr_hosp,
    prob_severe_multiplier = hr_severe_hosp,

    ve_infection = c(0, 0, 0.1, 0, 0),
    ve_disease = c(ve_d_2/change_delta, 0.002935928, ve_d_2, 0.043389992, 0.002869262),
    ve_transmission = 0.5,

    immune_escape = immune_escape,
    dur_ICU = 6,
    dur_ICU_death = 11.1,
    dur_hosp = 3,
    dur_hosp_death = 9,

    shift_duration = dur_s,
    
    dur_V = c(1/0.013728390, 1/0.013242636, 1/0.017599820 )
  )

#produce explainer

#WIP

#create list for each country
iso3cs <- unique(variant_start_dates$iso3c)
names(iso3cs) <- iso3cs
variant_characteristics <- map(iso3cs, function(iso3c){
  variant_characteristics$Delta$start_date <- variant_start_dates$Delta_start_date[
    variant_start_dates$iso3c == iso3c
  ]
  variant_characteristics$Omicron$start_date <- variant_start_dates$Omicron_start_date[
    variant_start_dates$iso3c == iso3c
  ]
  variant_characteristics
})

#save data
saveRDS(variant_characteristics, "variant_characteristics.rds")
