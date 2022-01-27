add_delta_characteristics <- function(data, adjust_omicron){
  #load delta data
  variant_characteristics <- readRDS("variant_characteristics.rds") %>% ungroup()

  #add to countries
  df <-  data %>% left_join(
      variant_characteristics %>%
        mutate(delta_shift_start = delta_start_date,
               delta_shift_end = delta_shift_start + delta_shift_duration) %>%
        select(iso3c, delta_shift_start, delta_shift_end, contains("delta_ve")),
      by = "iso3c"
      ) %>% #calculate where along the shift from wild to Delta we are
      mutate(shift = case_when(
        date_vaccine_change > delta_shift_end ~ 1,
        date_vaccine_change < delta_shift_start ~ 0,
        TRUE ~ as.numeric(date_vaccine_change - delta_shift_start)/
          as.numeric(delta_shift_end - delta_shift_start)
          )
        ) %>%
      mutate(
        ve_i_high = ve_i_high*(1 - shift) + delta_ve_i_high*shift,
        ve_d_high = ve_d_high*(1 - shift) + delta_ve_d_high*shift,
        ve_i_low = ve_i_low*(1 - shift) + delta_ve_i_low*shift,
        ve_d_low = ve_d_low*(1 - shift) + delta_ve_d_low*shift
      ) %>%
      select(!c(delta_ve_i_high, delta_ve_d_high, delta_ve_i_low, delta_ve_d_low, shift))

    #add omicron adjustments if also needed
    if(adjust_omicron){
      test <- df %>% left_join(
        variant_characteristics %>%
          filter(!omicron_imputed) %>%
          mutate(omicron_shift_start = omicron_start_date,
                 omicron_shift_end = omicron_shift_start + omicron_shift_duration) %>%
          select(iso3c, omicron_shift_start, omicron_shift_end, contains("omicron_ve")),
        by = "iso3c"
      ) %>% #calculate where along the shift from wild to Delta we are
        mutate(shift = case_when(
          is.na(omicron_shift_end) & is.na(omicron_shift_start) ~ as.numeric(NA),
          date_vaccine_change > omicron_shift_end ~ 1,
          date_vaccine_change < omicron_shift_start ~ 0,
          TRUE ~ as.numeric(date_vaccine_change - omicron_shift_start)/
            as.numeric(omicron_shift_end - omicron_shift_start)
        )
        ) %>%
        mutate(
          ve_i_high = if_else(is.na(shift), ve_i_high, ve_i_high*(1 - shift) + omicron_ve_i_high*shift),
          ve_i_low = if_else(is.na(shift), ve_i_low, ve_i_low*(1 - shift) + omicron_ve_i_low*shift)
        ) %>%
        select(!c(omicron_ve_i_high, omicron_ve_i_low, shift))
    }

  return(df)
}
