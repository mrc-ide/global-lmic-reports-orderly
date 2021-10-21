add_delta_characteristics <- function(data){
  #load delta data
  delta_characteristics <- readRDS("delta_characteristics.Rds") %>% ungroup()

  #add to countries
    data %>% left_join(
      delta_characteristics %>%
        mutate(shift_start = start_date,
               shift_end = shift_start + shift_duration) %>%
        select(iso3c, shift_start, shift_end, contains("ve")),
      by = "iso3c"
      ) %>% #calculate where along the shift from wild to Delta we are
      mutate(shift = case_when(
        date_vaccine_change > shift_end ~ 1,
        date_vaccine_change < shift_start ~ 0,
        TRUE ~ as.numeric(date_vaccine_change - shift_start)/
          as.numeric(shift_end - shift_start)
          )
        ) %>%
      mutate(
        ve_i_high = ve_i_high*(1 - shift) + ve_i_high_d*shift,
        ve_d_high = ve_d_high*(1 - shift) + ve_d_high_d*shift,
        ve_i_low = ve_i_low*(1 - shift) + ve_i_low_d*shift,
        ve_d_low = ve_d_low*(1 - shift) + ve_d_low_d*shift
      ) %>%
      select(!c(ve_i_high_d, ve_d_high_d, ve_i_low_d, ve_d_low_d, shift))
}
