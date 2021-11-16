r_list_format <- function(out, date_0) {

  df <- nim_sq_format(out,
                      var_select = c("infections","deaths","hospital_demand",
                                     "ICU_demand", "D", "hospital_incidence","ICU_incidence"),
                      date_0 = date_0)

  pr <- nim_sq_format(out, var_select = c("S","R","D"), date_0 = date_0) %>%
    na.omit %>%
    pivot_wider(names_from = compartment, values_from = y) %>%
    mutate(y = sum(out$parameters$population)-D-R-S,
           compartment = "prevalence") %>%
    select(replicate, compartment, t, y, date)

  rbind(df, pr)

}
