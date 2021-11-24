r_list_format <- function(out, date_0) {

  df <- nimue_format(out,
                      var_select = c("infections","deaths","hospital_demand",
                                     "ICU_demand", "D", "hospital_incidence","ICU_incidence"),
                      date_0 = date_0)

  pr <- nimue_format(out, var_select = c("S","R","D"), date_0 = date_0) %>%
    na.omit %>%
    pivot_wider(names_from = compartment, values_from = y) %>%
    mutate(y = sum(out$parameters$population)-D-R-S,
           compartment = "prevalence") %>%
    select(replicate, compartment, t, y, date)

  rbind(df, pr)

}
named_list <- function(...) {
  get <- as.character(match.call())
  l <- list(...)
  names(l) <- utils::tail(get, -1)
  return(l)
}

summarise_rt_futures <- function(Rt_futures){
  #for each trend type
  summarised <- lapply(names(Rt_futures), function(trend){
    #for each replicate
    do.call(
      rbind,
      lapply(seq_along(Rt_futures[[trend]]), function(rep){
        #just get the mean
        mean(Rt_futures[[trend]][[rep]])
      })
    )
  })
  #summarise and put into dataframe
  out <- data.frame(temp = 1
  )
  for(i in seq_along(Rt_futures)){
    out <- out %>%
      mutate(
        !! paste0(names(Rt_futures)[i], "_med") := median(summarised[[i]]),
        !! paste0(names(Rt_futures)[i], "_025") := quantile(summarised[[i]], 0.025),
        !! paste0(names(Rt_futures)[i], "_975") := quantile(summarised[[i]], 0.975)
      )
  }
  return(out %>%
           select(!temp))
}
