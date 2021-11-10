get_future_Rt = function(model_out, forcast_days){
  #cut overall Rt values up by Rt period and calculate all positive and negative
  #trends

  #get values used
  mobilities <- model_out$interventions$R0_change
  pars <- model_out$replicate_parameters
  R0 <- pars[["R0"]]
  date_R0_change <- model_out$interventions$date_R0_change
  Rt_args <- model_out$pmcmc_results$inputs$Rt_args
  start_date <- pars$start_date
  #get indexes of times post lockdown
  rts_start_index <- min(which(date_R0_change > Rt_args$date_Meff_change))
  Rw_duration <- Rt_args$Rt_rw_duration
  #function to get the Rt from the effects and vice versa
  f_rt <- function(x, replicate){
    return(
      2*R0[replicate]*plogis(x)
    )
  }
  inv_f_rt <- function(Rt, replicate){
    return(
      qlogis(Rt/(2*R0[replicate]))
    )
  }
  fitting <- TRUE
  while(fitting){
    #for each replicate
    trends <- tibble(
      replicate = seq_along(R0)
    ) %>%
      rowwise() %>%
      mutate(
        Rt_list = list({
          #1.get Rts for replicate
          Rt <- squire:::evaluate_Rt_pmcmc(
            mobilities, R0[replicate], date_R0_change, pars[replicate,], Rt_args
          )
          #reduce to releveant post lockdown
          Rt <- Rt[rts_start_index:length(Rt)]
          #2. reverse the logistic transformation to get the effects
          Rt_effects <- inv_f_rt(Rt, replicate)
          #get the final Rt effect value for use later
          final_Rt_effect <- tail(Rt_effects, 1)
          #3. estimate linear trend within each 14 day chunk
          gradients <- unlist(lapply(seq(1, floor(length(Rt)/Rw_duration)), function(period){
            Rt_effects_period <- Rt_effects[
              seq(1, Rw_duration)*period
            ]
            lm(y~x, data = data.frame(x = seq_along(Rt_effects_period),
                                      y = Rt_effects_period))$coefficients[2]
          }))
          #4. get the median postive (Rt increases) and negative trends (Rt decreases)
          #return NA if we don't have one or the other
          if(sum(gradients > 0) == 0 | sum(gradients < 0) == 0){
            c(NA, NA)
          }
          positive_trend <- quantile(gradients[gradients > 0], 0.90)
          negative_trend <- quantile(gradients[gradients < 0], 0.10)
          #5. project total effect over projection period
          positive_effects <- final_Rt_effect + positive_trend*seq(1, forcast_days)
          negative_effects <- final_Rt_effect + negative_trend*seq(1, forcast_days)
          #6. calculate the Rt values on each day of the forcast
          positive_Rts <- f_rt(positive_effects, replicate)
          negative_Rts <- f_rt(negative_effects, replicate)
          #7. for simplicity we'll use the mean value
          positive_Rt <- mean(positive_Rts)
          negative_Rt <- mean(negative_Rts)
          #6. return as vector
          c(
            positive_Rt,
            negative_Rt
          )
        }),
        pessimistic = Rt_list[[1]],
        optimistic = Rt_list[[2]]
      ) %>%
      ungroup() %>%
      arrange(replicate) %>%
      select(!c(replicate, Rt_list))
    #check for NA values
    if(any(is.na(trends$pessimistic))|any(is.na(trends$optimistic))){
      #reduce RW by half
      Rw_duration <- Rw_duration/2
      if(Rw_duration < 3){
        warning("Unable to estimate scenario effects correctly for some replicates,
                removing the problematic replicates")
        #just use the trends without the nas
        trends <- na.omit(trends)
        if(nrow(trends) == 0){
          stop("No non-problematic trends avaiable")
        }
      } else {
        Rw_duration <- round(Rw_duration)
      }
    } else {
      fitting <- FALSE
    }
  }
  return(trends)
}

update_Rt <- function(model_out, Rt){
  pars <- model_out$replicate_parameters
  R0_change <- model_out$interventions$R0_change
  date_R0_change <- model_out$interventions$date_R0_change
  Rt_args <- model_out$pmcmc_results$inputs$Rt_args
  R0 <- pars[["R0"]]
  spline_effects <- pars[grepl("Rt_rw", names(pars))]
  to_replace <- tail(colnames(spline_effects), 1)
  #calculate spline values by reverse calculating final Rt etc
  final_Rt <- unlist(lapply(1:length(R0), function(x){
    tail(squire:::evaluate_Rt_pmcmc(
      R0_change, R0[x], date_R0_change, pars[x,], Rt_args
    ), 1)
  }))
  effect <- qlogis(Rt/(R0*2)) -  qlogis(final_Rt/(R0*2))
  out <- model_out
  out$replicate_parameters[, to_replace] <- out$replicate_parameters[, to_replace] - effect
  return(out)
}



# get_future_Rt = function(model_out){
#   #cut overall Rt values up by Rt period and calculate all positive and negative
#   #trends
#
#   #1. calculate Rt values for each iteration
#   mobilities <- model_out$interventions$R0_change
#   pars <- model_out$pmcmc_results$results
#   R0 <- pars[["R0"]]
#   date_R0_change <- model_out$interventions$date_R0_change
#   Rt_args <- model_out$pmcmc_results$inputs$Rt_args
#   start_date <- pars$start_date
#   #get indexes of times post lockdown
#   rts_start_index <- min(which(date_R0_change > Rt_args$date_Meff_change))
#   Rts <- lapply(seq_along(R0), function(replicate){
#     values <- squire:::evaluate_Rt_pmcmc(
#       rep(mobility, length(raw_mobilities)), R0[replicate], date_R0_change, pars[replicate,], Rt_args
#     )
#     values[rts_start_index:length(values)]
#   })
#   #figure out which indexes the spline effects will kick in
#   rts_index <- seq(rts_start_index + Rt_args$plateau_duration, length(raw_mobilities) - 14, by = Rt_args$Rt_rw_duration)
#   #limit mobilities to this period
#   mobilities <- tail(raw_mobilities[rts_start_index:length(raw_mobilities)],1)
#   #calculate Rts over time for each replicate * mobility
#   res <- expand.grid(
#     replicate = 1:length(R0),
#     mobility = mobilities
#   )
#   #get the possible Rt values for each spline update
#   res <- res %>%
#     rowwise() %>%
#     mutate(
#       rts = list(
#         squire:::evaluate_Rt_pmcmc(
#           rep(mobility, length(raw_mobilities)), R0[replicate], date_R0_change, pars[replicate,], Rt_args
#         )[rts_index]
#       )
#     ) %>%
#     ungroup()
#
#   #get values out of list
#   for(i in 1:length(rts_index)){
#     name <- paste0("Rt_", i)
#     res <- res %>%
#       mutate(
#         !! name := unlist(rts)[((1:nrow(res)) - 1)*length(rts_index) + i]
#       )
#   }
#
#   #now extend out our values
#   res <- res %>%
#     select(!rts) %>%
#     pivot_longer(contains("Rt_"), names_to = "rt_index", values_to = "Rt")
#   res <- res %>%
#     mutate(
#       rt_index = as.numeric(
#         unlist(strsplit(rt_index, "_"))[(1:nrow(res))*2]
#       )
#     )
#
#   #split into those about or below current final Rt for each replicate
#   #calculate the final Rt values
#   res <- res %>%
#     group_by(replicate) %>%
#     left_join(
#       tibble(
#         replicate = 1:length(R0),
#         final_Rt = unlist(lapply(1:length(R0), function(x){
#           tail(squire:::evaluate_Rt_pmcmc(
#             raw_mobilities, R0[x], date_R0_change, pars[x,], Rt_args
#           ), 1)
#         }))))
#
#   res <- rbind(
#     filter(res, Rt <= final_Rt) %>%
#       mutate(optimistic = TRUE),
#     filter(res, Rt >= final_Rt) %>%
#       mutate(optimistic = FALSE),
#   ) %>%
#     group_by(replicate, optimistic)
#
#   #find the median and the correspond moblity and date of non-mobility shift
#   #not the true median in the case of an even dataset, then we just take the higher
#   #value
#   return(arrange(res, optimistic, Rt) %>%
#            mutate(
#              median = (1:length(Rt)) == ceiling(length(Rt)/2 + 1),
#              date_change = date_R0_change[rt_index * Rt_args$Rt_rw_duration]
#            ) %>%
#            filter(median) %>%
#            arrange(replicate) %>%
#            ungroup() %>%
#            select(Rt, mobility, date_change, optimistic))
# }
