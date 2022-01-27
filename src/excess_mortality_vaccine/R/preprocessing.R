#A function to pre-process the data for model fitting
preprocess_for_fitting <- function(data){
  excluded_deaths <- 0
  #weeks with positive deaths
  death_weeks <- which(data$deaths > 0)
  #time between those weeks in weeks
  time_to_next <- diff(death_weeks)
  if(any(time_to_next > 6*4)){
    #if any breaks are longer than 6 months we consider it a possible re-introduction
    location_of_first_epidemic_end <- which(time_to_next > 6*4)
    if(length(location_of_first_epidemic_end) > 1){
      stop("Multiple reintroductions possible please adjust preprocess_for_fitting().")
    }
    #check what proportion of deaths occured before this point
    actual_index <- death_weeks[location_of_first_epidemic_end]
    if(sum(data[1:actual_index,]$deaths)/sum(data$deaths) < .1){
      #if less than 10% of total deaths occur before this point, we remove them
      #and only fit to the later epidemic
      excluded_deaths <- sum(data[1:actual_index,]$deaths)
      data[1:actual_index,]$deaths <- 0
    }
  }
  #in all other situtations we leave the data as is
  return(list(data, excluded_deaths))
}
