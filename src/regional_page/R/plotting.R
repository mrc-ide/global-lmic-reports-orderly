
cumulative_deaths_plot_continent_projections <- function(continent, today, data, ecdc) {
  
  ## handle arugments coming in
  if(!continent %in% c("Asia","Europe","Africa","Americas","Oceania")) {
    stop("continent not matched")
  }
  today <- as.Date(today)
  
  # identify lmics
  rl <- readLines("_navbar.html")
  lmics <- gsub("(.*reports/)(\\w\\w\\w)(\".*)","\\2",grep("reports/(\\w\\w\\w)\"",rl, value =TRUE))
  
  # Are we presnting the surge if it's there
  data <- group_by(data, iso3c) %>% 
    filter(scenario == (if("Surged Maintain Status Quo" %in% unique(scenario)) {
      "Surged Maintain Status Quo"
    } else {
      "Maintain Status Quo"
    })) %>% ungroup()
  
  # create dataset
  slim <- data %>% 
    mutate(date = as.Date(.data$date)) %>% 
    filter(date > (today)) %>%
    filter(date < (today+28)) %>% 
    select(date, compartment, y_mean, y_025, y_975, country, iso3c) %>% 
    mutate(observed = FALSE) %>% 
    rename(y = y_mean)
  
  slim <- slim[,c("date", "y", "country", "iso3c", "observed", "compartment", "y_025", "y_975")]
  
  # handle ecdc
  names(ecdc)[names(ecdc) %in% c("dateRep","Region", "countryterritoryCode")] <- c("date", "country", "iso3c")
  ecdc <- ecdc %>% 
    mutate(date = as.Date(date),
           observed = TRUE,
           compartment = "deaths") %>% 
    rename(y = deaths) %>% 
    select(date, y, country, iso3c, observed, compartment) %>% 
    mutate(y_025 = NA, 
           y_975 = NA)
  
  
  df <- as.data.frame(do.call(rbind, list(as.data.frame(ecdc), slim)), stringsAsFactors = FALSE)
  
  d <- df
  d$country[d$country=="Congo"] <- "Republic of Congo"
  d$country[d$country=="United_Republic_of_Tanzania"] <- "Tanzania"
  d$country[d$country=="CuraÃ§ao"] <- "Curacao"
  start <- 10
  
  suppressWarnings(d$Continent <- countrycode::countrycode(d$country, origin = 'country.name', destination = 'continent'))
  d$Continent[d$country=="Eswatini"] <- "Africa"
  d$Continent[d$country=="United State of America"] <- "Americas"
  d$Continent[d$country=="Isle_of_Man"] <- "Europe"             
  d$Continent[d$country=="Kosovo"] <- "Europe"                  
  d$Continent[d$country=="Netherlands_Antilles"] <- "Americas"    
  d$Continent[d$country=="Saint_Lucia"] <- "Americas"             
  d$Continent[d$country=="South_Korea"] <- "Asia"             
  d$Continent[d$country=="United_States_of_America"] <- "Americas"
  
  doubling <- function(double = 2, start = 10, xmax = 100) {
    
    x <- seq(0, xmax, 0.1)
    y <- start * 2^(x/double) 
    return(data.frame(x= x, y = y, 
                      Doubling = paste0("Every ", double, " Days")))
  }
  
  d <- d[d$compartment=="deaths",]
  d$date <- as.Date(d$date)
  
  
  df <- group_by(d, iso3c) %>% 
    arrange(date) %>%  
    mutate(Cum_Deaths = cumsum(y))
  
  df$country <- gsub("_" ," ", df$country)
  df <- df[which(df$iso3c %in% unique(df$iso3c[which(df$Cum_Deaths>10 & df$observed)])), ]
  
  df_deaths <- df %>% 
    filter(Cum_Deaths > start) %>% 
    mutate(day_since = seq_len(n())-1)
  
  doubling_lines_deaths <- do.call(rbind, lapply(c(2, 3, 5, 7), function(x){
    doubling(x, start = start, xmax = max(df_deaths$day_since))
  }))
  
  df_deaths_latest <- df_deaths[df_deaths$date == max(df_deaths$date),]
  df_deaths$region <- countrycode::countrycode(df_deaths$iso3c, "iso3c", "region23")
  df_deaths_latest$region <- countrycode::countrycode(df_deaths_latest$iso3c, "iso3c", "region23")
  
  
  gg_deaths <- ggplot(df_deaths[which(df_deaths$Continent == continent), ], aes(x=day_since, y=Cum_Deaths, group = country)) + 
    #geom_line(data = doubling_lines_deaths, aes(x=x, y=y, linetype = Doubling), inherit.aes = FALSE, color = "black") +
    geom_line(show.legend = FALSE, color = "grey", alpha = 0.3) +
    geom_line(data = df_deaths[which(df_deaths$Continent %in% continent & df_deaths$iso3c %in% lmics & df_deaths$observed),], 
              mapping = aes(color = country), show.legend = FALSE, lwd = 0.75) +
    geom_line(data = df_deaths[which(df_deaths$Continent %in% continent & df_deaths$iso3c %in% lmics),], 
              linetype = "dashed", mapping = aes(color = country), show.legend = FALSE, lwd = 0.75) +
    #geom_point(data = df_deaths[which(df_deaths$country %in% country[1:7]),], mapping = aes(color = Continent)) +
    geom_point(data = df_deaths_latest[which(df_deaths_latest$Continent %in% continent & df_deaths_latest$iso3c %in% lmics), ], 
               mapping = aes(color = country), alpha = 0.5, show.legend = FALSE) + 
    ggrepel::geom_text_repel(data =  df_deaths_latest[which(df_deaths_latest$Continent %in% continent & df_deaths_latest$iso3c %in% lmics), ],
                             aes(label = country), show.legend = FALSE, min.segment.length = 0.1,nudge_x = 1,nudge_y = -0.1) + 
    scale_y_log10(limits=c(start, max(df_deaths$Cum_Deaths[df_deaths$Continent %in% continent & df_deaths$iso3c %in% lmics])), 
                  labels = scales::comma) +
    xlim(limits=c(0, max(df_deaths[which(df_deaths$Continent %in% continent & df_deaths$iso3c %in% lmics),]$day_since))) +
    theme_bw() +
    scale_color_hue("country",l = 50, c = 90) +
    #scale_linetype(name = "Doubling Time:") +
    ylab("Cumulative Deaths (Logarithmic Scale)") +
    xlab(paste("Days Since", start, "Deaths")) + 
    ggtitle(continent)
  
  gg_deaths
  
}

full_firework_plot <- function() {
  
  data <- readRDS("all_data.rds")
  
  ecdc <- readRDS("ecdc_all.rds")
  
  plots <- lapply(c("Asia","Europe","Africa","Americas","Oceania"), 
                  cumulative_deaths_plot_continent_projections, 
                  today = date, 
                  data = data, 
                  ecdc = ecdc)
  plotted <- lapply(plots[1:4], function(x){x+theme(legend.position = "none")})
  leg <- cowplot::get_legend(plots[[1]] + theme(legend.position = "top"))
  main <- cowplot::plot_grid(plotlist = plotted[1:4], ncol = 2)
  get <- cowplot::plot_grid(main,leg,ncol=1,rel_heights = c(1, 0.05))
  return(get)
}

one_firework_plot <- function(cont) {
  
  data <- readRDS("all_data.rds")
  
  ecdc <- readRDS("ecdc_all.rds")
  
  plot <- cumulative_deaths_plot_continent_projections(
    continent = cont, 
    today = date, 
    data = data, 
    ecdc = ecdc) + theme(plot.title = element_blank())
  
  return(plot)
  
}

forecasted_deaths_bar <- function(cont, today) {
    
    ## handle arugments coming in
    if(!continent %in% c("Asia","Europe","Africa","Americas","Oceania")) {
      stop("continent not matched")
    }
    today <- as.Date(today)
  
    data <- readRDS("all_data.rds")
    # create dataset
    slim <- data %>% 
      mutate(date = as.Date(.data$date)) %>% 
      filter(date == (today+28)) %>%
      filter(scenario == "Maintain Status Quo") %>% 
      select(date, compartment, y_mean, y_025, y_975, country, iso3c) %>% 
      mutate(observed = FALSE) %>% 
      rename(y = y_mean)
    
    d <- slim
    d$country[d$country=="Congo"] <- "Republic of Congo"
    d$country[d$country=="United_Republic_of_Tanzania"] <- "Tanzania"
    d$country[d$country=="CuraÃ§ao"] <- "Curacao"
    start <- 10
    
    suppressWarnings(d$Continent <- countrycode::countrycode(d$country, origin = 'country.name', destination = 'continent'))
    d$Continent[d$country=="Eswatini"] <- "Africa"
    d$Continent[d$country=="United State of America"] <- "Americas"
    d$Continent[d$country=="Isle_of_Man"] <- "Europe"             
    d$Continent[d$country=="Kosovo"] <- "Europe"                  
    d$Continent[d$country=="Netherlands_Antilles"] <- "Americas"    
    d$Continent[d$country=="Saint_Lucia"] <- "Americas"             
    d$Continent[d$country=="South_Korea"] <- "Asia"             
    d$Continent[d$country=="United_States_of_America"] <- "Americas"
    
    d$country <- as.character(d$country)
    d$country <- factor(d$country, levels = sort(unique(d$country)))
    d <- d %>% 
      filter(Continent == cont & compartment == "deaths")
    
    d <- filter(d, y>=sort(d$y,decreasing = TRUE)[10])

    ggplot(d, aes(x=as.factor(country), y = y)) + 
      geom_bar(stat = "identity") + 
      geom_errorbar(mapping = aes(ymin = y_025, ymax = y_975)) + 
      coord_flip() + 
      scale_y_continuous(labels = scales::comma) + 
      xlab("") + ylab("Projected Deaths in 28 days") + 
      theme_bw()
      
}

summaries_forecasts_plot <- function(sums, cont) {
  
  sums$country <- as.character(sums$country)
  sums$country <- factor(sums$country, levels = sort(unique(sums$country),decreasing = TRUE))
  sums <- mutate(sums) %>% 
    mutate(value = ceiling(value)) %>% 
  filter(continent == cont)

  gg <- ggplot(sums[sums$variable %in% c("hospital_28","icu_28"),] %>% 
                 filter(value > 1), 
               aes(x = country, y = value, color = variable, fill = variable)) + 
    geom_bar(stat="identity",position = "dodge", width = 0.5) + 
    scale_y_continuous(labels = scales::comma) + 
    scale_fill_manual("", labels = c("Estimated Hospital Beds\nNeeded in 28 days",
                                     "Estimated ICU Beds\nNeeded in 28 days"),
                      values = c("#440154FF", "#21908CFF")) +
    scale_color_manual("", labels = c("Estimated Hospital Beds\nNeeded in 28 days",
                                      "Estimated ICU Beds\nNeeded in 28 days"),
                       values = c("#440154FF", "#21908CFF")) + 
    theme_bw() +
    xlab("") +
    ylab("") + 
    theme(legend.key = element_rect(size = 2),
          legend.key.size = unit(2, 'lines')) + 
    scale_y_log10() +
    coord_flip()
  
  gg
}

rt_plot <- function(cont) {
  
  today <- date
  sum_rt <- readRDS("all_data.rds")
  sum_rt$continent <- countrycode::countrycode(sum_rt$iso3c, "iso3c", "continent")
  sum_rt <- sum_rt %>% filter(continent == cont)
  sum_rt <- sum_rt %>% filter(date <= today)
  
  ggplot(sum_rt[sum_rt$compartment == "Reff" & sum_rt$scenario == "Maintain Status Quo",], 
         aes(x=as.Date(date), ymin=y_025, ymax = y_975, group = iso3c, fill = iso3c)) +
    geom_line(aes(y = y_median), color = "#48996b") +
    geom_ribbon(fill = "#96c4aa") +
    geom_ribbon(mapping = aes(ymin = y_25, ymax = y_75), fill = "#48996b") +
    geom_hline(yintercept = 1, linetype = "dashed") +
    theme_bw() +
    theme(axis.text = element_text(size=12)) +
    xlab("") +
    facet_wrap(~country, ncol = 6) +
    scale_x_date(breaks = "3 week",
                 date_labels = "%d %b") + 
    theme(legend.position = "none") +
    theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1, colour = "black", size = 8),
          axis.text.y = ggplot2::element_text(colour = "black", size = 8),
          panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          strip.background = element_blank(), 
          panel.background = element_blank(), axis.line = element_line(colour = "black")
    )
}


rt_continental_plot <- function(cont) {
  
  today <- date
  sum_rt <- readRDS("all_data.rds")
  sum_rt$continent <- countrycode::countrycode(sum_rt$iso3c, "iso3c", "continent")
  sum_rt <- sum_rt %>% filter(continent == cont)
  sum_rt <- sum_rt %>% filter(date <= today)
  
  
  ggplot(sum_rt[sum_rt$compartment == "Reff" & sum_rt$scenario == "Maintain Status Quo",], 
         aes(x=as.Date(date), ymin=y_025, ymax = y_975, group = iso3c, fill = iso3c)) +
    geom_ribbon(fill = "#96c4aa") +
    geom_line(color = "#48996b") +
    geom_ribbon(mapping = aes(ymin = y_25, ymax = y_75), fill = "#48996b") +
    geom_hline(yintercept = 1, linetype = "dashed") +
    theme_bw() +
    theme(axis.text = element_text(size=12)) +
    xlab("") +
    facet_wrap(~country, ncol = 6) +
    scale_x_date(breaks = "3 week",
                 date_labels = "%d %b") + 
    theme(legend.position = "none") +
    theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1, colour = "black", size = 8),
          axis.text.y = ggplot2::element_text(colour = "black", size = 8),
          panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          strip.background = element_blank(), 
          panel.background = element_blank(), axis.line = element_line(colour = "black")
    ) + geofacet::facet_geo(~code, grid = "africa_countries_grid1",label = "name")
  
  
}