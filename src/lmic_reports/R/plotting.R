
cumulative_deaths_plot <- function(country) {
  
d <- readRDS("ecdc_all.rds")

d$Continent <- countrycode::countrycode(d$Region, origin = 'country.name', destination = 'continent')
d$Continent[d$Region=="Eswatini"] <- "Africa"
d$Continent[d$Region=="United State of America"] <- "Americas"
d$Continent[d$Region=="Isle_of_Man"] <- "Europe"             
d$Continent[d$Region=="Kosovo"] <- "Europe"                  
d$Continent[d$Region=="Netherlands_Antilles"] <- "Americas"    
d$Continent[d$Region=="Saint_Lucia"] <- "Americas"             
d$Continent[d$Region=="South_Korea"] <- "Asia"             
d$Continent[d$Region=="United_States_of_America"] <- "Americas"

doubling <- function(double = 2, start = 10, xmax = 100) {
  
  x <- seq(0, xmax, 0.1)
  y <- start * 2^(x/double) 
  return(data.frame(x= x, y = y, 
                    Doubling = paste0("Every ", double, " Days")))
}

df <- group_by(d, Region) %>% 
  arrange(dateRep) %>% 
  mutate(Cum_Deaths = cumsum(deaths),
         Cum_Cases = cumsum(cases))

df$Region <- gsub("_" ," ", df$Region)

df_deaths <- df %>% 
  filter(Cum_Deaths > start) %>% 
  mutate(day_since = seq_len(n())-1)

doubling_lines_deaths <- do.call(rbind, lapply(c(2, 3, 5, 7), function(x){
  doubling(x, start = start, xmax = max(df_deaths$day_since))
}))

df_deaths_latest <- df_deaths[df_deaths$dateRep == max(df_deaths$dateRep),]
continent <- unique(df$Continent[df$Region == country])

gg_deaths <- ggplot(df_deaths[which(df_deaths$Continent == continent), ], aes(x=day_since, y=Cum_Deaths, group = Region)) + 
  geom_line(data = doubling_lines_deaths, aes(x=x, y=y, linetype = Doubling),alpha=0.3, inherit.aes = FALSE) +
  geom_line(show.legend = FALSE, color = "grey", alpha = 0.6) +
  geom_line(data = df_deaths[which(df_deaths$Region == country), ], color = "red", lwd = 2) +
  geom_point(data = df_deaths_latest[which(df_deaths_latest$Continent == continent), ], alpha = 0.5, show.legend = FALSE) + 
  #geom_text(data = df_deaths_latest[df_deaths_latest$Continent == continent, ], aes(label = Region), hjust = 0, nudge_x = 0.2, show.legend = FALSE) + 
  ggrepel::geom_text_repel(data =  df_deaths_latest[which(df_deaths_latest$Continent == continent), ],
                           aes(label = Region), show.legend = FALSE, min.segment.length = 1,nudge_x = 1) + 
  scale_y_log10(limits=c(start, max(df_deaths$Cum_Deaths[df_deaths$Continent == continent]))) +
  xlim(limits=c(0, max(df_deaths$day_since[df_deaths$Continent == continent]))) +
  theme_bw() +
  ylab("Cumulative Deaths (Logarithmic Scale)") +
  xlab(paste("Days Since", start, "Deaths"))

gg_deaths

}