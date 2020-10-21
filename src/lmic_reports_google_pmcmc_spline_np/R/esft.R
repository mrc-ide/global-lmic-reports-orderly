income_Rt <- function(date_0) {
  
  # parse for latest date
  url <- "https://github.com/mrc-ide/global-lmic-reports/tree/master/data"
  html <- xml2::read_html(url)
  links <- rvest::html_nodes(html, ".js-navigation-open.link-gray-dark")
  text <- rvest::html_text(links, "title")
  latest <- which.max(as.Date(head(substr(text, 1, 10), -1)))
  d_link <- paste0("https://github.com/mrc-ide/global-lmic-reports/raw/master/data/",text[latest])
  
  # download and read in data file
  tf <- tempfile()
  suppressMessages(download.file(d_link, tf, quiet = TRUE))
  tf2 <- tempfile()
  extract <- unzip(tf, exdir = tf2)
  dat <- as.data.frame(data.table::fread(extract))
  
  # subset to Rt values
  dat <- dat[dat$compartment == "Rt" & dat$scenario == "Maintain Status Quo",]
  wb_metadata <- read.csv("World_Bank_Country_Metadata.csv", fileEncoding="UTF-8-BOM", stringsAsFactors = TRUE)
  dat$income <- wb_metadata$income_group[match(dat$iso3c, wb_metadata$country_code)]
  filt_date <- min(date_0, max(dat$date))
  income_rt <- group_by(dat %>% filter(date == filt_date), income) %>% summarise(Rt = mean(y_mean))
  return(income_rt)
}


income_R0 <- function() {
  
  pars_init <- readRDS("pars_init.rds")
  wb_metadata <- read.csv("World_Bank_Country_Metadata.csv", fileEncoding="UTF-8-BOM", stringsAsFactors = TRUE)
  
  R0s <- unlist(lapply(pars_init, "[[", "R0"))
  isos <- unlist(lapply(pars_init, "[[", "iso3c"))
  income <- wb_metadata$income_group[match(isos, wb_metadata$country_code)]
  
  inc_r0s <- data.frame("R0" = R0s, "iso3c" = isos, "income" = income)
  R0s_income <- group_by(inc_r0s, income) %>% summarise(R0 = median(R0))
  
  return(R0s_income)

}

init_state <- function(deaths_removed, iso3c, seeding_cases = 5) {
  
  # get an initial
  pop <- squire::get_population(iso3c = iso3c, simple_SEIR = FALSE)
  init <- squire:::init_check_explicit(NULL, pop$n, seeding_cases = seeding_cases)
  
  if(deaths_removed > 0) {
    
  # work out how many deaths and where
  probs <- (squire:::probs$prob_hosp * squire:::probs$prob_severe * squire:::probs$prob_severe_death_treatment) +
    (squire:::probs$prob_hosp * (1-squire:::probs$prob_severe * squire:::probs$prob_non_severe_death_treatment))
  probs <- probs*pop$n
  probs <- probs/sum(probs)
  deaths <- as.numeric(t(rmultinom(1, deaths_removed, probs)))
  
  # approximate IFR for income group
  wb_metadata <- read.csv("World_Bank_Country_Metadata.csv", fileEncoding="UTF-8-BOM", stringsAsFactors = TRUE)
  income <- wb_metadata$income_group[match(iso3c, wb_metadata$country_code)]
  ifrs <- data.frame("income" = c("Low income", "Lower middle income", "Upper middle income", "High income"),
                     "ifr" = c(0.17, 0.31, 0.51, 1.02))
  ifr <- ifrs$ifr[ifrs$income == income]
  R <- rpois(1, deaths_removed*1/ifr/0.01)
  R <- as.numeric(t(rmultinom(1, R, rep(1/length(probs), length(probs)))))
  R <- R - deaths
  
  # and update the inital to reflect
  init$D <- deaths
  init$S <- init$S - R - deaths
  init$R <- R
  
  }
  
  return(init)
}
