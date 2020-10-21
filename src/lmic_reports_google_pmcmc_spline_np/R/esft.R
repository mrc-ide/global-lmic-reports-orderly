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
  download.file(d_link, tf)  
  tf2 <- tempfile()
  extract <- unzip(tf, exdir = tf2)
  dat <- as.data.frame(data.table::fread(extract))
  
  # subset to Rt values
  dat <- dat[dat$compartment == "Rt" & dat$scenario == "Maintain Status Quo",]
  
  dat$income <- wb_metadata$income_group[match(dat$iso3c, wb_metadata$country_code)]
  income_rt <- group_by(dat %>% filter(date == date_0), income) %>% summarise(Rt = mean(y_mean))
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
