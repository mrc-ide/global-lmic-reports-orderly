get_dominant_vaccines <- function(iso3cs){
  url <- xml2::read_html("https://github.com/owid/covid-19-data/tree/master/public/data/vaccinations/country_data")
  ctries <- rvest::html_text(rvest::html_nodes(url, ".js-navigation-open.Link--primary"))
  links <- paste0("https://raw.githubusercontent.com/owid/covid-19-data/master/public/data/vaccinations/country_data/", ctries)
  links <- urltools::url_encode(links)
  dat <- lapply(links, read.csv)
  df <- data.frame("country" = gsub(".csv","", ctries, fixed = TRUE),
                   "vaccines" = unlist(lapply(dat, function(x){tail(x$vaccine,1)}))) %>%
    mutate(
      country = case_when(
        country %in% c("England", "Northern Ireland", "Scotland", "Wales") ~ "United Kingdom",
        country == "Timor" ~ "Timor Leste",
        TRUE ~ country
      ),
      country = if_else(
        country == "Kosovo",
        "XKX",
        countrycode(country, origin = "country.name", destination = "iso3c")
      )
    ) %>%
    filter(country %in% iso3cs) %>%
    rename(iso3c = country) %>%
    unique() %>%
    group_by(iso3c) %>%
    summarise(
      vaccines = paste0(vaccines, collapse = ", ")
    )

  #add missing values
  df <- df %>%
    rbind(
      tibble(
        iso3c = c("GUF", "GLP", "MTQ", "MYT", "REU"),
        vaccines = df %>%
          filter(iso3c == "FRA") %>%
          pull(vaccines)
      )
    )%>%
    rbind(
      vaccines = df %>%
        filter(iso3c == "GBR") %>%
        mutate(iso3c = "CHI")
    ) %>%
    rbind(
      tibble(
        iso3c = c("GUM", "FSM", "PRI", "VIR"),
        vaccines = c("Pfizer/BioNTech, Moderna, Johnson&Johnson",
                     "Moderna, Johnson&Johnson",
                     "Pfizer/BioNTech, Moderna, Johnson&Johnson",
                     "Pfizer/BioNTech, Moderna")
      )
    )

  #check for missing values

  missing <- setdiff(iso3cs, c(df$iso3c, "PRK", "ERI", "ESH"))
  if(length(missing) > 0){
    print(paste0(countrycode(missing, "iso3c", "country.name")), " is missing data")
  }
  #split into distinct categories
  vaccine_types <- unique(unlist(str_split(unique(df$vaccines), ", ")))
  #simplify to types
  vaccine_cat <- map_chr(vaccine_types, ~simplify_vaccine(.x))
  unique_cat <- unique(vaccine_cat)
  names(unique_cat) <- unique_cat
  #reduce
  df <- df %>%
    select(!vaccines) %>%
    cbind(
      map_dfc(unique_cat, function(category){
        reduce(vaccine_types[vaccine_cat == category],
               function(curr, new_vaccine_type){
                 curr + str_detect(df$vaccines, new_vaccine_type) > 0
               },
               .init = FALSE
        )
      })
    )#nK and Eritrea, Western Sahara have no vaccines (assume subunit vaccine)
  df <- df %>%
    rbind(
      df %>%
        filter(iso3c %in% c("GBR", "FRA", "DEU")) %>%
        #set everything to FALSE
        mutate(across(
          all_of(unique_cat),
          ~FALSE
        ),
        iso3c = c("PRK", "ERI", "ESH"),
        `Whole Virus` = c(TRUE, FALSE, FALSE),
        `Johnson&Johnson` = c(FALSE, TRUE, TRUE)
        )
    )


  #get the dominant type now
  dom <- readRDS("dominant_vaccines.rds") %>%
    mutate(iso3c = countrycode(country, "country.name", "iso3c")) %>%
    select(iso3c, dominant)
  df_dom <- df %>%
    left_join(dom, by = "iso3c") %>%
    mutate(
      dominant = if_else(
        iso3c == "CHN",
        "Sinovac",
        dominant
      )
    ) %>%
    rowwise() %>%
    mutate(
      dominant = simplify_vaccine(dominant)
    )
  #we assume all other nas are unknown
  df_dom
}
simplify_vaccine <- function(vaccine_type){
  #based upon https://covid19.trackvaccines.org/vaccines/
  if(vaccine_type %in% c("Pfizer/BioNTech", "Moderna")) {
    "mRNA"
  } else if (vaccine_type %in% c("Abdala", "Novavax", "ZF2001", "Soberana02", "Soberana Plus", "Razi Cov Pars", "SpikoGen", "EpiVacCorona", "Medigen", "Abdala (Subunit like novavax)", "Anhui ZL (Recombinant Platform)")) {
    "Subunit"
  } else if (vaccine_type %in% c("Oxford/AstraZeneca", "Sputnik V", "Sputnik Light", "CanSino", "Covishield")) {
    "Adenovirus"
  } else if (vaccine_type %in% c("Sinopharm/Wuhan", "Sinopharm/Beijing", "Sinovac", "Covaxin", "COVIran Barekat", "FAKHRAVAC", "QazVac", "Turkovac")) {
    "Whole Virus"
  } else {
    vaccine_type
  }
}
