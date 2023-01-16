get_platforms <- function(iso3cs, vacc_types, vdm, who_vacc, who_vacc_meta){
  #also get the ones scrape from OWID + other sources (i.e. UNICEF Dashboard)
  unicef <- readRDS("dominant_vaccines.Rds") %>%
    mutate(iso3c = countrycode(country, "country.name", "iso3c")) %>%
    select(iso3c, dominant)
  owid_2 <- get_owid_types(iso3cs)
  #create one large dataset
  all_df <- vacc_types %>%
    unnest(vaccine_types) %>%
    rbind(
      vdm %>%
        transmute(iso3c = iso3c, vaccine_types = vaccine) %>%
        unique()
    ) %>%
    rbind(
      who_vacc %>%
        transmute(iso3c = iso3c, vaccine_types = VACCINES_USED) %>%
        unnest(vaccine_types)
    ) %>%
    rbind(
      who_vacc_meta %>%
        transmute(iso3c = iso3c, vaccine_types = VACCINE_NAME)
    ) %>%
    rbind(
      unicef %>%
        transmute(iso3c = iso3c, vaccine_types = dominant)
    ) %>%
    rbind(
      owid_2 %>%
        transmute(
          iso3c = iso3c,
          vaccine_types = strsplit(vaccines, ", ")
        ) %>%
        unnest(vaccine_types)
    ) %>%
    unique()
  #detect any outlier
  not_matched <- all_df %>%
    transmute(platform = types_to_platforms(vaccine_types)) %>%
    pull(platform) %>%
    unique() %>%
    na.omit %>%
    setdiff(c("Adenovirus", "Whole Virus", "Single-Dose", "mRNA", "Subunit"))
  if(length(not_matched) > 0){
    stop(paste0("Failed to match: ", paste0(not_matched, collapse = ", "),
                " please modify types_to_platforms in R/get_platforms.R!"))
  }
  #now convert to platforms
  all_df %>%
    transmute(iso3c, platform = types_to_platforms(vaccine_types)) %>%
    unique() %>%
    filter(!is.na(platform)) %>%
    mutate(value = TRUE) %>%
    complete(iso3c = iso3cs, platform = unique(platform), fill = list(value = FALSE)) %>%
    pivot_wider(names_from = platform, values_from = value) %>%
  #if no vaccines we'll assume subunit
    mutate(
      Subunit = if_else(
        !Adenovirus & !mRNA & !`Single-Dose` & !`Whole Virus`,
        TRUE,
        Subunit
      )
    )


}

types_to_platforms <- function(vaccine_types){
  vaccine_types <- stringr::str_trim(vaccine_types)
  case_when(
    vaccine_types %in% c("Pfizer/BioNTech", "Moderna", "Moderna - Spikevax", "Pfizer BioNTech - Comirnaty", "Moderna - mRNA-1273",
                         "Zydus - ZyCov-D" ###REALLY DNA
                        ) ~ "mRNA",
    vaccine_types %in% c("Medicago", "Abdala", "Novavax", "ZF2001", "Soberana02", "Soberana Plus", "Razi Cov Pars", "SpikoGen", "EpiVacCorona", "Medigen", "Abdala (Subunit like novavax)", "Anhui ZL (Recombinant Platform)", "Corbevax", "Novavax-NUVAXOVID", "Anhui ZL - Zifivax", "CIGB - CIGB-66",
                         "Finlay - Soberana Plus", "Finlay - Soberana-02", "Biological E - Corbevax", "SII - Covovax", "Novavax - Covavax", "SRCVB - EpiVacCorona", "Sanofi/GSK","SKYCovione") ~ "Subunit",
    vaccine_types %in% c("Oxford/AstraZeneca", "Sputnik V", "Sputnik Light", "Gamaleya - Gam-Covid-Vac", "Gamaleya - Sputnik-Light", "CanSino", "Covishield", "AstraZeneca - Vaxzevria", "CanSino - Convidecia", "SII - Covishield", "AstraZeneca - AZD1222", "Gamaleya - Sputnik V",
                         "Shenzhen - LV-SMENP-DC" ###NOT REALLY AN ADENOVIRUS, but similar?
                         ) ~ "Adenovirus",
    vaccine_types %in% c("Sinopharm/Wuhan", "Sinopharm/Beijing", "Sinovac", "Covaxin", "Bharat - Covaxin", "COVIran Barekat", "FAKHRAVAC", "QazVac", "Turkovac", "KoviVac/Chumakov", "IMBCAMS", "KCONVAC", "Beijing CNBG - BBIBP-CorV", "Sinovac - CoronaVac", "Wuhan CNBG - Inactivated", "Chumakov - Covi-Vac", "IMB - Covidful", "Shifa - COVIran Barakat", "RIBSP - QazVac", "Julphar - Hayat-Vax", "Valneva", "Valneva - VLA2001") ~ "Whole Virus",
    vaccine_types %in% c("Janssen - Ad26.COV 2-S", "Johnson&Johnson") ~ "Single-Dose",
    vaccine_types %in% c("Moderna – Spikevax Bivalent Original/Omicron BA.1", "Pfizer BioNTech - Comirnaty Bivalent Original/Omicron BA.1", "Pfizer BioNTech - Comirnaty Bivalent Original/Omicron BA.4/BA.5", "Moderna – Spikevax Bivalent Original/Omicron  - Generic", "Pfizer BioNTech - Comirnaty Bivalent Original/Omicron - Generic") ~ "mRNA", #for now but need to add new VEs #"mRNA - Omicron",
    vaccine_types %in% c("Unknown Vaccine") ~ as.character(NA),
    TRUE ~ vaccine_types
  )
}

get_owid_types <- function(iso3cs){
  url <- xml2::read_html("https://github.com/owid/covid-19-data/tree/master/public/data/vaccinations/country_data")
  ctries <- rvest::html_text(rvest::html_nodes(url, ".js-navigation-open.Link--primary"))
  links <- paste0("https://raw.githubusercontent.com/owid/covid-19-data/master/public/data/vaccinations/country_data/", ctries)
  links <- urltools::url_encode(links)
  dat <- lapply(links, read.csv)
  data.frame("country" = gsub(".csv","", ctries, fixed = TRUE),
                   "vaccines" = unlist(lapply(dat, function(x){tail(x$vaccine,1)}))) %>%
    mutate(
      country = case_when(
        country %in% c("England", "Northern Ireland", "Scotland", "Wales") ~ "United Kingdom",
        country == "Timor" ~ "Timor Leste",
        TRUE ~ country
      ),
      country =
        countrycode(country, origin = "country.name", destination = "iso3c",
                    custom_match = c(Kosovo = "XKX"))
    ) %>%
    filter(country %in% iso3cs) %>%
    rename(iso3c = country) %>%
    group_by(iso3c) %>%
    summarise(vaccines = paste0(vaccines, collapse = ", "))
}
