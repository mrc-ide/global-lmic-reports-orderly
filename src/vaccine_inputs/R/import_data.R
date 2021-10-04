import_vaccine_agreements <- function(){
  vacc_loc_url <- "https://raw.githubusercontent.com/owid/covid-19-data/master/public/data/vaccinations/locations.csv"
  vacc_loc <- readr::read_csv(vacc_loc_url)
  vacc_loc$vaccine_types <- strsplit(vacc_loc$vaccines, ", ")
  dplyr::rename(vacc_loc, iso3c = iso_code)
}

import_vaccine_doses_by_manufacturer <- function(){
  vacc_by_type_url <- "https://raw.githubusercontent.com/owid/covid-19-data/master/public/data/vaccinations/vaccinations-by-manufacturer.csv"
  vacc_by_type <- readr::read_csv(vacc_by_type_url)
  dplyr::mutate(
    dplyr::filter(vacc_by_type, location != "European Union"),
                iso3c = countrycode::countrycode(
                  location, "country.name.en", "iso3c"
                )
    )
}

import_who_vaccination <- function(){
  who_vacc_url <- "https://covid19.who.int/who-data/vaccination-data.csv"
  who_vacc <- readr::read_csv(who_vacc_url)
  dplyr::mutate(
    dplyr::rename(who_vacc, iso3c = ISO3),
    VACCINES_USED = strsplit(VACCINES_USED, ",")
  )
}

import_who_vaccination_meta <- function(){
  who_vacc_meta_url <- "https://covid19.who.int/who-data/vaccination-metadata.csv"
  who_vacc_meta <- readr::read_csv(who_vacc_meta_url)
  dplyr::mutate(
    dplyr::rename(who_vacc_meta,
                  iso3c = ISO3),
         START_DATE = dplyr::if_else(
           START_DATE == "",
           as.Date(NA),
           START_DATE
         ))
}

import_owid <- function(){
  owid_url <- "https://covid.ourworldindata.org/data/owid-covid-data.csv"
  owid <- readr::read_csv(owid_url)
  dplyr::filter(
    dplyr::select(
    dplyr::rename(owid, iso3c = iso_code),
    c(iso3c, date, tidyselect::contains("vacc"))
  ),
  nchar(iso3c) == 3
  )
}
