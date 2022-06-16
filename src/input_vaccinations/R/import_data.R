import_vaccine_agreements <- function(){
  read_csv("https://raw.githubusercontent.com/owid/covid-19-data/master/public/data/vaccinations/locations.csv") %>%
    rename(iso3c = iso_code) %>%
    mutate(vaccine_types = strsplit(vaccines, ", ")) %>%
    select(iso3c, vaccine_types)
}

import_vaccine_doses_by_manufacturer <- function(){
  read_csv("https://raw.githubusercontent.com/owid/covid-19-data/master/public/data/vaccinations/vaccinations-by-manufacturer.csv") %>%
    filter(location != "European Union") %>%
    mutate(iso3c = countrycode::countrycode(location, "country.name", "iso3c")) %>%
    select(!location)
}

import_who_vaccination <- function(){
  read_csv("https://covid19.who.int/who-data/vaccination-data.csv") %>%
    rename(iso3c = ISO3) %>%
    mutate(
      VACCINES_USED = strsplit(VACCINES_USED, ",")
    )
}

import_who_vaccination_meta <- function(){
  read_csv("https://covid19.who.int/who-data/vaccination-metadata.csv") %>%
    rename(iso3c = ISO3)
}

import_owid <- function(){
  read_csv("https://covid.ourworldindata.org/data/owid-covid-data.csv") %>%
    rename(iso3c = iso_code) %>%
    select(iso3c, date, tidyselect::contains("vacc"), tidyselect::contains("boost")) %>%
    filter(nchar(iso3c) == 3)
}
