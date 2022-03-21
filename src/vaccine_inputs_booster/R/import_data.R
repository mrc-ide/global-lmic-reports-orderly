import_vaccine_agreements <- function(){

  vacc_loc <- readRDS("vaccine_agreements.rds")
  dplyr::rename(vacc_loc, iso3c = iso_code)
}

import_vaccine_doses_by_manufacturer <- function(){
  vacc_by_type <- readRDS("vaccine_doses_by_manufacturer.rds")
  dplyr::filter(
    dplyr::rename(vacc_by_type, iso3c = countryterritoryCode),
    !is.na(iso3c)
  )
}

import_who_vaccination <- function(){
  who_vacc <- readRDS("who_vacc.rds")
  dplyr::mutate(
    dplyr::rename(who_vacc, iso3c = countryterritoryCode),
    VACCINES_USED = strsplit(VACCINES_USED, ",")
  )
}

import_who_vaccination_meta <- function(){
  who_vacc_meta <- readRDS("who_vacc_meta.rds")
    dplyr::rename(who_vacc_meta,
                  iso3c = countryterritoryCode)
}

import_owid <- function(){
  owid <- readRDS("owid.rds")
  dplyr::filter(
    dplyr::select(
    dplyr::rename(owid, iso3c = countryterritoryCode),
    c(iso3c, date, tidyselect::contains("vacc"), tidyselect::contains("boost"))
  ),
  nchar(iso3c) == 3
  )
}
