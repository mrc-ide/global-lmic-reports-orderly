#this repo
repo <- here::here()
#my one drive folder
destination <- file.path(
  repo, "gh-fits",
  ifelse(excess_mortality, "excess", "standard")
)
#get the fits
message("Gathering Fits")
fits <- squire.page::get_fits(repo = repo, date = date, iso3cs = NULL, excess = excess_mortality)

#temp
fits <- fits[!names(fits) %in% c("BRA", "IRQ", "KEN", "MEX", "MWI", "MYS", "PER", "ROU", "RUS",
                                 "SRB", "TUR", "UKR", "VEN", "URY", "ARE", "AUT", "BHR", "CHE",
                                 "EST", "HUN", "ISR", "POL", "PRT", "SVN", "SWE", "USA", "TCD",
                                 "DEU", "ITA", "MDV")]

#upload to folder replacing existing files
dir.create(destination, recursive = TRUE, showWarnings = FALSE)
for (iso in names(fits)) {
  saveRDS(fits[[iso]],
          file.path(
            destination,
            paste0(iso, ".Rds")
          ))
}
