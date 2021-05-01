
staging_navbar <- function() {
  
  country_root <- list.files("gh-pages", full.names = TRUE)
  sapply(grep("html", country_root, value = TRUE),
         xfun::gsub_file,"global-lmic-reports","global-lmic-reports-staging")
  country_root <- country_root[dir.exists(country_root)]
  sapply(country_root, function(x){list.files(x)})
  country_root <- country_root[-grep("data", country_root)]
  done <- sapply(file.path(country_root, "index.html"), function(x) {
    if(file.exists(x)) {
      xfun::gsub_file(x,"global-lmic-reports","global-lmic-reports-staging")
    }
  })
  
}


if (!interactive()) {
  usage <- "Usage:\n./staging_navbar.R"
  args <- docopt::docopt(usage)
  staging_navbar()
}