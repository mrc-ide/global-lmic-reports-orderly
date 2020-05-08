#!/usr/bin/env Rscript

remove_old_reports <- function(date = NULL, keep = 3) {
  
  cwd <- getwd()
  wd <- "/home/oj/net/lmic/home/lmic/global-lmic-reports-orderly/"
  setwd(wd)
  country_root <- list.files("gh-pages", full.names = TRUE)
  country_root <- country_root[dir.exists(country_root)]
  country_root <- country_root[-grep("data", country_root)]
  
  for(i in seq_along(country_root)) {
    message(i)
  pages <- list.files(country_root[i])
  pages_full <- list.files(country_root[i], full.names = TRUE)
  to_clean <- which(as.Date(pages) < as.Date(date)-keep)
  
  files <- list.files(pages_full[to_clean], full.names = TRUE)
  unlink(files[-grep("json$|pdf$", files)], recursive = TRUE, force = TRUE)
  }
  setwd(cwd)
}


if (!interactive()) {
  usage <- "Usage:\n./remove_old_reports.R [<date>]"
  args <- docopt::docopt(usage)
  remove_old_reports(args$date)
}
