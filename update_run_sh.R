#!/usr/bin/env Rscript

update_run_sh <- function(date, HICs = FALSE) {
rl <- readLines(file.path(here::here(),"countries"))

currently <- seq_along(rl)[-grep("#", rl)]
currently_iso <- rl[currently]
not <- seq_along(rl)[grepl("#", rl) & nchar(rl)<6]
not <- not[-1]
not_iso <- gsub("# ", "", rl[not])

db <- orderly::orderly_db("destination")
if (is.null(date)) {
  date <- as.character(Sys.Date())
}

## First find the id corresponding to the ecdc report with data.  If
## there are more than one, it's not totally clear what you want to
## do as you might want to take the earliest or the latest.
## Probably we want to take *all* and do the join over that, which
## is easy enough to do if you replace the '= $1' and replace with
## 'IN (%s)' and interpolate 'paste(sprintf('"%s"', id), collapse = ", ")'
sql <- 'SELECT report_version.id
            FROM report_version
            JOIN parameters
              ON parameters.report_version = report_version.id
           WHERE report_version.report = "ecdc"
             AND parameters.value = $1'
id <- DBI::dbGetQuery(db, sql, date)$id
if (length(id) == 0L) {
  stop(sprintf("No 'ecdc' report for '%s'", as.character(date)))
} else if (length(id) > 1) {
  message(sprintf("Multiple 'ecdc' reports for '%s'", as.character(date)))
}

ecdc <- readRDS(paste0(here::here(),"/archive/ecdc/",tail(id,1),"/ecdc_all.rds"))
with_deaths <- unique(ecdc$countryterritoryCode[ecdc$deaths>0])

# any to change
to_uncomment <- which(not_iso %in% with_deaths)
to_comment <- which(!currently_iso %in% with_deaths)

# do we remove the extra HICs
if(HICs) {
  also <- rl[(grep("Other HICs", rl)+1):length(rl)]
  to_comment <- unique(c(to_comment, which(currently_iso %in% also)))
}

# change them 
if(length(to_uncomment) > 0) {
rl[not[to_uncomment]] <- not_iso[to_uncomment]
}

if(length(to_comment) > 0) {
rl[currently[to_comment]] <- paste0("# ", currently_iso[to_comment])
}

writeLines(rl, file.path(here::here(),"countries"))

}

if(!interactive()) {
  usage <- "Usage:\n./update_run_sh.R [<date>]"
  args <- docopt::docopt(usage)
  update_run_sh(args$date)
}
