# inputs
date <- "2021-08-16"
short_run <- TRUE
full_scenarios <- FALSE
gibbs_sampling <- FALSE
n_mcmc <- 20
countries <- "countries"
parallel <- FALSE

message("*** ECDC data")
ecdc_id <- orderly::orderly_run("ecdc", parameters = list(date=date), echo = FALSE)
orderly::orderly_commit(ecdc_id)

message("*** Oxford GRT data")
oxford_id <- orderly::orderly_run("oxford_grt", parameters = list(date=date), echo = FALSE)
orderly::orderly_commit(oxford_id)

message("*** Google BRT data")
google_id <- orderly::orderly_run("brt_google_mobility", parameters = list(date=date, short_run = short_run), echo = FALSE)
orderly::orderly_commit(google_id)

message("*** Running country reports")

# get the isos
iso3cs <- grep('^[A-Z]{3}\\s*', readLines(countries), value = TRUE)

# make the orderly bundles to be run on the cluster
path_bundles <- file.path("L:/OJ/glodide/analysis/data/raw/", date)
dir.create(path_bundles, showWarnings = FALSE)
file.remove(list.files(path_bundles, full.names = TRUE))

# bundle these up
bundles <- lapply(
  iso3cs, function(x) {
    orderly::orderly_bundle_pack(
      path = path_bundles,
      name = "lmic_reports_vaccine",
      parameters = list(
        iso3c = x,
        date=date,
        short_run=short_run,
        parallel=parallel,
        full_scenarios=full_scenarios,
        gibbs_sampling=gibbs_sampling,
        n_mcmc=n_mcmc
      )
    )
  }
)

# pull these back into orderly
path_outs <- list.files(file.path("L:/OJ/glodide/analysis/data/derived", date), full.names = TRUE)
import <- lapply(path_outs, orderly::orderly_bundle_import)

# now we run this script to output out combined pdf of fits for us to manually check in the fits folder
system(
  command = "cmd.exe",
  input = paste(
    '"C:\\Program Files\\R\\R-4.0.2\\bin\\i386\\Rscript.exe" help/fits_for_checking.R',
    date,
    "lmic_reports_vaccine"
  )
)

# TODO:

# 1. Get gh-pages
# If gh-pages is not in the repo, i.e. it has just been called then you need to run following in terminal
# mkdir gh-pages
# git clone git@github.com:mrc-ide/global-lmic-reports.git gh-pages

# 2. The collation step for run_collate_vaccine
# Jump to the terminal and execute the bash (probably a way to do this from R but still working on it)
# alt + shift + m to go to terminal and then run
# ./run/run_collate_vaccine.sh 2021-08-16

# 3. Push the results to the repository if we are happy with them again from bash
# N.B. This is done in two steps as the git pack size if all done at once is too big
# VERSION=$(git rev-parse --short HEAD)
# git add data
# git commit --no-verify -m "Update data for version ${VERSION}"
# git push
# git add .
# git commit --no-verify -m "Update pages for version ${VERSION}"
# git push

# Alternatively, this is bundled in this one bash:
# ./scripts/publish_website.sh production

# Lastly we also have a staging site, so if you ever want to check that the pages etc look okay
# the set the remote to there and push to there
# N.B. The hyperlinks in the page build refer to the main site but the URLs are correct so you can check an individual page by going to the URL.
# git remote set-url origin git@github.com:mrc-ide/global-lmic-reports-staging.git






#
# if we are happy with these
