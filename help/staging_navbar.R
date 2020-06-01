wd <- "/home/oj/net/lmic_new/datadrive/lmic/testing/"
setwd(wd)

country_root <- list.files("gh-pages", full.names = TRUE)
sapply(grep("html", country_root, value = TRUE),
       xfun::gsub_file,"global-lmic-reports","global-lmic-reports-staging")
country_root <- country_root[dir.exists(country_root)]
country_root <- country_root[-grep("data", country_root)]
done <- sapply(file.path(country_root, "index.html"),
               xfun::gsub_file, "global-lmic-reports","global-lmic-reports-staging")
