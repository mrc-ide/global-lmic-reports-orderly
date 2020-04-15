yml_path <- file.path(here::here(),"src/lmic_reports/_site.yml")
yml <- yaml::read_yaml(yml_path)

rl <- readLines(file.path(here::here(),"run.sh"))
currently <- strsplit(grep("iso3c",rl,value=TRUE),"") %>%
  lapply(function(x) {paste0(tail(x,3),collapse="")}) %>%
  unlist

continents <- countrycode::countrycode(currently, origin = 'iso3c', destination = 'continent')

yml$navbar$left[[2]]$menu[[1]] <- list()
yml$navbar$left[[2]]$menu[[1]]$text <- "Africa"
yml$navbar$left[[2]]$menu[[2]] <- list()
yml$navbar$left[[2]]$menu[[2]]$text <- "Americas"
yml$navbar$left[[2]]$menu[[3]] <- list()
yml$navbar$left[[2]]$menu[[3]]$text <- "Asia"
yml$navbar$left[[2]]$menu[[4]] <- list()
yml$navbar$left[[2]]$menu[[4]]$text <- "Europe"
yml$navbar$left[[2]]$menu[[5]] <- list()
yml$navbar$left[[2]]$menu[[5]]$text <- "Oceania"

yml$navbar$left[[2]]$menu[[1]]$menu <- lapply(currently[continents == "Africa"],function(x){
  list("text" = unique(squire::population$country[squire::population$iso3c==x]),
       "href" = paste0("https://mrc-ide.github.io/global-lmic-reports/", x))
})
yml$navbar$left[[2]]$menu[[2]]$menu <- lapply(currently[continents == "Americas"],function(x){
  list("text" = unique(squire::population$country[squire::population$iso3c==x]),
       "href" = paste0("https://mrc-ide.github.io/global-lmic-reports/", x))
})
yml$navbar$left[[2]]$menu[[3]]$menu <- lapply(currently[continents == "Asia"],function(x){
  list("text" = unique(squire::population$country[squire::population$iso3c==x]),
       "href" = paste0("https://mrc-ide.github.io/global-lmic-reports/", x))
})
yml$navbar$left[[2]]$menu[[4]]$menu <- lapply(currently[continents == "Europe"],function(x){
  list("text" = unique(squire::population$country[squire::population$iso3c==x]),
       "href" = paste0("https://mrc-ide.github.io/global-lmic-reports/", x))
})
yml$navbar$left[[2]]$menu[[5]]$menu <- lapply(currently[continents == "Oceania"],function(x){
  list("text" = unique(squire::population$country[squire::population$iso3c==x]),
       "href" = paste0("https://mrc-ide.github.io/global-lmic-reports/", x))
})

yaml::write_yaml(yml, yml_path)