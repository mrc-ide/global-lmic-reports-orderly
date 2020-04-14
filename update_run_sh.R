
update_run_sh <- function() {
rl <- readLines(file.path(here::here(),"run.sh"))

currently <- strsplit(grep("iso3c",rl,value=TRUE),"") %>% 
  lapply(function(x) {paste0(tail(x,3),collapse="")}) %>% 
  unlist

ecdc <- readRDS(file.path(here::here(),"archive","ecdc",tail(list.files(file.path(here::here(),"archive/ecdc/")),1),"ecdc_all.rds"))

with_deaths <- ecdc$countryterritoryCode[ecdc$deaths>0] %>% unique()

comment <- grep(paste0(currently[which(!currently %in% with_deaths)],collapse="|"),rl)

comment_func <- function(x, comment = TRUE) {
  
  if(substr(x,1,1) == "#") {
  x_comment <- x  
  x_no <- substr(x,2,nchar(x))
  } else {
    x_no <- x
    x_comment <- paste0("#", x, collapse = "")
  }
  
  if(comment) {
    return(x_comment)
  } else  {
    return(x_no)
  }
  
}

lines <- grep("iso3c",rl)
rl[lines] <- sapply(lines, function(x) {
  comment_func(rl[x], x %in% comment)
}) 

writeLines(rl, file.path(here::here(),"run.sh"))

}

update_for_dummy <- function(x){
  
  rl <- readLines(file.path(here::here(),"run.sh"))
  
  currently <- strsplit(grep("iso3c",rl,value=TRUE),"") %>% 
    lapply(function(x) {paste0(tail(x,3),collapse="")}) %>% 
    unlist
  
  ecdc <- readRDS(file.path(here::here(),"archive","ecdc",tail(list.files(file.path(here::here(),"archive/ecdc/")),1),"ecdc_all.rds"))
  
  with_deaths <- ecdc$countryterritoryCode[ecdc$deaths>0] %>% unique()
  
  comment <- grep(paste0(currently[which(!currently %in% with_deaths)],collapse="|"),rl)
  
  comment_func <- function(x, comment = TRUE) {
    
    if(substr(x,1,1) == "#") {
      x_comment <- x  
      x_no <- substr(x,2,nchar(x))
    } else {
      x_no <- x
      x_comment <- paste0("#", x, collapse = "")
    }
    
    if(comment) {
      return(x_comment)
    } else  {
      return(x_no)
    }
    
  }
  
  lines <- grep("iso3c",rl)
  rl[lines] <- sapply(lines, function(x) {
    comment_func(rl[x], x > 8)
  }) 
  
  writeLines(rl, file.path(here::here(),"run.sh"))
  
  
}

if(!interactive()) {
  update_run_sh()
}


if(!interactive()) {
  update_for_dummy()
}