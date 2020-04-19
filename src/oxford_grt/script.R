#############

# get user provided date
date <- as.Date(date, "%Y-%m-%d")
if (is.na(date)) {
  stop("Date must be provided in ISO format (i.e., YYYY-MM-DD)")
}

# download
url <- "https://ocgptweb.azurewebsites.net/CSVDownload"
tf <- tempfile()
download.file(url, tf)
d <- read.csv(tf)

# max date
max_date <- as.Date(as.character(max(d$Date)), "%Y%m%d")
if(max_date != date) {
  stop("Oxford Goverment Tracker Dataset Not Updated for Today")
}

# indicator function
ind <- function(x) {
  if(x %in% 1:2) {
    return(1)
  } else {
    return(0)
  }
}

# overall movement reduction
s <- group_by(d, CountryCode, Date) %>% 
  summarise(C_SW = 1-(((0.15*ind(S1_School.closing)) + ((0.75*0.6)*ind(S2_Workplace.closing)))),
            C_GM = 1-(if(ind(S6_Restrictions.on.internal.movement)) {0.75} else {0.1*ind(S3_Cancel.public.events)}),
            C_H = if(ind(S6_Restrictions.on.internal.movement)) { 
              1.2 
            } else {
              1 + 0.2*((0.15*ind(S1_School.closing)) + (0.6*ind(S2_Workplace.closing)))
            },
            C = (C_SW + C_GM + C_H)/3)

# summarise these per country
tt <- lapply(unique(s$CountryCode), function(x){
  
  if(length(unique(s$C[s$CountryCode==x]))) {
    return(data.frame("iso3c"=x, 
                      tt_R0 = 0, 
                      R0 = 3))  
  } else {
  r <- rle(s$C[s$CountryCode==x])
  d <- s$Date[s$CountryCode==x][cumsum(r$lengths)[-1]]
  c <- s$C[s$CountryCode==x][cumsum(r$lengths)[-1]]
  return(data.frame("iso3c"=x, 
             tt_R0 = as.numeric(as.Date(as.character(d),"%Y%m%d")-as.Date("20200201","%Y%m%d")), 
             R0 = c*3))
  }
})

# add blanks 
nms <- unique(squire::population$iso3c)
tt_no <- lapply(nms[!nms %in% names(tt)], function(x) {
  return(data.frame("iso3c"=x, 
                    tt_R0 = 0, 
                    R0 = 3))  
})
names(tt_no) <- nms[!nms %in% names(tt)]

res <- append(tt, tt_no)
saveRDS(res, "oxford_grt.rds")

