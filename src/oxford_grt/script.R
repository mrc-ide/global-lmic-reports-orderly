#############

# get user provided date
date <- as.Date(date, "%Y-%m-%d")
if (is.na(date)) {
  stop("Date must be provided in ISO format (i.e., YYYY-MM-DD)")
}

# download
url <- "https://ocgptweb.azurewebsites.net/CSVDownload"
tf <- tempfile(fileext = )
download.file(url, tf)
d <- read.csv(tf)

# max date
max_date <- as.Date(as.character(max(d$Date)), "%Y%m%d")
if(max_date < date) {
  stop("Oxford Goverment Tracker Dataset Not Updated for Today")
}

names(d)[names(d) %in% c("S1_School.closing",
                         "S2_Workplace.closing",
                         "S3_Cancel.public.events",
                         "S6_Restrictions.on.internal.movement")] <- c("S1","S2","S3","S6")
vars <- names(d)[c(1,2,3,4,7,10,13,16,19,38,39)]

# cleaning



# indicator function
ind <- function(x) {
  as.integer(x %in% 1:2)
}

# overall movement reduction
s <- group_by(d, CountryCode, Date) %>% 
  summarise(S1 = S1,
            S2 = S2,
            S3 = S3, 
            S6 = S6,
            C_SW = 1-(((0.15/0.75*ind(S1)) + ((0.75*0.6/0.75)*ind(S2)))),
            C_GM = 1-(if(ind(S6)) {0.75} else {0.1*ind(S3)}),
            C = (C_SW + C_GM)/2)


## assumed starting
start <- "20200201"
early_start <- "20200115"

## different starts
earlys <- c("BRA","DOM","DZA","ECU","IDN","IND","MEX","PER","PHL","ROU","RUS","TUR", "COL","IRQ","JAM","LKA", "PAK", "PSE")

# summarise these per country to make an assumed start of 1st February.
tt <- lapply(unique(s$CountryCode), function(x){
  
  if(x %in% earlys) {
    correct <- early_start
  } else {
    correct <- start
  }
  
  empty <- data.frame(iso3c=x, 
                      date = correct,
                      S1 = 0,
                      S2 = 0,
                      S3 = 0,
                      S6 = 0,
                      C_SW = 1,
                      C_GM = 1,
                      C = 1,
                      tt_R0 = 0,
                      R0 = 3,
                      stringsAsFactors = FALSE)
  
  if(length(unique(s$C[s$CountryCode==x])) == 1) {
    return(empty)  
  } else {
    
    d <- s[s$CountryCode==x, ]
    d <- d[d$Date>=20200201,]
    d$tt_R0 <- as.numeric(as.Date(as.character(d$Date),"%Y%m%d")-as.Date(correct,"%Y%m%d"))
    d$R0 <- d$C * 3
    names(d)[1:2] <- c("iso3c","date")
    
    if(d$date[1] != correct) {
     d <- rbind(empty, as.data.frame(d)) 
    }
    
    return(as.data.frame(d))
  }
})
names(tt) <- unique(s$CountryCode)

## data cleaning 
today <- as.numeric(Sys.Date() - as.Date(start,"%Y%m%d"))
early_today <- as.numeric(Sys.Date() - as.Date(early_start,"%Y%m%d"))
tt2 <- tt
for(j in seq_along(tt2)) {
  
  m <- tt[[j]]
  
  if(as.character(m$iso3c[1]) %in% earlys) {
    real_today <- early_today
  } else {
    real_today <- today
  }
  
  
  lates <- unique(c(which(!m$tt_R0 <= (real_today-23))))
  if(length(lates) > 0) {
    for(i in seq_along(lates)) {
      if(is.na(m$S1[lates[i]])) {
        m$S1[lates[i]] <- m$S1[lates[i]-1]
      }
      if(is.na(m$S2[lates[i]])) {
        m$S2[lates[i]] <- m$S2[lates[i]-1]
      }
      if(is.na(m$S3[lates[i]])) {
        m$S3[lates[i]] <- m$S3[lates[i]-1]
      }
      if(is.na(m$S6[lates[i]])) {
        m$S6[lates[i]] <- m$S6[lates[i]-1]
      }
    }
  }
  tt2[[j]] <- m
}

# manual correction after reading
tt2$ALB$S6[tt2$ALB$date > c(20200312)] <- 2
tt2$ALB$S3[tt2$ALB$date > c(20200312)] <- 2
tt2$ALB$S1[tt2$ALB$date > c(20200312)] <- 2
tt2$TZA$S6[tt2$TZA$date %in% c(20200327)] <- 2
tt2$BMU$S6[tt2$BMU$date %in% c(20200323)] <- 0
tt2$NOR$S6[tt2$NOR$date %in% c(20200323)] <- 0
tt2$PHL$S3[tt2$PHL$date %in% c(20200323)] <- 0
tt2$PSE$S3[tt2$PSE$date %in% c(20200320,20200321)] <- 2

# anything NA now is a 0
for(i in seq_along(tt2)) {
  tt2[[i]][is.na(tt2[[i]])] <- 0
}

# add blanks 
nms <- unique(squire::population$iso3c)
tt_no <- lapply(nms[!nms %in% names(tt)], function(x) {
  return(data.frame(iso3c=x, 
                    date = start,
                    S1 = 0,
                    S2 = 0,
                    S3 = 0,
                    S6 = 0,
                    C_SW = 1,
                    C_GM = 1,
                    C = 1,
                    tt_R0 = 0,
                    R0 = 3,
                    stringsAsFactors = FALSE))  
})
names(tt_no) <- nms[!nms %in% names(tt)]

# group them and apply the change calculations again after the data cleaning
res <- append(tt2, tt_no)
res <- lapply(res, function(x) {
  x %>% group_by(iso3c, tt_R0, date) %>% 
  summarise(S1 = S1,
            S2 = S2,
            S3 = S3,
            S6 = S6,
            C_SW = 1-(((0.15/0.75*ind(S1)) + ((0.75*0.6/0.75)*ind(S2)))),
            C_GM = 1-(if(ind(S6)) {0.75} else {0.1*ind(S3)}),
            C = (C_SW + C_GM)/2,
            R0 = 3*C) %>% 
    as.data.frame(stringsAsFactors = FALSE)
})

# remove any dates in the future
res <- lapply(res,function(x){
  
  x[as.Date(as.character(x$date), "%Y%m%d") <= date,]
  
})

saveRDS(res, "oxford_grt.rds")

