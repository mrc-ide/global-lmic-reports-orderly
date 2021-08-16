# copy data
file_copy <- function(from, to) {
  ok <- file.copy(from, to, overwrite = TRUE, recursive = TRUE)
  if (any(!ok)) {
    stop("There was an error copying files")
  }
}

file_copy_transfer <- function(x) {
  
  x2 <- gsub("testing", "global-lmic-reports-orderly", x)
  file_copy(x, x2)
  
}

file_copy_dir <- function(x) {
  
  x2 <- gsub("testing", "global-lmic-reports-orderly", x)
  dir.create(x2, FALSE, TRUE)
  to_copy <- list.files(x)
  to_copy <- file.path(x, to_copy)
  fz <- file.size(to_copy)
  to_copy <- to_copy[fz>0]
  file_copy(to_copy, x2)
  
  # also to copy
  also <- file.path(dirname(x), grep(".", list.files(dirname(x)), fixed = TRUE, value = TRUE))
  file_copy(also, dirname(x2))
}

roots <- list.files("/home/oj/net/lmic_new/datadrive/lmic/testing/gh-pages", full.names = TRUE)
roots <- roots[-grep(".", roots, fixed = TRUE)]
roots <- roots[-grep("data$", roots)]
roots <- file.path(roots, "2021-01-18")

for(i in 13:197) {
  
  message(i)
  file_copy_dir(roots[i])
  
}

file_copy_transfer("/home/oj/net/lmic_new/datadrive/lmic/testing/gh-pages/combined_reports.pdf")
file_copy_transfer("/home/oj/net/lmic_new/datadrive/lmic/testing/gh-pages/404.html")
file_copy_transfer("/home/oj/net/lmic_new/datadrive/lmic/testing/gh-pages/FAQ.html")
file_copy_transfer("/home/oj/net/lmic_new/datadrive/lmic/testing/gh-pages/index.html")
file_copy_transfer("/home/oj/net/lmic_new/datadrive/lmic/testing/gh-pages/News.html")
file_copy_transfer("/home/oj/net/lmic_new/datadrive/lmic/testing/gh-pages/parameters.html")

file_copy_transfer("/home/oj/net/lmic_new/datadrive/lmic/testing/gh-pages/data/2021-01-18_v7.csv.zip")
