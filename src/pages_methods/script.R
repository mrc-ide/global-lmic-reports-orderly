# Git and Orderly tagged footer

# orderly_id <- tryCatch(orderly::orderly_run_info()$id,
#                        error = function(e) "<id>") # bury this in the html, docx

# git <- system("git rev-parse --short HEAD", intern = TRUE)
# git_url <- "https://github.com/mrc-ide/global-lmic-reports-orderly/tree/"
# footer_rl <- readLines("footer.html")
# footer_rl[3] <- gsub("</p>",paste0(" | ",
#                                    "Git SHA: <a href=\"",git_url,git,"\">", git, "</a></p>"),
#                      footer_rl[3])
# writeLines(footer_rl, "footer.html")

rmarkdown::render("parameters.Rmd", output_format = c("html_document"),
                  output_options = list(pandoc_args = c("--metadata=title:COVID-19 Methods and Parameters")))
