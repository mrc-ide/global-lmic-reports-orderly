orderly_id <- tryCatch(orderly::orderly_run_info()$id,
                       error = function(e) "<id>") # bury this in the html, docx

rmarkdown::render("404.Rmd", output_format = c("html_document"))