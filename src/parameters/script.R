orderly_id <- tryCatch(orderly::orderly_run_info()$id,
                       error = function(e) "<id>") # bury this in the html, docx

rmarkdown::render("parameters.Rmd", output_format = c("html_document"),
                  output_options = list(pandoc_args = c("--metadata=title:\"COVID-19 Methods and Parameters\"")))
