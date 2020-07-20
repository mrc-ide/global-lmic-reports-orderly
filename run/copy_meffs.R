#!/usr/bin/env Rscript

file_copy <- function(from, to) {
  ok <- file.copy(from, to, overwrite = TRUE, recursive = TRUE)
  if (any(!ok)) {
    stop("There was an error copying files")
  }
}

## Possibly useful:
## dir("gh-pages", pattern = "index\\.html$", recursive = TRUE)
copy_meffs <- function(date = NULL, what = "both", dic_only = TRUE, is_latest = TRUE) {
  
  dic_only <- as.logical(dic_only)
  if(!((what) %in% c("both", "3p", "4p"))) {
    stop("what is wrong")
  }
  
  #system("echo pre-DB")
  db <- orderly::orderly_db("destination")
  
  # first set up results dir
  target <- "gh-meffs/analysis/data/raw_data/server_results"
  #unlink("gh-results", recursive = TRUE, force = TRUE)
  
  # copy over orderly utilities
  dir.create(target, showWarnings = FALSE, recursive = TRUE)
  file_copy("orderly.sqlite", target)
  file_copy("orderly_config.yml", target)
  file_copy("global/", target)
  
  # now copy over ecdc, brt and reports
  if (is.null(date)) {
    date <- as.character(Sys.Date())
  }
  
  #system("echo pre-ecdc")
  ## First find the id corresponding to the ecdc report with data.  If
  ## there are more than one, it's not totally clear what you want to
  ## do as you might want to take the earliest or the latest.
  ## Probably we want to take *all* and do the join over that, which
  ## is easy enough to do if you replace the '= $1' and replace with
  ## 'IN (%s)' and interpolate 'paste(sprintf('"%s"', id), collapse = ", ")'
  sql <- 'SELECT report_version.id
            FROM report_version
            JOIN parameters
              ON parameters.report_version = report_version.id
           WHERE report_version.report = "ecdc"
             AND parameters.value = $1'
  id <- DBI::dbGetQuery(db, sql, date)$id
  if (length(id) == 0L) {
    stop(sprintf("No 'ecdc' report for '%s'", as.character(date)))
  } else if (length(id) > 1) {
    message(sprintf("Multiple 'ecdc' reports for '%s'", as.character(date)))
  }
  id_ecdc_max <- max(id)
  
  # copy lmic_reports_google
  src <- file.path("archive", "ecdc", id_ecdc_max)
  dest <- file.path(target,sprintf("%s/%s/%s", "archive", "ecdc",id_ecdc_max))
  worked <- vapply(dest, dir.create, logical(1), recursive = TRUE)
  worked <- mapply(file_copy, from = src, to = dirname(dest))
  
  
  ##  --------------------------------------------------------------------------
  ##  LMIC 3 Parameter ---------------------------------------------------------
  ##  --------------------------------------------------------------------------
  
  if(what %in% c("both", "3p")) {
    
    #system("echo pre-report")
    ## Then find all lmic_reports reports that use files from this ecdc
    ## report.  This is a bit awful and I might add direct link or a
    ## view to make this easier at some point.
    sql <- 'SELECT report_version.id, parameters.value as country
            FROM report_version_artefact
            JOIN file_artefact
              ON file_artefact.artefact = report_version_artefact.id
            JOIN depends
              ON depends.use = file_artefact.id
            JOIN report_version
              ON report_version.id = depends.report_version
            JOIN parameters
              ON parameters.report_version = report_version.id
           WHERE report_version_artefact.report_version IN (%s)
             AND report = "lmic_reports_google_pmcmc_no_decouple"
             AND parameters.name = "iso3c"
           ORDER BY country, report_version.id'
    sql <- sprintf(sql, paste(sprintf('"%s"', id), collapse = ", "))
    reports <- DBI::dbGetQuery(db, sql)
    
    if (any(duplicated(reports$country))) {
      keep <- tapply(seq_len(nrow(reports)), reports$country, max)
      reports <- reports[keep, ]
      rownames(reports) <- NULL
    }
    
    if(nrow(reports) > 0) {
      
      reports$date <- as.character(date)
      
      # copy lmic_reports_google key bits
      if (as.logical(dic_only)) {
        
        src <- file.path("archive", "lmic_reports_google_pmcmc_no_decouple", reports$id)
        dest <- file.path(target,sprintf("%s/%s/%s", "archive", "lmic_reports_google_pmcmc_no_decouple",reports$id))
        worked <- vapply(dest, dir.create, logical(1), recursive = TRUE)
        
        worked <- mapply(function(from, to){
          append <- c("grid_out.rds")
          file_copy(file.path(from, append), to)
        }, from = src, to = dest)
        
        for(i in seq_len(length(dest))) {
          out <- readRDS(file.path(dest[i], "grid_out.rds"))
          
          # sort deaths out in the output
          index <- squire:::odin_index(out$model)
          for(i in seq_along(out$parameters$population)) {
            collect <- vapply(1:out$parameters$replicates, function(j) {
              pos <- seq(i, length(index$D), by = length(out$parameters$population))
              pos <- index$D[pos]
              diff(out$output[,pos,j])
            }, FUN.VALUE = numeric(nt-1))
            out$output[1+seq_len(nt-1),index$delta_D[i],] <- collect
          }
          nt <-  nrow(out$output)
          deaths <- squire:::odin_sv(out$output[,index$delta_D,], replicates = dim(out$output)[3], nt, reduce_age = TRUE)
          df <- data.frame("t" = as.numeric(out$output[,index$time,]), 
                           "replicate" = as.numeric(mapply(rep, seq_len(out$parameters$replicates), nt)),
                           "compartment" = "deaths",
                           "y" = deaths)
          df$date <- as.Date(df$t + as.Date(date),
                              format = "%Y-%m-%d")
          out$output <- df
          
          # clear out what's not needed
          out$pmcmc_results$inputs$squire_model <- NULL
          for(j in seq_len(length(out$pmcmc_results$chains))) {
            out$pmcmc_results$chains[[j]]$covariance_matrix <- NULL
            out$pmcmc_results$chains[[j]]$scaling_factor <- NULL
            out$pmcmc_results$chains[[j]]$acceptances <- NULL
            out$pmcmc_results$chains[[j]]$results <- tail(out$pmcmc_results$chains[[j]]$results,5000)
          }
          saveRDS(out, file.path(dest[i], "grid_out.rds"))
        }
        
        
      } else {
        
        src <- file.path("archive", "lmic_reports_google_pmcmc_no_decouple", reports$id)
        dest <- file.path(target,sprintf("%s/%s/%s", "archive", "lmic_reports_google_pmcmc_no_decouple",reports$id))
        worked <- vapply(dest, dir.create, logical(1), recursive = TRUE)
        
        worked <- mapply(function(from, to){
          append <- c("fitting.pdf", "grid_out.rds", "projections.csv", "orderly_run.rds")
          file_copy(file.path(from, append), to)
        }, from = src, to = dest)
        
        
      }
      
    }
    
  }
  ##  --------------------------------------------------------------------------
  ##  LMIC PMCMC COPY ----------------------------------------------------------
  ##  --------------------------------------------------------------------------
  
  if(what %in% c("both", "4p")) {
    
    #system("echo pre-report")
    ## Then find all lmic_reports reports that use files from this ecdc
    ## report.  This is a bit awful and I might add direct link or a
    ## view to make this easier at some point.
    sql <- 'SELECT report_version.id, parameters.value as country
            FROM report_version_artefact
            JOIN file_artefact
              ON file_artefact.artefact = report_version_artefact.id
            JOIN depends
              ON depends.use = file_artefact.id
            JOIN report_version
              ON report_version.id = depends.report_version
            JOIN parameters
              ON parameters.report_version = report_version.id
           WHERE report_version_artefact.report_version IN (%s)
             AND report = "lmic_reports_google_pmcmc"
             AND parameters.name = "iso3c"
           ORDER BY country, report_version.id'
    sql <- sprintf(sql, paste(sprintf('"%s"', id), collapse = ", "))
    reports <- DBI::dbGetQuery(db, sql)
    
    if (any(duplicated(reports$country))) {
      keep <- tapply(seq_len(nrow(reports)), reports$country, max)
      reports <- reports[keep, ]
      rownames(reports) <- NULL
    }
    
    if(nrow(reports) > 0) {
      
      reports$date <- as.character(date)
      
      # copy lmic_reports_google key bits
      if (as.logical(dic_only)) {
        
        src <- file.path("archive", "lmic_reports_google_pmcmc", reports$id)
        dest <- file.path(target, sprintf("%s/%s/%s", "archive", "lmic_reports_google_pmcmc", reports$id))
        worked <- vapply(dest, dir.create, logical(1), recursive = TRUE)
        
        worked <- mapply(function(from, to){
          append <- c("grid_out.rds")
          file_copy(file.path(from, append), to)
        }, from = src, to = dest)
        
        for(i in seq_len(length(dest))) {
          out <- readRDS(file.path(dest[i], "grid_out.rds"))
          
          # sort deaths out in the output
          index <- squire:::odin_index(out$model)
          for(i in seq_along(out$parameters$population)) {
            collect <- vapply(1:out$parameters$replicates, function(j) {
              pos <- seq(i, length(index$D), by = length(out$parameters$population))
              pos <- index$D[pos]
              diff(out$output[,pos,j])
            }, FUN.VALUE = numeric(nt-1))
            out$output[1+seq_len(nt-1),index$delta_D[i],] <- collect
          }
          nt <-  nrow(out$output)
          deaths <- squire:::odin_sv(out$output[,index$delta_D,], replicates = dim(out$output)[3], nt, reduce_age = TRUE)
          df <- data.frame("t" = as.numeric(out$output[,index$time,]), 
                           "replicate" = as.numeric(mapply(rep, seq_len(out$parameters$replicates), nt)),
                           "compartment" = "deaths",
                           "y" = deaths)
          df$date <- as.Date(df$t + as.Date(date),
                             format = "%Y-%m-%d")
          out$output <- df
          
          # clear out what's not needed
          out$pmcmc_results$inputs$squire_model <- NULL
          for(j in seq_len(length(out$pmcmc_results$chains))) {
            out$pmcmc_results$chains[[j]]$covariance_matrix <- NULL
            out$pmcmc_results$chains[[j]]$scaling_factor <- NULL
            out$pmcmc_results$chains[[j]]$acceptances <- NULL
            out$pmcmc_results$chains[[j]]$results <- tail(out$pmcmc_results$chains[[j]]$results,5000)
          }
          saveRDS(out, file.path(dest[i], "grid_out.rds"))
        }
        
        
      } else {
        
        src <- file.path("archive", "lmic_reports_google_pmcmc", reports$id)
        dest <- file.path(target,sprintf("%s/%s/%s", "archive", "lmic_reports_google_pmcmc",reports$id))
        worked <- vapply(dest, dir.create, logical(1), recursive = TRUE)
        
        worked <- mapply(function(from, to){
          append <- c("fitting.pdf", "grid_out.rds", "projections.csv", "orderly_run.rds")
          file_copy(file.path(from, append), to)
        }, from = src, to = dest)
        
        
      }
      
    }
   
  }
  
  ##  --------------------------------------------------------------------------
  ## BRT COPY ------------------------------------------------------------------
  ##  --------------------------------------------------------------------------
  
  #system("echo pre-brt")
  sql <- 'SELECT report_version.id
            FROM report_version
            JOIN parameters
              ON parameters.report_version = report_version.id
           WHERE report_version.report = "brt_google_mobility"
             AND parameters.value = $1'
  sql <- sprintf(sql, paste(sprintf('"%s"', id), collapse = ", "))
  reports <- DBI::dbGetQuery(db, sql, date)
  brt_id_max <- max(reports$id)
  
  # copy lmic_reports_google
  src <- file.path("archive", "brt_google_mobility", brt_id_max)
  dest <- file.path(target,sprintf("%s/%s/%s", "archive", "brt_google_mobility",reports$id))
  worked <- vapply(dest, dir.create, logical(1), recursive = TRUE)
  worked <- mapply(file_copy, from = src, to = dirname(dest))
  
}


if (!interactive()) {
  usage <- "Usage:\n./copy_meffs.R [<date>] [<what>] [<dic_only>]"
  args <- docopt::docopt(usage)
  copy_meffs(args$date, args$what, args$dic_only)
}
