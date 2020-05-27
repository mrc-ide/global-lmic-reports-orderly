
cumulative_deaths_plot <- function(country) {
  
  d <- readRDS("ecdc_all.rds")
  start <- 10
  
  d$Continent <- countrycode::countrycode(d$Region, origin = 'country.name', destination = 'continent')
  d$Continent[d$Region=="Eswatini"] <- "Africa"
  d$Continent[d$Region=="United State of America"] <- "Americas"
  d$Continent[d$Region=="Isle_of_Man"] <- "Europe"             
  d$Continent[d$Region=="Kosovo"] <- "Europe"                  
  d$Continent[d$Region=="Netherlands_Antilles"] <- "Americas"    
  d$Continent[d$Region=="Saint_Lucia"] <- "Americas"             
  d$Continent[d$Region=="South_Korea"] <- "Asia"             
  d$Continent[d$Region=="United_States_of_America"] <- "Americas"
  
  doubling <- function(double = 2, start = 10, xmax = 100) {
    
    x <- seq(0, xmax, 0.1)
    y <- start * 2^(x/double) 
    return(data.frame(x= x, y = y, 
                      Doubling = paste0("Every ", double, " Days")))
  }
  
  df <- group_by(d, Region) %>% 
    arrange(dateRep) %>% 
    mutate(Cum_Deaths = cumsum(deaths),
           Cum_Cases = cumsum(cases))
  
  df$Region <- gsub("_" ," ", df$Region)
  
  df_deaths <- df %>% 
    filter(Cum_Deaths > start) %>% 
    mutate(day_since = seq_len(n())-1)
  
  doubling_lines_deaths <- do.call(rbind, lapply(c(2, 3, 5, 7), function(x){
    doubling(x, start = start, xmax = max(df_deaths$day_since))
  }))
  
  df_deaths_latest <- df_deaths[df_deaths$dateRep == max(df_deaths$dateRep),]
  continent <- unique(df$Continent[df$Region %in% country])
  
  gg_deaths <- ggplot(df_deaths[df_deaths$Region %in% country,], aes(x=day_since, y=Cum_Deaths, group = Region)) + 
    geom_line(data = doubling_lines_deaths, aes(x=x, y=y, linetype = Doubling), inherit.aes = FALSE, color = "black") +
    geom_line(data = df_deaths, show.legend = FALSE, color = "grey", alpha = 0.6) +
    geom_line(data = df_deaths[which(df_deaths$Region %in% country[1:7]),], mapping = aes(color = Continent)) +
    #geom_point(data = df_deaths[which(df_deaths$Region %in% country[1:7]),], mapping = aes(color = Continent)) +
    geom_point(data = df_deaths_latest[which(df_deaths_latest$Region %in% country), ], alpha = 0.5, show.legend = FALSE) + 
    ggrepel::geom_text_repel(data =  df_deaths_latest[which(df_deaths_latest$Region %in% country), ],
                             aes(label = Region), show.legend = FALSE, min.segment.length = 2,nudge_x = 1) + 
    scale_y_log10(limits=c(start, max(df_deaths$Cum_Deaths[df_deaths$Continent %in% continent]))) +
    xlim(limits=c(0, max(df_deaths$day_since[df_deaths$Continent %in% continent]))) +
    theme_bw() +
    ylab("Cumulative Deaths (Logarithmic Scale)") +
    xlab(paste("Days Since", start, "Deaths"))
  
  gg_deaths
  
}

cumulative_deaths_plot_continent <- function(continent) {
  
  if(!continent %in% c("Asia","Europe","Africa","Americas","Oceania")) {
    stop("continent not matched")
  }
  
  rl <- readLines("_navbar.html")
  lmics <- gsub("(.*reports/)(\\w\\w\\w)(\".*)","\\2",grep("reports/(\\w\\w\\w)\"",rl, value =TRUE))
  
  colors <- c("#003b73","#e4572e","#BB750D","#003844","#925e78")
  col <- colors[match(continent, c("Asia","Europe","Africa","Americas","Oceania"))]
  
  d <- readRDS("ecdc_all.rds")
  d$Region[d$Region=="Congo"] <- "Republic of Congo"
  d$Region[d$Region=="United_Republic_of_Tanzania"] <- "Tanzania"
  d$Region[d$Region=="CuraÃ§ao"] <- "Curacao"
  start <- 10
  
  suppressWarnings(d$Continent <- countrycode::countrycode(d$Region, origin = 'country.name', destination = 'continent'))
  d$Continent[d$Region=="Eswatini"] <- "Africa"
  d$Continent[d$Region=="United State of America"] <- "Americas"
  d$Continent[d$Region=="Isle_of_Man"] <- "Europe"             
  d$Continent[d$Region=="Kosovo"] <- "Europe"                  
  d$Continent[d$Region=="Netherlands_Antilles"] <- "Americas"    
  d$Continent[d$Region=="Saint_Lucia"] <- "Americas"             
  d$Continent[d$Region=="South_Korea"] <- "Asia"             
  d$Continent[d$Region=="United_States_of_America"] <- "Americas"
  
  doubling <- function(double = 2, start = 10, xmax = 100) {
    
    x <- seq(0, xmax, 0.1)
    y <- start * 2^(x/double) 
    return(data.frame(x= x, y = y, 
                      Doubling = paste0("Every ", double, " Days")))
  }
  
  df <- group_by(d, Region) %>% 
    arrange(dateRep) %>% 
    mutate(Cum_Deaths = cumsum(deaths),
           Cum_Cases = cumsum(cases))
  
  df$Region <- gsub("_" ," ", df$Region)
  
  df_deaths <- df %>% 
    filter(Cum_Deaths > start) %>% 
    mutate(day_since = seq_len(n())-1)
  
  doubling_lines_deaths <- do.call(rbind, lapply(c(2, 3, 5, 7), function(x){
    doubling(x, start = start, xmax = max(df_deaths$day_since))
  }))
  
  df_deaths_latest <- df_deaths[df_deaths$dateRep == max(df_deaths$dateRep),]
  
  
  gg_deaths <- ggplot(df_deaths, aes(x=day_since, y=Cum_Deaths, group = Region)) + 
    geom_line(data = doubling_lines_deaths, aes(x=x, y=y, linetype = Doubling), inherit.aes = FALSE, color = "black") +
    geom_line(show.legend = FALSE, color = "grey", alpha = 0.3) +
    geom_line(data = df_deaths[which(df_deaths$Continent %in% continent & df_deaths$countryterritoryCode %in% lmics),], color = col) +
    #geom_point(data = df_deaths[which(df_deaths$Region %in% country[1:7]),], mapping = aes(color = Continent)) +
    geom_point(data = df_deaths_latest[which(df_deaths_latest$Continent %in% continent & df_deaths_latest$countryterritoryCode %in% lmics), ], alpha = 0.5, show.legend = FALSE) + 
    ggrepel::geom_text_repel(data =  df_deaths_latest[which(df_deaths_latest$Continent %in% continent & df_deaths_latest$countryterritoryCode %in% lmics), ],
                             aes(label = Region), show.legend = FALSE, min.segment.length = 0.1,nudge_x = 1,nudge_y = -0.1) + 
    scale_y_log10(limits=c(start, max(df_deaths$Cum_Deaths[df_deaths$Continent %in% continent & df_deaths$countryterritoryCode %in% lmics])), labels = scales::comma) +
    xlim(limits=c(0, max(df_deaths$day_since[df_deaths$Continent %in% continent & df_deaths$countryterritoryCode %in% lmics]))) +
    theme_bw() +
    scale_linetype(name = "Doubling Time:") +
    ylab("Cumulative Deaths (Logarithmic Scale)") +
    xlab(paste("Days Since", start, "Deaths")) + 
    ggtitle(continent)
  
  gg_deaths
  
}

cumulative_deaths_plot_continent_projections <- function(continent, today, data, ecdc) {
  
  ## handle arugments coming in
  if(!continent %in% c("Asia","Europe","Africa","Americas","Oceania")) {
    stop("continent not matched")
  }
  today <- as.Date(today)
  
  # identify lmics
  rl <- readLines("_navbar.html")
  lmics <- gsub("(.*reports/)(\\w\\w\\w)(\".*)","\\2",grep("reports/(\\w\\w\\w)\"",rl, value =TRUE))
  
  # set up colors
  colors <- c("#003b73","#e4572e","#BB750D","#003844","#925e78")
  col <- colors[match(continent, c("Asia","Europe","Africa","Americas","Oceania"))]
  
  # create dataset
  slim <- data %>% 
    mutate(date = as.Date(.data$date)) %>% 
    filter(date > (today)) %>%
    filter(date < (today+28)) %>% 
    filter(scenario == "Maintain Status Quo") %>% 
    select(date, compartment, y_mean, y_025, y_975, country, iso3c) %>% 
    mutate(observed = FALSE) %>% 
    rename(y = y_mean)
  
  slim <- slim[,c("date", "y", "country", "iso3c", "observed", "compartment", "y_025", "y_975")]
  
  # handle ecdc
  names(ecdc)[names(ecdc) %in% c("dateRep","Region", "countryterritoryCode")] <- c("date", "country", "iso3c")
  ecdc <- ecdc %>% 
    mutate(date = as.Date(date),
           observed = TRUE,
           compartment = "deaths") %>% 
    rename(y = deaths) %>% 
    select(date, y, country, iso3c, observed, compartment) %>% 
    mutate(y_025 = NA, 
           y_975 = NA)
  
  
  df <- as.data.frame(do.call(rbind, list(as.data.frame(ecdc), slim)), stringsAsFactors = FALSE)
  
  d <- df
  d$country[d$country=="Congo"] <- "Republic of Congo"
  d$country[d$country=="United_Republic_of_Tanzania"] <- "Tanzania"
  d$country[d$country=="CuraÃ§ao"] <- "Curacao"
  start <- 10
  
  suppressWarnings(d$Continent <- countrycode::countrycode(d$country, origin = 'country.name', destination = 'continent'))
  d$Continent[d$country=="Eswatini"] <- "Africa"
  d$Continent[d$country=="United State of America"] <- "Americas"
  d$Continent[d$country=="Isle_of_Man"] <- "Europe"             
  d$Continent[d$country=="Kosovo"] <- "Europe"                  
  d$Continent[d$country=="Netherlands_Antilles"] <- "Americas"    
  d$Continent[d$country=="Saint_Lucia"] <- "Americas"             
  d$Continent[d$country=="South_Korea"] <- "Asia"             
  d$Continent[d$country=="United_States_of_America"] <- "Americas"
  
  doubling <- function(double = 2, start = 10, xmax = 100) {
    
    x <- seq(0, xmax, 0.1)
    y <- start * 2^(x/double) 
    return(data.frame(x= x, y = y, 
                      Doubling = paste0("Every ", double, " Days")))
  }
  
  d <- d[d$compartment=="deaths",]
  d$date <- as.Date(d$date)
  
  
  df <- group_by(d, iso3c) %>% 
    arrange(date) %>%  
    mutate(Cum_Deaths = cumsum(y))
  
  df$country <- gsub("_" ," ", df$country)
  df <- df[which(df$iso3c %in% unique(df$iso3c[which(df$Cum_Deaths>10 & df$observed)])), ]
  
  df_deaths <- df %>% 
    filter(Cum_Deaths > start) %>% 
    mutate(day_since = seq_len(n())-1)
  
  doubling_lines_deaths <- do.call(rbind, lapply(c(2, 3, 5, 7), function(x){
    doubling(x, start = start, xmax = max(df_deaths$day_since))
  }))
  
  df_deaths_latest <- df_deaths[df_deaths$date == max(df_deaths$date),]
  
  
  gg_deaths <- ggplot(df_deaths, aes(x=day_since, y=Cum_Deaths, group = country)) + 
    #geom_line(data = doubling_lines_deaths, aes(x=x, y=y, linetype = Doubling), inherit.aes = FALSE, color = "black") +
    geom_line(show.legend = FALSE, color = "grey", alpha = 0.3) +
    geom_line(data = df_deaths[which(df_deaths$Continent %in% continent & df_deaths$iso3c %in% lmics & df_deaths$observed),], 
              color = col) +
    geom_line(data = df_deaths[which(df_deaths$Continent %in% continent & df_deaths$iso3c %in% lmics),], 
              linetype = "dashed", color = col) +
    #geom_point(data = df_deaths[which(df_deaths$country %in% country[1:7]),], mapping = aes(color = Continent)) +
    geom_point(data = df_deaths_latest[which(df_deaths_latest$Continent %in% continent & df_deaths_latest$iso3c %in% lmics), ], 
               alpha = 0.5, show.legend = FALSE) + 
    ggrepel::geom_text_repel(data =  df_deaths_latest[which(df_deaths_latest$Continent %in% continent & df_deaths_latest$iso3c %in% lmics), ],
                             aes(label = country), show.legend = FALSE, min.segment.length = 0.1,nudge_x = 1,nudge_y = -0.1) + 
    scale_y_log10(limits=c(start, max(df_deaths$Cum_Deaths[df_deaths$Continent %in% continent & df_deaths$iso3c %in% lmics])), 
                  labels = scales::comma) +
    xlim(limits=c(0, max(df_deaths$day_since[df_deaths$Continent %in% continent & df_deaths$iso3c %in% lmics])+15)) +
    theme_bw() +
    #scale_linetype(name = "Doubling Time:") +
    ylab("Cumulative Deaths (Logarithmic Scale)") +
    xlab(paste("Days Since", start, "Deaths"))
  
  gg_deaths
  
}

full_firework_plot <- function() {
  
  data <- read.csv("all_data.csv", stringsAsFactors = FALSE)
  
  ecdc <- readRDS("ecdc_all.rds")
  
  plots <- lapply(c("Asia","Europe","Africa","Americas","Oceania"), 
                  cumulative_deaths_plot_continent_projections, 
                  today = date, 
                  data = data, 
                  ecdc = ecdc)
  plotted <- lapply(plots[1:4], function(x){x+theme(legend.position = "none")})
  leg <- cowplot::get_legend(plots[[1]] + theme(legend.position = "top"))
  main <- cowplot::plot_grid(plotlist = plotted[1:4], ncol = 2)
  get <- cowplot::plot_grid(main,leg,ncol=1,rel_heights = c(1, 0.05))
  get
}

full_plot <- function() {

plots <- lapply(c("Asia","Europe","Africa","Americas","Oceania"), cumulative_deaths_plot_continent)
plotted <- lapply(plots[1:4], function(x){x+theme(legend.position = "none")})
leg <- cowplot::get_legend(plots[[1]] + theme(legend.position = "top"))
main <- cowplot::plot_grid(plotlist = plotted[1:4], ncol = 2)
get <- cowplot::plot_grid(main,leg,ncol=1,rel_heights = c(1, 0.05))
get
}

# define facet_zoom2 function to use FacetZoom2 instead of FacetZoom
# (everything else is the same as facet_zoom)
facet_zoom2 <- function(x, y, xy, zoom.data, xlim = NULL, ylim = NULL, 
                        split = FALSE, horizontal = TRUE, zoom.size = 2, 
                        show.area = TRUE, shrink = TRUE) {
  x <- if (missing(x)) if (missing(xy)) NULL else lazyeval::lazy(xy) else lazyeval::lazy(x)
  y <- if (missing(y)) if (missing(xy)) NULL else lazyeval::lazy(xy) else lazyeval::lazy(y)
  zoom.data <- if (missing(zoom.data)) NULL else lazyeval::lazy(zoom.data)
  if (is.null(x) && is.null(y) && is.null(xlim) && is.null(ylim)) {
    stop("Either x- or y-zoom must be given", call. = FALSE)
  }
  if (!is.null(xlim)) x <- NULL
  if (!is.null(ylim)) y <- NULL
  ggproto(NULL, FacetZoom2,
          shrink = shrink,
          params = list(
            x = x, y = y, xlim = xlim, ylim = ylim, split = split, zoom.data = zoom.data,
            zoom.size = zoom.size, show.area = show.area,
            horizontal = horizontal
          )
  )
}

# define FacetZoom as a ggproto object that inherits from FacetZoom,
# with a modified draw_panels function. the compute_layout function references
# the version currently on GH, which is slightly different from the CRAN
# package version.
FacetZoom2 <- ggproto(
  "FacetZoom2",
  ggforce::FacetZoom,
  
  compute_layout = function(data, params) {
    layout <- rbind( # has both x & y dimension
      data.frame(name = 'orig', SCALE_X = 1L, SCALE_Y = 1L),
      data.frame(name = 'x', SCALE_X = 2L, SCALE_Y = 1L),
      data.frame(name = 'y', SCALE_X = 1L, SCALE_Y = 2L),
      data.frame(name = 'full', SCALE_X = 2L, SCALE_Y = 2L),
      data.frame(name = 'orig_true', SCALE_X = 1L, SCALE_Y = 1L),
      data.frame(name = 'zoom_true', SCALE_X = 1L, SCALE_Y = 1L)
    )
    if (is.null(params$y) && is.null(params$ylim)) { # no y dimension
      layout <- layout[c(1,2, 5:6),]
    } else if (is.null(params$x) && is.null(params$xlim)) { # no x dimension
      layout <- layout[c(1,3, 5:6),]
    }
    layout$PANEL <- seq_len(nrow(layout))
    layout
  },
  
  draw_panels = function(panels, layout, x_scales, y_scales, ranges, coord,
                         data, theme, params) {
    
    if (is.null(params$x) && is.null(params$xlim)) {
      params$horizontal <- TRUE
    } else if (is.null(params$y) && is.null(params$ylim)) {
      params$horizontal <- FALSE
    }
    if (is.null(theme[['zoom']])) {
      theme$zoom <- theme$strip.background
    }
    if (is.null(theme$zoom.x)) {
      theme$zoom.x <- theme$zoom
    }
    if (is.null(theme$zoom.y)) {
      theme$zoom.y <- theme$zoom
    }
    axes <- render_axes(ranges, ranges, coord, theme, FALSE)
    panelGrobs <- ggforce:::create_panels(panels, axes$x, axes$y)
    panelGrobs <- panelGrobs[seq_len(length(panelGrobs) - 2)]
    if ('full' %in% layout$name && !params$split) {
      panelGrobs <- panelGrobs[c(1, 4)]
    }
    
    # changed coordinates in indicator / lines to zoom from 
    # the opposite horizontal direction
    if ('y' %in% layout$name) {
      if (!inherits(theme$zoom.y, 'element_blank')) {
        zoom_prop <- scales::rescale(
          y_scales[[2]]$dimension(ggforce:::expansion(y_scales[[2]])),
          from = y_scales[[1]]$dimension(ggforce:::expansion(y_scales[[1]])))
        indicator <- polygonGrob(
          x = c(0, 0, 1, 1), # was x = c(1, 1, 0, 0), 
          y = c(zoom_prop, 1, 0), 
          gp = gpar(col = NA, fill = alpha(theme$zoom.y$fill, 0.5)))
        lines <- segmentsGrob(
          x0 = c(1, 1), x1 = c(0, 0), # was x0 = c(0, 0), x1 = c(1, 1)
          y0 = c(0, 1), y1 = zoom_prop,
          gp = gpar(col = theme$zoom.y$colour,
                    lty = theme$zoom.y$linetype,
                    lwd = theme$zoom.y$size,
                    lineend = 'round'))
        indicator_h <- grobTree(indicator, lines)
      } else {
        indicator_h <- zeroGrob()
      }
    }
    
    if ('x' %in% layout$name) {
      if (!inherits(theme$zoom.x, 'element_blank')) {
        zoom_prop <- scales::rescale(x_scales[[2]]$dimension(ggforce:::expansion(x_scales[[2]])),
                                     from = x_scales[[1]]$dimension(ggforce:::expansion(x_scales[[1]])))
        indicator <- polygonGrob(c(zoom_prop, 1, 0), c(1, 1, 0, 0), 
                                 gp = gpar(col = NA, fill = alpha(theme$zoom.x$fill, 0.5)))
        lines <- segmentsGrob(x0 = c(0, 1), y0 = c(0, 0), x1 = zoom_prop, y1 = c(1, 1), 
                              gp = gpar(col = theme$zoom.x$colour,
                                        lty = theme$zoom.x$linetype,
                                        lwd = theme$zoom.x$size,
                                        lineend = 'round'))
        indicator_v <- grobTree(indicator, lines)
      } else {
        indicator_v <- zeroGrob()
      }
    }
    
    if ('full' %in% layout$name && params$split) {
      space.x <- theme$panel.spacing.x
      if (is.null(space.x)) space.x <- theme$panel.spacing
      space.x <- unit(5 * as.numeric(convertUnit(space.x, 'cm')), 'cm')
      space.y <- theme$panel.spacing.y
      if (is.null(space.y)) space.y <- theme$panel.spacing
      space.y <- unit(5 * as.numeric(convertUnit(space.y, 'cm')), 'cm')
      
      # change horizontal order of panels from [zoom, original] to [original, zoom]
      # final <- gtable::gtable_add_cols(panelGrobs[[3]], space.x)
      # final <- cbind(final, panelGrobs[[1]], size = 'first')
      # final_tmp <- gtable::gtable_add_cols(panelGrobs[[4]], space.x)
      # final_tmp <- cbind(final_tmp, panelGrobs[[2]], size = 'first')
      final <- gtable::gtable_add_cols(panelGrobs[[1]], space.x)
      final <- cbind(final, panelGrobs[[3]], size = 'first')
      final_tmp <- gtable::gtable_add_cols(panelGrobs[[2]], space.x)
      final_tmp <- cbind(final_tmp, panelGrobs[[4]], size = 'first')
      
      final <- gtable::gtable_add_rows(final, space.y)
      final <- rbind(final, final_tmp, size = 'first')
      final <- gtable::gtable_add_grob(final, list(indicator_h, indicator_h),
                                       c(2, 6), 3, c(2, 6), 5,
                                       z = -Inf, name = "zoom-indicator")
      final <- gtable::gtable_add_grob(final, list(indicator_v, indicator_v), 
                                       3, c(2, 6), 5, 
                                       z = -Inf, name = "zoom-indicator")
      heights <- unit.c(
        unit(max_height(list(axes$x[[1]]$top, axes$x[[3]]$top)), 'cm'),
        unit(1, 'null'),
        unit(max_height(list(axes$x[[1]]$bottom, axes$x[[3]]$bottom)), 'cm'),
        space.y,
        unit(max_height(list(axes$x[[2]]$top, axes$x[[4]]$top)), 'cm'),
        unit(params$zoom.size, 'null'),
        unit(max_height(list(axes$x[[2]]$bottom, axes$x[[4]]$bottom)), 'cm')
      )
      
      # swop panel width specifications according to the new horizontal order
      widths <- unit.c(
        # unit(max_width(list(axes$y[[3]]$left, axes$y[[4]]$left)), 'cm'),
        # unit(params$zoom.size, 'null'),
        # unit(max_height(list(axes$y[[3]]$right, axes$y[[4]]$right)), 'cm'),
        # space.x,
        # unit(max_width(list(axes$y[[1]]$left, axes$y[[2]]$left)), 'cm'),
        # unit(1, 'null'),
        # unit(max_height(list(axes$y[[1]]$right, axes$y[[2]]$right)), 'cm')        
        unit(max_width(list(axes$y[[1]]$left, axes$y[[2]]$left)), 'cm'),
        unit(1, 'null'),
        unit(max_height(list(axes$y[[1]]$right, axes$y[[2]]$right)), 'cm'),
        space.x,
        unit(max_width(list(axes$y[[3]]$left, axes$y[[4]]$left)), 'cm'),
        unit(params$zoom.size, 'null'),
        unit(max_height(list(axes$y[[3]]$right, axes$y[[4]]$right)), 'cm')
        
      )
      final$heights <- heights
      final$widths <- widths
    } else {
      if (params$horizontal) {
        space <- theme$panel.spacing.x
        if (is.null(space)) space <- theme$panel.spacing
        space <- unit(5 * as.numeric(convertUnit(space, 'cm')), 'cm')
        heights <- unit.c(
          unit(max_height(list(axes$x[[1]]$top, axes$x[[2]]$top)), 'cm'),
          unit(1, 'null'),
          unit(max_height(list(axes$x[[1]]$bottom, axes$x[[2]]$bottom)), 'cm')
        )
        
        # change horizontal order of panels from [zoom, original] to [original, zoom]
        # first <- gtable::gtable_add_cols(panelGrobs[[2]], space)
        # first <- cbind(final, panelGrobs[[1]], size = 'first')
        final <- gtable::gtable_add_cols(panelGrobs[[1]], space) 
        final <- cbind(final, panelGrobs[[2]], size = "first") 
        
        final$heights <- heights
        
        # swop panel width specifications according to the new horizontal order
        # unit(c(params$zoom.size, 1), 'null')
        final$widths[panel_cols(final)$l] <- unit(c(1, params$zoom.size), 'null') 
        
        final <- gtable::gtable_add_grob(final, indicator_h, 2, 3, 2, 5, 
                                         z = -Inf, name = "zoom-indicator")
      } else {
        space <- theme$panel.spacing.y
        if (is.null(space)) space <- theme$panel.spacing
        space <- unit(5 * as.numeric(convertUnit(space, 'cm')), 'cm')
        widths <- unit.c(
          unit(max_width(list(axes$y[[1]]$left, axes$y[[2]]$left)), 'cm'),
          unit(1, 'null'),
          unit(max_height(list(axes$y[[1]]$right, axes$y[[2]]$right)), 'cm')
        )
        final <- gtable::gtable_add_rows(panelGrobs[[1]], space)
        final <- rbind(final, panelGrobs[[2]], size = 'first')
        final$widths <- widths
        final$heights[panel_rows(final)$t] <- unit(c(1, params$zoom.size), 'null')
        final <- gtable::gtable_add_grob(final, indicator_v, 3, 2, 5, 
                                         z = -Inf, name = "zoom-indicator")
      }
    }
    final
  }
)


summaries_cases_continets_plot <- function(summaries) {
  
  sub <- summaries[!summaries$variable %in% c("hospital_14","icu_14",
                                              "hospital_14_mit","icu_14_mit",
                                              "hospital_14_rev","icu_14_rev",
                                              "report_deaths","report_deaths"),]
  
  sub$continent[is.na(sub$continent)] <- "Africa"
  sub <- group_by(sub, continent,variable) %>% 
    summarise(value=sum(value))
  
  ggplot(sub, 
         aes(x = continent, y = value, color = variable, fill = variable)) + 
    geom_bar(stat="identity",position = position_dodge2(preserve = "single"), width = 0.4) + 
    scale_y_log10(labels = scales::comma,breaks = c(10*(10^(0:8))),limits=c(1,10^7)) + 
    scale_fill_manual("", labels = c("Estimated Infections", "Reported Infections"),
                      values = viridis::viridis(3)) +
    scale_color_manual("", labels = c("Estimated Infections", "Reported Infections"),
                       values = viridis::viridis(3)) + 
    theme_bw() +
    xlab("") +
    ylab("") + 
    theme(legend.key = element_rect(size = 5),
          legend.key.size = unit(2, 'lines')) + 
    coord_flip() 
  
}

summaries_cases_plot <- function(summaries) {
  
  sub <- summaries[!summaries$variable %in% c("hospital_14","icu_14", "hospital_14_mit","icu_14_mit","report_deaths"),]
  levels(sub$country) <- rev(levels(sub$country))
  ggplot(sub, 
         aes(x = country, y = value, color = variable, fill = variable)) + 
    geom_bar(stat="identity",position = position_dodge2(preserve = "single"), width = 0.4) + 
    scale_y_log10(labels = scales::comma) + 
    scale_fill_manual("", labels = c("Estimated Infections", "Reported Infections", "Reported Deaths"),
                      values = viridis::viridis(3)) +
    scale_color_manual("", labels = c("Estimated Infections", "Reported Infections", "Reported Deaths"),
                      values = viridis::viridis(3)) + 
    theme_bw() +
    xlab("") +
    ylab("") + 
    facet_wrap(~continent, scales = "free") +
    theme(legend.key = element_rect(size = 5),
      legend.key.size = unit(2, 'lines')) + 
    coord_flip() 
  
}

summaries_forecasts_plot <- function(summaries) {
  
  summaries <- mutate(summaries, country = fct_rev(country)) %>% 
    mutate(value = ceiling(value)) %>% 
    mutate(group_large = "Africa & Europe") 
  
  summaries$group_large[summaries$continent %in% c("Asia", "Americas")] <- "Asia & Americas"
  
  gg <- ggplot(summaries[summaries$variable %in% c("hospital_14","icu_14"),], 
         aes(x = country, y = value, color = variable, fill = variable)) + 
    geom_bar(stat="identity",position = "dodge", width = 0.5) + 
    scale_y_continuous(labels = scales::comma) + 
    scale_fill_manual("", labels = c("Estimated Hospital Beds\nNeeded in 14 days",
                                     "Estimated ICU Beds\nNeeded in 14 days"),
                      values = viridis::viridis(3)) +
    scale_color_manual("", labels = c("Estimated Hospital Beds\nNeeded in 14 days",
                                      "Estimated ICU Beds\nNeeded in 14 days"),
                       values = viridis::viridis(3)) + 
    theme_bw() +
    xlab("") +
    ylab("") + 
    facet_grid(continent~., scales = "free", space = "free",) +
    theme(legend.key = element_rect(size = 5),
          legend.key.size = unit(2, 'lines')) + 
    scale_y_log10() +
    coord_flip()
  
  return(gg)
}

deaths_plot <- function(out, data, date = Sys.Date()) {
  
  o1 <- squire:::calibrate_output_parsing(
    out, 
    date_0 = as.Date(data$date[max(which(data$deaths == max(data$deaths)))])
  )
  
  gg_cases <- squire:::plot_calibration_healthcare_barplot(o1, data = data, forecast = 14) 
  gg_cases + geom_label(
    data = data.frame(x = c(as.Date(data$date[max(which(data$deaths == max(data$deaths)))]),date),
                      y = c(max(o1$y[o1$compartment == "deaths" & o1$date < (date+14)])*0.95,
                            max(o1$y[o1$compartment == "deaths" & o1$date < (date+14)])*0.75),
                      label=c("Calibration Date",as.character(date))), 
    aes(x=x, y=y, label=label), inherit.aes = FALSE)
  
  
}

healthcare_plot <- function(out, data) {
  
  o1 <- squire:::calibrate_output_parsing(
    out, 
    date_0 = as.Date(data$date[max(which(data$deaths == max(data$deaths)))])
  )
  
  cowplot::plot_grid(squire:::plot_calibration_healthcare_individual_barplot(df = o1, data = data, what = "hospital"),
                     squire:::plot_calibration_healthcare_individual_barplot(df = o1, data = data, what = "ICU"),
                     ncol=2)
  
}
