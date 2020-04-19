
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
  continent <- unique(df$Continent[df$Region == country])
  
  gg_deaths <- ggplot(df_deaths[which(df_deaths$Continent == continent), ], aes(x=day_since, y=Cum_Deaths, group = Region)) + 
    geom_line(data = doubling_lines_deaths, aes(x=x, y=y, linetype = Doubling),alpha=0.3, inherit.aes = FALSE) +
    geom_line(show.legend = FALSE, color = "grey", alpha = 0.6) +
    geom_line(data = df_deaths[which(df_deaths$Region == country), ], color = "red", lwd = 1) +
    geom_point(data = df_deaths[which(df_deaths$Region == country), ], color = "red") +
    geom_point(data = df_deaths_latest[which(df_deaths_latest$Continent == continent), ], alpha = 0.5, show.legend = FALSE) + 
    ggrepel::geom_text_repel(data =  df_deaths_latest[which(df_deaths_latest$Continent == continent), ],
                             aes(label = Region), show.legend = FALSE, min.segment.length = 1,nudge_x = 1) + 
    scale_y_log10(limits=c(start, max(df_deaths$Cum_Deaths[df_deaths$Continent == continent]))) +
    xlim(limits=c(0, max(df_deaths$day_since[df_deaths$Continent == continent]))) +
    theme_bw() +
    ylab("Cumulative Deaths (Logarithmic Scale)") +
    xlab(paste("Days Since", start, "Deaths"))
  
  gg_deaths
  
}


plotly_style <- function(country) {
  
  d <- readRDS("ecdc_all.rds")
  
  d$Continent <- countrycode::countrycode(d$Region, origin = 'country.name', destination = 'continent')
  d$Continent[d$Region=="Eswatini"] <- "Africa"
  d$Continent[d$Region=="United State of America"] <- "Americas"
  d$Continent[d$Region=="Isle_of_Man"] <- "Europe"             
  d$Continent[d$Region=="Kosovo"] <- "Europe"                  
  d$Continent[d$Region=="Netherlands_Antilles"] <- "Americas"    
  d$Continent[d$Region=="Saint_Lucia"] <- "Americas"             
  d$Continent[d$Region=="South_Korea"] <- "Asia"             
  d$Continent[d$Region=="United_States_of_America"] <- "Americas"
  d$Region <- gsub("_" ," ", d$Region)
  
  # to reduce need for user inputs, precalculate for set options
  d <- group_by(d, Region) %>% 
    arrange(dateRep) %>% 
    mutate(Cumulative_Deaths = cumsum(deaths),
           Cumulative_Cases = cumsum(cases),
           cc_10 = Cumulative_Cases > 10,
           cc_100 = Cumulative_Cases > 100,
           cc_1000 = Cumulative_Cases > 1000,
           cd_10 = Cumulative_Deaths > 10,
           cd_100 = Cumulative_Deaths > 100,
           cd_1000 = Cumulative_Deaths > 1000) %>%
    ungroup()
  
  d <- d %>% group_by(Region, cc_10) %>% 
    mutate(day_since_cc_10 = seq_len(n())-1) %>% ungroup %>%
    group_by(Region, cc_100) %>% 
    mutate(day_since_cc_100 = seq_len(n())-1) %>% ungroup %>%
    group_by(Region, cc_1000) %>% 
    mutate(day_since_cc_1000 = seq_len(n())-1) %>% ungroup %>%
    
    group_by(Region, cd_10) %>% 
    mutate(day_since_cd_10 = seq_len(n())-1) %>% ungroup %>%
    group_by(Region, cd_100) %>% 
    mutate(day_since_cd_100 = seq_len(n())-1) %>% ungroup %>%
    group_by(Region, cd_1000) %>% 
    mutate(day_since_cd_1000 = seq_len(n())-1) %>% ungroup 
  
  # doubling function
  doubling <- function(double = 2, start = 10, xmax = 100, step = 0.1) {
    
    x <- seq(0, xmax, step)
    y <- start * 2^(x/double) 
    return(data.frame(x= x, y = y, 
                      Doubling = paste0("Every ", double, " Days"),
                      Country = NA))
  }
  
  doubling_segment <- function(double = 2, start = 10, xmax = 100, step = 0.1) {
    
    x <- seq(0, xmax, step)
    y <- start * 2^(x/double) 
    return(data.frame(x= x[1], xend = x[length(x)],
                      y = y[1], yend = y[length(y)], 
                      Doubling = paste0("Every ", double, " Days"),
                      Country = NA))
  }
  
  
  doubling_lines_10 <- do.call(rbind, lapply(c(1:3,7), function(x){
    doubling(x, start = 10, xmax = 100)
  }))
  
  
  gg <- d %>% 
    group_by(Region) %>% filter(cd_10) %>% ungroup %>% 
    rename(Days = day_since_cd_10,
           Deaths = Cumulative_Deaths,
           Country = Region) %>% 
    highlight_key(~Country, "Select Countries") %>% 
    ggplot(aes(x=Days, y = Deaths, color = Continent, group = Country)) +
    geom_line(data = doubling_lines_10, aes(x=x, y=y, group = Doubling), 
              alpha = 0.5, color = "grey", linetype = "dashed", inherit.aes = FALSE) +
    geom_line() +
    scale_y_log10(limits = c(10, max(d$Cumulative_Deaths)*1.2), 
                  breaks = c(10,20,100, 200, 500, 1000, 2000, 5000, 10000, 20000),
                  labels = scales::comma) + 
    xlab("Days Since 10 Deaths") +
    ylab("Cumulative Cases Since 10 Deaths") +
    xlim(c(0, max(d$day_since_cd_10))) +
    theme_bw() + 
    scale_color_discrete(name="") +
    theme(axis.title.y = ggplot2::element_text(margin=ggplot2::margin(10, 200, 10, 10)),
          axis.title.x = ggplot2::element_text(margin=ggplot2::margin(10, 10, 10, 10)),
          axis.line = element_line(),
          panel.border = element_blank())
  
  p <- plotly::ggplotly(gg, tooltip = c("Country", "Region", "Cases", "Days"))
  
  p_10 <- p %>% 
    group_by(Country) %>% 
    plotly::highlight(on = c("plotly_click"), off = "plotly_doubleclick",debounce = 100, 
                      selected = plotly::attrs_selected(showlegend = FALSE, mode = "lines+markers" ), 
                      persistent = TRUE,
                      defaultValues = country,
                      selectize = TRUE)
  
  p_10
  
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



cases_plot <- function(out, data, date = Sys.Date(), date_0) {

  o1 <- squire:::calibrate_output_parsing(
    out, 
    date_0 = date_0
  )
  
  
  
  gg_cases <- squire:::plot_calibration_cases_barplot(o1, data = data, forecast = 0) + 
    ggplot2::xlim(c(date - 28, date))
  
  gg_cases + ggplot2::theme(legend.position = c(0,1), 
                            legend.justification = c(0,1), 
                            legend.direction = "horizontal") + 
    facet_zoom2(ylim = c(0, max(abs(diff(data$cases))*2)), zoom.size = 0.5) +
    ggtitle("Plot on right zoomed in on reported cases")
  
}



deaths_plot <- function(out, data,date_0, date = Sys.Date()) {

  o1 <- squire:::calibrate_output_parsing(
    out, 
    date_0 = date_0
  )
  
  gg_deaths <- squire:::plot_calibration_deaths_barplot(o1, data = data, forecast = 14, cumulative = FALSE) 
  gg_deaths 
    
  # geom_label(
  #   data = data.frame(x = c(as.Date(data$date[max(which(data$deaths == max(data$deaths)))]),date),
  #                     y = c(max(o1$y[o1$compartment == "deaths" & o1$date < (date+14)])*0.85,
  #                           max(o1$y[o1$compartment == "deaths" & o1$date < (date+14)])*0.75),
  #                     label=c("Calibration Date",as.character(date))), 
  #   aes(x=x, y=y, label=label), inherit.aes = FALSE) +
  #   ggplot2::theme(legend.position = c(0,1), 
  #                  legend.justification = c(0,1), 
  #                  legend.direction = "horizontal") + 
  #   
  #   
  
  
}


healthcare_plot <- function(out, data) {
  
  o1 <- squire:::calibrate_output_parsing(
    out, 
    date_0 = date_0
  )
  
  cowplot::plot_grid(squire:::plot_calibration_healthcare_barplot(df = o1, data = data, what = "hospital"),
                     squire:::plot_calibration_healthcare_barplot(df = o1, data = data, what = "ICU"),
                     ncol=2)
  
}
