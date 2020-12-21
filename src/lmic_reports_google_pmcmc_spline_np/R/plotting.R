## colors
first <- "#3f8ea7"
second <- "#c59e96"
third <- "#9eeccd"


match_clean <- function(a,b, quiet=TRUE){
  a <- gsub("[[:punct:][:space:]]","",tolower(stringi::stri_trans_general(a, "latin-ascii")))
  b <- gsub("[[:punct:][:space:]]","",tolower(stringi::stri_trans_general(b, "latin-ascii")))
  ret <- match(a,b)
  if(sum(is.na(ret)>0)){
    dists <- stringdist::seq_distmatrix(lapply(a,utf8ToInt),lapply(b,utf8ToInt))
    ret[is.na(ret)] <- apply(dists[which(is.na(ret)),,drop=FALSE],1,which.min)
    if(!quiet){
      return(unique(cbind(a,b[ret])))
    }
  }
  return(ret)
}

continent_match <- function(region) {
  
  suppressWarnings(continent <- countrycode::countrycode(region, origin = 'country.name', destination = 'continent'))
  
  if(region=="Eswatini") continent <- "Africa"
  if(region=="United State of America") continent <- "Americas"
  if(region=="Isle_of_Man") continent <- "Europe"             
  if(region=="Kosovo") continent <- "Europe"                  
  if(region=="Netherlands_Antilles") continent <- "Americas"    
  if(region=="Saint_Lucia") continent <- "Americas"             
  if(region=="South_Korea") continent <- "Asia"             
  if(region=="United_States_of_America") continent <- "Americas"
  
  return(continent)
  
}

cumulative_deaths_plot <- function(country) {
  
  d <- readRDS("ecdc_all.rds")
  d <- readRDS("jhu_all.rds")
  # d <- readRDS("worldometers_all.rds")
  d$Region[d$Region=="Congo"] <- "Republic of Congo"
  d$Region[d$Region=="United_Republic_of_Tanzania"] <- "Tanzania"
  d$Region[d$Region=="CuraÃ§ao"] <- "Curacao"
  start <- 10
  
  country <- gsub("[[:punct:]]", "", country)
  country <- gsub(" ", "_", country)
  country <- d$Region[match_clean(country, d$Region)]
  
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
  
  country <- gsub("_" ," ", country)
  df$Region <- gsub("_" ," ", df$Region)
  
  df_deaths <- df %>% 
    filter(Cum_Deaths > start) %>% 
    mutate(day_since = seq_len(n())-1)
  
  doubling_lines_deaths <- do.call(rbind, lapply(c(2, 3, 5, 7), function(x){
    doubling(x, start = start, xmax = max(df_deaths$day_since))
  }))
  
  df_deaths_latest <- df_deaths[df_deaths$dateRep == max(df_deaths$dateRep),]
  continent <- na.omit(unique(df$Continent[df$Region == country]))
  
  gg_deaths <- ggplot(df_deaths[which(df_deaths$Continent == continent), ], aes(x=day_since, y=Cum_Deaths, group = Region)) + 
    geom_line(data = doubling_lines_deaths, aes(x=x, y=y, linetype = Doubling),alpha=0.7, inherit.aes = FALSE,color = "black") +
    geom_line(show.legend = FALSE, color = "grey", alpha = 0.6) +
    geom_line(data = df_deaths[which(df_deaths$Region == country), ], color = "red", lwd = 1) +
    geom_point(data = df_deaths[which(df_deaths$Region == country), ], color = "red") +
    geom_point(data = df_deaths_latest[which(df_deaths_latest$Continent == continent), ], alpha = 0.5, show.legend = FALSE) + 
    ggrepel::geom_text_repel(data =  df_deaths_latest[which(df_deaths_latest$Continent == continent), ],
                             aes(label = Region), show.legend = FALSE, min.segment.length = 2,nudge_x = 1) + 
    scale_y_log10(limits=c(start, max(df_deaths$Cum_Deaths[df_deaths$Continent == continent]))) +
    xlim(limits=c(0, max(df_deaths$day_since[df_deaths$Continent == continent]))) +
    theme_bw() +
    ylab("Cumulative Deaths (Logarithmic Scale)") +
    xlab(paste("Days Since", start, "Deaths"))
  
  gg_deaths
  
}


plotly_style <- function(country) {
  
  #ecdc <- readRDS("ecdc_all.rds")
  ecdc <- readRDS("jhu_all.rds")
  
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



cases_plot <- function(df, data, date = Sys.Date(), date_0) {
  
  # day
  df$day <- as.Date(as.character(df$date))
  
  # split to correct dates
  sub <- df[df$compartment %in% c("infections") &
              df$date <=  date + 1,] %>%
    dplyr::group_by(.data$day, .data$replicate) %>%
    dplyr::summarise(y = sum(.data$y)) %>%
    dplyr::filter(.data$day <= date) %>% 
    dplyr::filter(!is.na(.data$y))
  
  
  pd_group <- dplyr::group_by(sub, .data$day) %>%
    dplyr::summarise(quants = list(quantile(.data$y, c(0.025, 0.25, 0.5, 0.75, 0.975))),
                     ymin = .data$quants[[1]][1],
                     ymax = .data$quants[[1]][5],
                     yinner_min = .data$quants[[1]][2],
                     yinner_max = .data$quants[[1]][4],
                     y = mean(.data$y))
  
  # format cases
  data$cases <- rev(c(tail(data$cases,1), diff(rev(data$cases))))
  
  # Plot
  gg_cases <- ggplot2::ggplot(sub, ggplot2::aes(x = .data$day,
                                                y = .data$y,
                                                col = .data$compartment)) +
    ggplot2::geom_ribbon(data = pd_group,
                         mapping = ggplot2::aes(ymin = .data$ymin,
                                                ymax = .data$ymax,
                                                fill = "Estimated"),
                         color = "white",
                         alpha = 0.2,
                         size = 0,
                         show.legend = TRUE) +
    ggplot2::geom_ribbon(data = pd_group,
                         mapping = ggplot2::aes(ymin = .data$yinner_min,
                                                ymax = .data$yinner_max,
                                                fill = "Estimated"),
                         color = "white",
                         alpha = 0.8,
                         size = 0,
                         show.legend = TRUE) +
    ggplot2::geom_bar(data = data,
                      mapping = ggplot2::aes(x = .data$date, y = .data$cases,
                                             fill = "Reported"),
                      stat = "identity",
                      show.legend = TRUE,
                      inherit.aes = FALSE) +
    ggplot2::ylab("Daily Number of Infections") +
    ggplot2::theme_bw()  +
    ggplot2::scale_y_continuous(expand = c(0,0)) +
    ggplot2::scale_fill_manual(name = "", labels = (c("Estimated", "Reported")),
                               values = (c("#3f8ea7","#c59e96"))) +
    ggplot2::scale_x_date(date_breaks = "1 month", date_labels = "%b %d") +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1, colour = "black"),
                   axis.title.x = ggplot2::element_blank(),
                   panel.grid.major.x = ggplot2::element_blank(),
                   panel.grid.minor.x = ggplot2::element_blank(),
                   panel.border = ggplot2::element_blank(),
                   panel.background = ggplot2::element_blank(),
                   axis.line = ggplot2::element_line(colour = "black")
    )
  
  gg_cases + ggplot2::theme(legend.position = "top", 
                            legend.justification = c(0,1),
                            legend.direction = "horizontal") + 
    facet_zoom2(ylim = c(0, max((data$cases)*1)), zoom.size = 0.5) +
    ggtitle("Plot on right zoomed in on reported cases") +
    geom_vline(xintercept = date, linetype = "dashed")
  
}

cases_plot_single <- function(df, data, date = Sys.Date(), date_0) {
  
  # day
  df$day <- as.Date(as.character(df$date))
  
  # split to correct dates
  sub <- df[df$compartment %in% c("infections") &
              df$date <=  date + 1,] %>%
    dplyr::group_by(.data$day, .data$replicate) %>%
    dplyr::summarise(y = sum(.data$y)) %>%
    dplyr::filter(.data$day <= date) %>% 
    dplyr::filter(!is.na(.data$y))
  
  
  pd_group <- dplyr::group_by(sub, .data$day) %>%
    dplyr::summarise(quants = list(quantile(.data$y, c(0.025, 0.25, 0.5, 0.75, 0.975))),
                     ymin = .data$quants[[1]][1],
                     ymax = .data$quants[[1]][5],
                     yinner_min = .data$quants[[1]][2],
                     yinner_max = .data$quants[[1]][4],
                     y = mean(.data$y))
  
  # format cases
  #data$cases <- rev(c(tail(data$cases,1), diff(rev(data$cases))))
  
  # Plot
  gg_cases <- ggplot2::ggplot(sub, ggplot2::aes(x = .data$day,
                                                y = .data$y,
                                                col = .data$compartment)) +
    ggplot2::geom_ribbon(data = pd_group,
                         mapping = ggplot2::aes(ymin = .data$ymin,
                                                ymax = .data$ymax,
                                                fill = "Estimated"),
                         color = "white",
                         alpha = 0.2,
                         size = 0,
                         show.legend = TRUE) +
    ggplot2::geom_ribbon(data = pd_group,
                         mapping = ggplot2::aes(ymin = .data$yinner_min,
                                                ymax = .data$yinner_max,
                                                fill = "Estimated"),
                         color = "white",
                         alpha = 0.8,
                         size = 0,
                         show.legend = TRUE) +
    ggplot2::geom_bar(data = data,
                      mapping = ggplot2::aes(x = .data$date, y = .data$cases,
                                             fill = "Reported"),
                      stat = "identity",
                      show.legend = TRUE,
                      inherit.aes = FALSE) +
    ggplot2::ylab("Daily Number of Infections") +
    ggplot2::theme_bw()  +
    ggplot2::scale_y_continuous(expand = c(0,0)) +
    ggplot2::scale_fill_manual(name = "", labels = (c("Estimated", "Reported")),
                               values = (c("#3f8ea7","#c59e96"))) +
    ggplot2::scale_x_date(date_breaks = "1 month", date_labels = "%b %d") +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1, colour = "black"),
                   axis.title.x = ggplot2::element_blank(),
                   panel.grid.major.x = ggplot2::element_blank(),
                   panel.grid.minor.x = ggplot2::element_blank(),
                   panel.border = ggplot2::element_blank(),
                   panel.background = ggplot2::element_blank(),
                   axis.line = ggplot2::element_line(colour = "black")
    )
  
  gg_cases + ggplot2::theme(legend.position = "top", 
                            legend.justification = c(0,1),
                            legend.direction = "horizontal") +
    geom_vline(xintercept = date, linetype = "dashed")
  
}

deaths_plot_single <- function(out, data, date_0, date = Sys.Date(), 
                               forecast = 14, single = FALSE, full = TRUE) {
  
  
  if(full) {
    start_date <- min(as.Date(out$replicate_parameters$start_date))
  } else {
    start_date <- min(data$date[which(data$deaths>0)])
  }
    
  
  if("pmcmc_results" %in% names(out)) {
    wh <- "pmcmc_results"
  } else {
    wh <- "scan_results"
  }
  
  date <- as.Date(date)
  gg <- plot(out, "deaths", date_0 = date_0, x_var = "date", summary_f = median) 
  ymax <- max(out[[wh]]$inputs$data$deaths, gg$layers[[1]]$data$ymax[gg$layers[[1]]$data$x<=(as.Date(date)+forecast)])
  
  gg <- gg + 
    geom_point(data = out[[wh]]$inputs$data, mapping = aes(x=date, y=deaths,shape="Reported")) +
    ggplot2::geom_vline(xintercept = date, linetype = "dashed") +
    ggplot2::theme_bw()  +
    ggplot2::scale_y_continuous(limits = c(0, ymax+1)) +
    ggplot2::scale_x_date(date_breaks = "2 weeks", date_labels = "%b %d",
                          limits = c(start_date, date + forecast),
                          expand = c(0, 0)) +
    ggplot2::scale_fill_manual(name = "", labels = rev(c("Estimated")),
                               values = (c("#c59e96"))) +
    ggplot2::scale_color_manual(name = "", labels = rev(c("Estimated")),
                                values = (c("#c59e96"))) +
    ggplot2::ylab("") +
    ggplot2::scale_shape_manual(name = "", values = 20) +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1, colour = "black"),
                   axis.title.x = ggplot2::element_blank(),
                   panel.grid.major.x = ggplot2::element_blank(),
                   panel.grid.minor.x = ggplot2::element_blank(),
                   panel.border = ggplot2::element_blank(),
                   panel.background = ggplot2::element_blank(),
                   axis.line = ggplot2::element_line(colour = "black")) 
  
  rg <- gg + ylab("Daily Deaths") +
    ylim(c(0, max(gg$layers[[1]]$data[gg$layers[[1]]$data$x < date+1,]$ymax,out[[wh]]$inputs$data$deaths+1))) +
    xlim(c(start_date, date)) +
    theme(legend.position = "none") +
    ggtitle("Model Fit up to Current Day")
  
  if(single) {
    return(rg)
  } else {
    
    gg <- gg + theme(legend.position = "none") + ylab("") + ggtitle("Model Fit & 28 Day Projection")
    
    leg <- cowplot::get_legend(gg)
    return(cowplot::plot_grid(leg, 
                              cowplot::plot_grid(rg, gg+theme(legend.position = "none"),rel_widths=c(0.66,1)),
                              ncol = 1, rel_heights = c(0.1,1)))
  }
  
  
}

deaths_plot_single_surge <- function(out, out2, data, date_0, date = Sys.Date(), 
                                     forecast = 14) {
  
  
  if("pmcmc_results" %in% names(out)) {
    wh <- "pmcmc_results"
  } else {
    wh <- "scan_results"
  }
  
  # build it from scratch
  r_list <- list(out, out2)
  
  pd_list <- lapply(r_list, FUN = squire:::squire_simulation_plot_prep,
                    var_select = "deaths",
                    x_var = "date", 
                    q = c(0.025, 0.975),
                    summary_f = median,
                    date_0 = date_0)
  
  # append scenarios
  scenarios <- c("Current healthcare", "Surge in healthcare")
  for(i in seq_along(scenarios)) {
    pd_list[[i]]$pd$Scenario <- scenarios[i]
    pd_list[[i]]$pds$Scenario <- scenarios[i]
  }
  
  pds <- do.call(rbind, lapply(pd_list, "[[", "pds")) %>% ungroup
  pd <- do.call(rbind, lapply(pd_list, "[[", "pd"))
  
  # Plot
  pds <- pds %>% filter(x <= date_0 + forecast)
  p <- ggplot2::ggplot(data = pds, 
                       ggplot2::aes(x = .data$x, y = .data$y, col = Scenario))
  
  p <- p + ggplot2::geom_line(data = pds,
                              ggplot2::aes(x = .data$x, y = .data$y,
                                           col = .data$Scenario,
                                           linetype = .data$compartment))
  p <- p + ggplot2::geom_ribbon(data = pds,
                                ggplot2::aes(x = .data$x,
                                             ymin = .data$ymin,
                                             ymax = .data$ymax,
                                             fill = .data$Scenario,
                                             linetype = .data$compartment),
                                alpha = 0.25, col = NA)
  
  # Add remaining formatting
  gg <- p +
    ggplot2::scale_color_discrete(name = "") +
    ggplot2::scale_fill_discrete(guide = FALSE) +
    ggplot2::xlab("Time") +
    ggplot2::ylab("N") +
    ggplot2::theme_bw() + 
    theme(legend.position = "top") + 
    guides(linetype = FALSE) 
  
  date <- as.Date(date)
  ymax <- max(out[[wh]]$inputs$data$deaths, 
              gg$layers[[1]]$data$ymax[gg$layers[[1]]$data$x<=(as.Date(date)+forecast)])
  
  gg2 <- gg + 
    geom_point(data = out[[wh]]$inputs$data, mapping = aes(x=date, y=deaths,shape="Reported Deaths"), inherit.aes = FALSE) +
    ggplot2::geom_vline(xintercept = date, linetype = "dashed") +
    ggplot2::theme_bw()  +
    ggplot2::scale_x_date(date_breaks = "2 week", date_labels = "%b %d",
                          limits = c(min(data$date[which(data$deaths>0)]), date + forecast),
                          expand = c(0, 0)) +
    ggplot2::scale_fill_manual(name = "", labels = (c("Estimated with Current Healthcare Capacity", 
                                                      "Estimated with Surge in Healthcare Capacity")),
                               values = (c("#c59e96","#3f8ea7"))) +
    ggplot2::scale_color_manual(name = "", labels = (c("Estimated with Current Healthcare Capacity", 
                                                       "Estimated with Surge in Healthcare Capacity")),
                                values = (c("#c59e96","#3f8ea7"))) +
    ggplot2::ylab("Daily Deaths") +
    ggplot2::scale_shape_manual(name = "", values = 20) +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1, colour = "black"),
                   axis.title.x = ggplot2::element_blank(),
                   panel.grid.major.x = ggplot2::element_blank(),
                   panel.grid.minor.x = ggplot2::element_blank(),
                   panel.border = ggplot2::element_blank(),
                   panel.background = ggplot2::element_blank(),
                   axis.line = ggplot2::element_line(colour = "black")) 
  
  gg2 <- gg2 + ggplot2::theme(legend.position = "top", 
                              legend.justification = c(0.5,1),
                              legend.direction = "horizontal") +
    geom_vline(xintercept = date, linetype = "dashed")
  
  rg <- gg2 + ylab("Daily Deaths") +
    ylim(c(0, max(gg2$data[gg2$data$x < date+1,]$ymax, out[[wh]]$inputs$data$deaths+1))) +
    xlim(c(min(data$date[which(data$deaths>0)]), date)) +
    theme(legend.position = "none") +
    ggtitle("Model Fit up to Current Day")
  
  gg2 <- gg2 + theme(legend.position = "none") + ylab("") + ggtitle("Model Fit & 28 Day Projection")
  
  leg <- cowplot::get_legend(gg2)
  cowplot::plot_grid(leg, 
                     cowplot::plot_grid(rg, gg2+theme(legend.position = "none"),rel_widths=c(0.66,1)),
                     ncol = 1, rel_heights = c(0.1,1))
  
  
}


deaths_plot_contrast_triple <- function(o1, o2, o3, data, date_0, date = Sys.Date(), 
                                        forecast = 14, cumulative = TRUE) {
  
  o1$Scenario <- "No"
  o2$Scenario <- "Yes"
  o3$Scenario <- "Worse"
  df <- rbind(o1,o2,o3)
  
  # day
  df$day <- as.Date(as.character(df$date))
  
  # split to correct dates
  if(!cumulative) {
    sub <- df[df$compartment == "deaths" &
                df$date <=  date + forecast + 1,]  %>%
      dplyr::group_by(.data$day, .data$replicate, .data$Scenario) %>%
      dplyr::summarise(y = mean(.data$y), n=dplyr::n()) %>%
      dplyr::filter(.data$day <= date + forecast) %>% 
      dplyr::filter(!is.na(.data$y))
    
    title <- "Daily Deaths"
    
    # format deaths
    data$deaths <- rev(c(tail(data$deaths,1), diff(rev(data$deaths))))
    
  } else {
    sub <- df[df$compartment == "D" &
                df$date <=  date + forecast + 1,]  %>%
      dplyr::group_by(.data$day, .data$replicate, .data$Scenario) %>%
      dplyr::summarise(y = mean(.data$y), n=dplyr::n()) %>%
      dplyr::filter(.data$day <= date + forecast) %>% 
      dplyr::filter(!is.na(.data$y))
    
    title <- "Cumulative Deaths"
  }
  
  pd_group <- dplyr::group_by(sub, .data$day, .data$Scenario) %>%
    dplyr::summarise(quants = list(quantile(.data$y, c(0.025, 0.25, 0.5, 0.75, 0.975))),
                     ymin = round(.data$quants[[1]][1]),
                     ymax = round(.data$quants[[1]][5]),
                     yinner_min = round(.data$quants[[1]][2]),
                     yinner_max = round(.data$quants[[1]][4]),
                     y = mean(.data$y),
                     n = dplyr::n())
  
  
  # Plot
  gg_healthcare <- ggplot2::ggplot(sub,
                                   ggplot2::aes(x = .data$day,
                                                y = .data$y,
                                                fill = .data$compartment)) +
    ggplot2::geom_ribbon(data = pd_group,
                         mapping = ggplot2::aes(ymin = .data$ymin,
                                                ymax = .data$ymax,
                                                fill = Scenario),
                         color = "white",
                         alpha = 0.2,
                         size = 0,
                         show.legend = TRUE) +
    ggplot2::geom_ribbon(data = pd_group,
                         mapping = ggplot2::aes(ymin = .data$yinner_min,
                                                ymax = .data$yinner_max,
                                                fill = Scenario),
                         color = "white",
                         alpha = 0.8,
                         size = 0,
                         show.legend = TRUE) +
    ggplot2::geom_vline(xintercept = date, linetype = "dashed") +
    ggplot2::theme_bw()  +
    ggplot2::ylab(title) +
    ggplot2::scale_y_continuous(expand = c(0,0)) +
    ggplot2::scale_x_date(date_breaks = "2 weeks", date_labels = "%b %d",
                          limits = c(date - 7, date + forecast)) +
    ggplot2::scale_fill_manual(name = "", labels = (c( "Maintain Status Quo","Relax Interventions 50%","Additional 50% Reduction")),
                               values = (c("#9eeccd","#c59e96","#3f8ea7"))) +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1, colour = "black"),
                   axis.title.x = ggplot2::element_blank(),
                   panel.grid.major.x = ggplot2::element_blank(),
                   panel.grid.minor.x = ggplot2::element_blank(),
                   panel.border = ggplot2::element_blank(),
                   panel.background = ggplot2::element_blank(),
                   axis.line = ggplot2::element_line(colour = "black"),
                   legend.position = "top", 
                   legend.justification = c(0,1), 
                   legend.direction = "horizontal") 
  
  gg_healthcare
  
  
}

cases_contrast_triple <- function(o1, o2, o3, data, date_0, date = Sys.Date(), 
                                  forecast = 14) {
  
  o1$Scenario <- "No"
  o2$Scenario <- "Yes"
  o3$Scenario <- "Worse"
  df <- rbind(o1,o2,o3)
  
  # day
  df$day <- as.Date(as.character(df$date))
  
  # split to correct dates
  sub <- df[df$compartment == "infections" &
              df$date <=  date + forecast + 1,]  %>%
    dplyr::group_by(.data$day, .data$replicate, .data$Scenario) %>%
    dplyr::summarise(y = mean(.data$y), n=dplyr::n()) %>%
    dplyr::filter(.data$day <= date + forecast) %>% 
    dplyr::filter(!is.na(.data$y))
  
  title <- "Daily Cases"
  
  
  pd_group <- dplyr::group_by(sub, .data$day, .data$Scenario) %>%
    dplyr::summarise(quants = list(quantile(.data$y, c(0.025, 0.25, 0.5, 0.75, 0.975))),
                     ymin = round(.data$quants[[1]][1]),
                     ymax = round(.data$quants[[1]][5]),
                     yinner_min = round(.data$quants[[1]][2]),
                     yinner_max = round(.data$quants[[1]][4]),
                     y = mean(.data$y),
                     n = dplyr::n())
  ymax <- max(pd_group$ymax[pd_group$day<(date+forecast) & pd_group$day>(date-forecast)])
  
  # Plot
  gg_healthcare <- ggplot2::ggplot(sub,
                                   ggplot2::aes(x = .data$day,
                                                y = .data$y,
                                                color = .data$Scenario)) +
    # ggplot2::geom_line(sub, 
    #                    mapping = ggplot2::aes(group = interaction(Scenario, replicate),
    #                                           y = .data$y, x=.data$day),
    #                    show.legend = TRUE, se = FALSE,alpha=0.1) +
    ggplot2::geom_smooth(data = pd_group,
                         mapping = ggplot2::aes(y = .data$y, 
                                                ymin = .data$ymin,
                                                ymax = .data$ymax,
                                                color = Scenario),
                         show.legend = FALSE, se = FALSE) +
    ggplot2::geom_smooth(data = pd_group,
                         mapping = ggplot2::aes(y = .data$ymin, 
                                                color = Scenario),
                         linetype = "dashed",
                         show.legend = FALSE, se = FALSE, alpha = 0.2) +
    ggplot2::geom_smooth(data = pd_group,
                         mapping = ggplot2::aes(y = .data$ymax, 
                                                color = Scenario),
                         linetype = "dashed",
                         show.legend = FALSE, se = FALSE, alpha = 0.2) +
    ggplot2::geom_vline(xintercept = date, linetype = "dashed") +
    ggplot2::theme_bw()  +
    ggplot2::ylab(title) +
    ggplot2::scale_y_continuous(expand = c(0,0),limits=c(0,ymax)) +
    ggplot2::scale_x_date(date_breaks = "2 weeks", date_labels = "%b %d",
                          limits = c(date - 7, date + forecast)) +
    ggplot2::scale_color_manual(name = "", labels = (c( "Maintain Status Quo","Relax Interventions 50%","Additional 50% Reduction")),
                                values = (c("#9eeccd","#c59e96","#3f8ea7"))) +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1, colour = "black"),
                   axis.title.x = ggplot2::element_blank(),
                   panel.grid.major.x = ggplot2::element_blank(),
                   panel.grid.minor.x = ggplot2::element_blank(),
                   panel.border = ggplot2::element_blank(),
                   panel.background = ggplot2::element_blank(),
                   axis.line = ggplot2::element_line(colour = "black"),
                   legend.position = "top", 
                   legend.justification = c(0.5,1), 
                   legend.direction = "horizontal") + 
    guides(colour = guide_legend(override.aes = list(alpha = 1)))
  
  gg_healthcare
  
  
}

cases_contrast_triple_bars <- function(o1, o2, o3, data, date_0, date = Sys.Date(), forecast = 14, what = "ICU_demand") {
  
  o1$Scenario <- "No"
  o2$Scenario <- "Yes"
  o3$Scenario <- "Worse"
  df <- rbind(o1,o2,o3)
  
  # day
  df$day <- as.Date(as.character(df$date))
  
  # split to correct dates
  sub <- df[df$compartment %in% "infections" &
              df$date <=  date + forecast + 1,] %>%
    dplyr::group_by(.data$day, .data$replicate, .data$Scenario) %>%
    dplyr::summarise(y = mean(.data$y), n=dplyr::n()) %>%
    dplyr::filter(.data$day <= date + forecast) %>% 
    dplyr::filter(.data$day >= date - 7) %>% 
    dplyr::filter(!is.na(.data$y))
  
  sub$Scenario <- factor(sub$Scenario, levels = c( "No","Worse", "Yes"))
  pd_group <- dplyr::group_by(sub, .data$day, .data$Scenario) %>%
    dplyr::summarise(quants = list(quantile(.data$y, c(0.25, 0.5, 0.75))),
                     ymin = .data$quants[[1]][1],
                     y = mean(.data$y),
                     ymax = .data$quants[[1]][3])
  ymax <- max(pd_group$y[pd_group$day<=(date+forecast) & pd_group$day>(date-forecast)])
  
  # Plot
  suppressMessages(suppressWarnings(
    gg_healthcare <- ggplot2::ggplot(sub, ggplot2::aes(x = .data$day,
                                                       y = .data$y, 
                                                       fill = .data$Scenario)) +
      ggplot2::geom_bar(data = pd_group[pd_group$Scenario=="Worse",],
                        mapping = ggplot2::aes(x = .data$day, y = .data$y, fill = .data$Scenario),
                        stat = "identity",
                        position = "identity",
                        show.legend = TRUE,
                        inherit.aes = FALSE) +
      ggplot2::geom_bar(data = pd_group[pd_group$Scenario=="No",],
                        mapping = ggplot2::aes(x = .data$day, y = .data$y, fill = .data$Scenario),
                        stat = "identity",
                        position = "identity",
                        show.legend = TRUE,
                        inherit.aes = FALSE) +
      ggplot2::geom_bar(data = pd_group[pd_group$Scenario=="Yes",],
                        mapping = ggplot2::aes(x = .data$day, y = .data$y, fill = .data$Scenario),
                        stat = "identity",
                        position = "identity",
                        show.legend = TRUE,
                        inherit.aes = FALSE) +
      ggplot2::theme_bw()  +
      ggplot2::ylab("Daily Infections") +
      ggplot2::scale_y_continuous(expand = c(0,0), limits = c(0, ymax)) +
      ggplot2::scale_fill_manual(name = "", labels = (c( "Maintain Status Quo","Relax Interventions 50%","Additional 50% Reduction")),
                                 values = (c("#9eeccd","#c59e96","#3f8ea7"))) +
      ggplot2::scale_x_date(date_breaks = "2 weeks", date_labels = "%b %d", 
                            expand = c(0, 0)) +
      ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1, colour = "black"),
                     axis.title.x = ggplot2::element_blank(),
                     panel.grid.major.x = ggplot2::element_blank(),
                     panel.grid.minor.x = ggplot2::element_blank(),
                     legend.position = "top", 
                     legend.justification = c(0.5,1),
                     legend.direction = "horizontal",
                     panel.border = ggplot2::element_blank(),
                     panel.background = ggplot2::element_blank(),
                     axis.line = ggplot2::element_line(colour = "black")) +
      geom_vline(xintercept = date, linetype = "dashed")
  ))
  
  gg_healthcare
  
}

healthcare_plot_contrast <- function(o1, o2, data, date_0, date = Sys.Date(), forecast = 14, what = "ICU_demand") {
  
  o1$Scenario <- "No"
  o2$Scenario <- "Yes"
  df <- rbind(o1,o2)
  
  # day
  df$day <- as.Date(as.character(df$date))
  
  # split to correct dates
  sub <- df[df$compartment %in% what &
              df$date <=  date + forecast + 1,] %>%
    dplyr::group_by(.data$day, .data$replicate, .data$Scenario) %>%
    dplyr::summarise(y = mean(.data$y), n=dplyr::n()) %>%
    dplyr::filter(.data$day <= date_0 + forecast) %>% 
    dplyr::filter(!is.na(.data$y))
  
  sub$Scenario <- factor(sub$Scenario, levels = c("No","Yes"))
  pd_group <- dplyr::group_by(sub, .data$day, .data$Scenario) %>%
    dplyr::summarise(quants = list(quantile(.data$y, c(0.25, 0.5, 0.75))),
                     ymin = .data$quants[[1]][1],
                     y = mean(.data$y),
                     ymax = .data$quants[[1]][3])
  
  
  # y axis
  if (what == "ICU_demand") {
    title <- "ICU Demand"
  } else if(what == "hospital_demand") {
    title <- "Hospital Bed Demand"
  }
  
  # Plot
  suppressMessages(suppressWarnings(
    gg_healthcare <- ggplot2::ggplot(sub, ggplot2::aes(x = .data$day,
                                                       y = .data$y, 
                                                       fill = .data$Scenario)) +
      ggplot2::geom_bar(data = pd_group[pd_group$Scenario=="No",],
                        mapping = ggplot2::aes(x = .data$day, y = .data$y, fill = .data$Scenario),
                        stat = "identity",
                        position = "identity",
                        show.legend = TRUE,
                        inherit.aes = FALSE) +
      ggplot2::geom_bar(data = pd_group[pd_group$Scenario=="Yes",],
                        mapping = ggplot2::aes(x = .data$day, y = .data$y, fill = .data$Scenario),
                        stat = "identity",
                        position = "identity",
                        show.legend = TRUE,
                        inherit.aes = FALSE) +
      ggplot2::ylab(title) +
      ggplot2::theme_bw()  +
      ggplot2::scale_y_continuous(expand = c(0,0)) +
      ggplot2::scale_fill_manual(name = "", labels = (c("Maintain Status Quo", "Additional 50% Reduction")),
                                 values = rev(c("#3f8ea7","#c59e96"))) +
      ggplot2::scale_x_date(date_breaks = "2 weeks", date_labels = "%b %d", limits = c(date_0-7, date_0 + forecast)) +
      ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1, colour = "black"),
                     axis.title.x = ggplot2::element_blank(),
                     panel.grid.major.x = ggplot2::element_blank(),
                     panel.grid.minor.x = ggplot2::element_blank(),
                     legend.position = "top", 
                     legend.justification = c(0,1),
                     legend.direction = "horizontal",
                     panel.border = ggplot2::element_blank(),
                     panel.background = ggplot2::element_blank(),
                     axis.line = ggplot2::element_line(colour = "black")) +
      geom_vline(xintercept = date, linetype = "dashed")
  ))
  gg_healthcare
  
}

healthcare_plot_contrast_lines <- function(o1, o2, data, date_0, date = Sys.Date(), forecast = 14, what = "ICU_demand") {
  
  o1$Scenario <- "No change to epidemic"
  o2$Scenario <- "Mitigation (50% reduction from today)"
  df <- rbind(o1,o2)
  
  # day
  df$day <- as.Date(as.character(df$date))
  
  # split to correct dates
  sub <- df[df$compartment %in% what &
              df$date <=  date + forecast + 1,] %>%
    dplyr::group_by(.data$day, .data$replicate, .data$Scenario) %>%
    dplyr::summarise(y = mean(.data$y), n=dplyr::n()) %>%
    dplyr::filter(.data$day <= date_0 + forecast) %>% 
    dplyr::filter(!is.na(.data$y))
  
  pd_group <- dplyr::group_by(sub[sub$day>date_0-7,], .data$day, .data$Scenario) %>%
    dplyr::summarise(quants = list(quantile(.data$y, c(0.25, 0.5, 0.75))),
                     ymin = .data$quants[[1]][1],
                     ymin_t = t.test(.data$y)$conf.int[1],
                     y = mean(.data$y),
                     ymax = .data$quants[[1]][3],
                     ymin_t = t.test(.data$y)$conf.int[2])
  
  # y axis
  if (what == "ICU_demand") {
    title <- "ICU Demand"
  } else if(what == "hospital_demand") {
    title <- "Hospital Bed Demand"
  }
  
  # Plot
  suppressMessages(suppressWarnings(
    gg_healthcare <- ggplot2::ggplot(sub, ggplot2::aes(x = .data$day,
                                                       y = .data$y, 
                                                       fill = .data$Scenario)) +
      ggplot2::geom_line(data = pd_group,
                         mapping = ggplot2::aes(x = .data$day, y = .data$y, fill = .data$Scenario),
                         show.legend = TRUE,
                         inherit.aes = FALSE) +
      ggplot2::geom_errorbar(data = pd_group,
                             mapping = ggplot2::aes(x = .data$day, ymin = .data$ymin,ymax=.data$ymax, fill = .data$Scenario),
                             stat = "identity",
                             show.legend = TRUE,
                             inherit.aes = FALSE) +
      ggplot2::geom_vline(xintercept = date, linetype = "dashed") +
      ggplot2::ylab(title) +
      ggplot2::theme_bw()  +
      ggplot2::scale_y_continuous(expand = c(0,0)) +
      ggplot2::scale_fill_manual(name = "", labels = rev(c("Maintain Status Quo", "Additional 50% Reduction")),
                                 values = rev(c("#3f8ea7","#c59e96"))) +
      ggplot2::scale_x_date(date_breaks = "2 weeks", date_labels = "%b %d", limits = c(date_0-7, date_0 + forecast)) +
      ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1, colour = "black"),
                     axis.title.x = ggplot2::element_blank(),
                     panel.grid.major.x = ggplot2::element_blank(),
                     panel.grid.minor.x = ggplot2::element_blank(),
                     legend.position = "top", 
                     legend.justification = c(0,1),
                     legend.direction = "horizontal",
                     panel.border = ggplot2::element_blank(),
                     panel.background = ggplot2::element_blank(),
                     axis.line = ggplot2::element_line(colour = "black"))
  ))
  gg_healthcare
  
}

healthcare_plot_contrast_triple <- function(o1, o2, o3, data, date_0, date = Sys.Date(), forecast = 14, what = "ICU_demand") {
  
  o1$Scenario <- "No"
  o2$Scenario <- "Yes"
  o3$Scenario <- "Worse"
  df <- rbind(o1,o2,o3)
  
  # day
  df$day <- as.Date(as.character(df$date))
  
  # split to correct dates
  sub <- df[df$compartment %in% what &
              df$date <=  date + forecast + 1,] %>%
    dplyr::group_by(.data$day, .data$replicate, .data$Scenario) %>%
    dplyr::summarise(y = mean(.data$y), n=dplyr::n()) %>%
    dplyr::filter(.data$day <= date + forecast) %>%
    dplyr::filter(.data$day >= date - 7) %>% 
    dplyr::filter(!is.na(.data$y))
  
  sub$Scenario <- factor(sub$Scenario, levels = c( "No","Worse", "Yes"))
  pd_group <- dplyr::group_by(sub, .data$day, .data$Scenario) %>%
    dplyr::summarise(quants = list(quantile(.data$y, c(0.25, 0.5, 0.75))),
                     ymin = .data$quants[[1]][1],
                     y = mean(.data$y),
                     ymax = .data$quants[[1]][3])
  
  
  # y axis
  if (what == "ICU_demand") {
    title <- "ICU Demand"
  } else if(what == "hospital_demand") {
    title <- "Hospital Bed Demand"
  }
  
  # Plot
  suppressMessages(suppressWarnings(
    gg_healthcare <- ggplot2::ggplot(sub, ggplot2::aes(x = .data$day,
                                                       y = .data$y, 
                                                       fill = .data$Scenario)) +
      ggplot2::geom_bar(data = pd_group[pd_group$Scenario=="Worse",],
                        mapping = ggplot2::aes(x = .data$day, y = .data$y, fill = .data$Scenario),
                        stat = "identity",
                        position = "identity",
                        show.legend = TRUE,
                        inherit.aes = FALSE) +
      ggplot2::geom_bar(data = pd_group[pd_group$Scenario=="No",],
                        mapping = ggplot2::aes(x = .data$day, y = .data$y, fill = .data$Scenario),
                        stat = "identity",
                        position = "identity",
                        show.legend = TRUE,
                        inherit.aes = FALSE) +
      ggplot2::geom_bar(data = pd_group[pd_group$Scenario=="Yes",],
                        mapping = ggplot2::aes(x = .data$day, y = .data$y, fill = .data$Scenario),
                        stat = "identity",
                        position = "identity",
                        show.legend = TRUE,
                        inherit.aes = FALSE) +
      ggplot2::ylab(title) +
      ggplot2::theme_bw()  +
      ggplot2::scale_y_continuous(expand = c(0,0)) +
      ggplot2::scale_fill_manual(name = "", labels = (c( "Maintain Status Quo","Relax Interventions 50%","Additional 50% Reduction")),
                                 values = (c("#9eeccd","#c59e96","#3f8ea7"))) +
      ggplot2::scale_x_date(date_breaks = "2 weeks", date_labels = "%b %d",
                            expand = c(0, 0)) +
      ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1, colour = "black"),
                     axis.title.x = ggplot2::element_blank(),
                     panel.grid.major.x = ggplot2::element_blank(),
                     panel.grid.minor.x = ggplot2::element_blank(),
                     legend.position = "top", 
                     legend.justification = c(0,1),
                     legend.direction = "horizontal",
                     panel.border = ggplot2::element_blank(),
                     panel.background = ggplot2::element_blank(),
                     axis.line = ggplot2::element_line(colour = "black")) +
      geom_vline(xintercept = date, linetype = "dashed")
  ))
  
  gg_healthcare
  
}

plot_scan <- function(x, what = "likelihood", log = TRUE) {
  
  if (what == "likelihood") {
    
    # create df
    df <- data.frame("z" = as.numeric(x$mat_log_ll))
    df$x <- x$x
    df$y <- sort(rep(x$y, length(x$x)))
    
    # make plot
    gg <- ggplot(data=df,aes(x=x,y=y,fill=-z)) +
      geom_tile(color = "white") + xlab("R0") + ylab("Date") +
      scale_x_continuous( expand = c(0, 0)) +
      scale_y_date( expand = c(0, 0)) +
      scale_fill_gradient( name = "-Log L.") +
      theme_bw() +
      theme(panel.grid = element_blank(),
            panel.border = element_blank(),
            plot.title = element_text(hjust = 0.5)) + 
      ggtitle("Likelihood")
    
    if(log) {
      gg <- gg + scale_fill_gradient( name = "-Log L.", trans = 'log' )
    }
    
    
  } else if (what == "probability") {
    
    # create df
    df <- data.frame("z" = as.numeric(x$renorm_mat_LL))
    df$x <- x$x
    df$y <- sort(rep(x$y, length(x$x)))
    
    # make plot
    gg <- ggplot(data=df,aes(x=x,y=y,fill=z)) +
      geom_tile(color = "white") + xlab("R0") + ylab("Date") +
      scale_x_continuous( expand = c(0, 0)) +
      scale_fill_gradient( name = "Probability") + 
      scale_y_date( expand = c(0, 0)) +
      theme_bw() +
      theme(panel.grid = element_blank(),
            panel.border = element_blank(),
            plot.title = element_text(hjust = 0.5)) + 
      ggtitle("Probability")
    
    if(log) {
      gg <- gg + scale_fill_gradient( name = "Probability", trans = 'log' )
    }
  }
  return(gg)
}

intervention_plot <- function(res, date) {
  
  if (nrow(res) == 0) {
    res <- data.frame(date = as.Date((as.Date(date)-56):(as.Date(date)),lubridate::origin),
                      "S1" = FALSE,
                      "S2" = FALSE,
                      "S3" = FALSE,
                      "S6" = FALSE)
  } 
  
  res %>% dplyr::rename("School Closure" = S1,
                        "Work Closure" = S2,
                        "Public Events\nBanned" = S3,
                        "Lockdown" = S6) %>%  
    tidyr::pivot_longer(cols = c("School Closure", "Work Closure", 
                                 "Public Events\nBanned", "Lockdown")) %>% 
    ggplot(aes(x=date, y=as.factor(name),fill=as.factor(value>0))) + 
    geom_tile(color="black") + 
    ylab("") + 
    xlab("Date") + 
    scale_fill_manual(name="Intervention Active", values = rev(c("#3f8ea7","#c59e96"))) +
    ggplot2::scale_x_date(date_breaks = "2 weeks", date_labels = "%b %d",
                          expand = c(0,0)) +
    scale_y_discrete(expand = c(0,0)) +
    theme_bw() + 
    theme(legend.position = "top")
  
  
}

intervention_plot_google <- function(res, date, data, forecast, start_date) {
  
  if(nrow(res) == 0) {
    res <- data.frame(date = data$date, C = 1, observed = FALSE)
  }
  
  date <- as.Date(date)
  
  ggplot(res[res$date <= date+forecast & res$date >= start_date,], 
         aes(x = date, y = C, color = observed)) + 
    geom_point() + 
    scale_color_discrete(name = "Observed") +
    theme_bw() + 
    theme(legend.position = "top") + 
    ylab("% Mobility") + 
    xlab("Date")
  
}

rt_plot <- function(out) {
  
  if("pmcmc_results" %in% names(out)) {
    wh <- "pmcmc_results"
  } else {
    wh <- "scan_results"
  }
  
  date <- as.Date(date)
  date_0 <- date
  
  # create the Rt data frame
  rts <- lapply(seq_len(length(out$replicate_parameters$R0)), function(y) {
    
    tt <- squire:::intervention_dates_for_odin(dates = out$interventions$date_R0_change, 
                                               change = out$interventions$R0_change, 
                                               start_date = out$replicate_parameters$start_date[y],
                                               steps_per_day = 1/out$parameters$dt)
    
    if(wh == "scan_results") {
      Rt <- c(out$replicate_parameters$R0[y], 
              vapply(tt$change, out[[wh]]$inputs$Rt_func, numeric(1), 
                     R0 = out$replicate_parameters$R0[y], Meff = out$replicate_parameters$Meff[y])) 
    } else {
      Rt <- squire:::evaluate_Rt_pmcmc(
        R0_change = tt$change, 
        date_R0_change = tt$dates, 
        R0 = out$replicate_parameters$R0[y], 
        pars = as.list(out$replicate_parameters[y,]),
        Rt_args = out$pmcmc_results$inputs$Rt_args) 
    }
    
    df <- data.frame(
      "Rt" = Rt,
      "date" = tt$dates,
      "iso" = iso3c,
      rep = y,
      stringsAsFactors = FALSE)
    df$pos <- seq_len(nrow(df))
    return(df)
  } )
  
  rt <- do.call(rbind, rts)
  rt$date <- as.Date(rt$date)
  
  rt <- rt[,c(3,2,1,4,5)]
  
  new_rt_all <- rt %>%
    group_by(iso, rep) %>% 
    arrange(date) %>% 
    complete(date = seq.Date(min(rt$date), date_0, by = "days")) 
  
  column_names <- colnames(new_rt_all)[-c(1,2,3)]
  new_rt_all <- fill(new_rt_all, all_of(column_names), .direction = c("down"))
  new_rt_all <- fill(new_rt_all, all_of(column_names), .direction = c("up"))
  
  suppressMessages(sum_rt <- group_by(new_rt_all, iso, date) %>% 
                     summarise(Rt_min = quantile(Rt, 0.025),
                               Rt_q25 = quantile(Rt, 0.25),
                               Rt_q75 = quantile(Rt, 0.75),
                               Rt_max = quantile(Rt, 0.975),
                               Rt_median = median(Rt),
                               Rt = mean(Rt)))
  
  min_date <- min(as.Date(out$replicate_parameters$start_date))
  
  country_plot <- function(vjust = -1.2) {
    ggplot(sum_rt %>% filter(
      date > min_date & date <= as.Date(as.character(date_0+as.numeric(lubridate::wday(date_0))))), 
      aes(x=date, y = Rt, ymin=Rt_min, ymax = Rt_max, group = iso, fill = iso)) +
      geom_ribbon(fill = "#96c4aa") +
      geom_line(color = "#48996b") +
      geom_ribbon(mapping = aes(ymin = Rt_q25, ymax = Rt_q75), fill = "#48996b") +
      geom_hline(yintercept = 1, linetype = "dashed") +
      theme_bw() +
      theme(axis.text = element_text(size=12)) +
      xlab("") +
      scale_x_date(breaks = "2 weeks",
                   limits = as.Date(c(as.character(min_date),
                                      as.character(date_0+as.numeric(lubridate::wday(date_0))))), 
                   date_labels = "%d %b",
                   expand = c(0,0)) + 
      theme(legend.position = "none") +
      theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1, colour = "black"),
            panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
            panel.border = element_blank(),
            panel.background = element_blank(), axis.line = element_line(colour = "black")
      )
  }
  
  res <- list("plot" = suppressWarnings(country_plot()), "rts" = sum_rt)
  return(res)  
}


get_immunity_ratios <- function(out, max_date = NULL) {
  
  mixing_matrix <- squire:::process_contact_matrix_scaled_age(
    out$pmcmc_results$inputs$model_params$contact_matrix_set[[1]],
    out$pmcmc_results$inputs$model_params$population
  )
  
  dur_ICase <- out$parameters$dur_ICase
  dur_IMild <- out$parameters$dur_IMild
  prob_hosp <- out$parameters$prob_hosp
  
  # assertions
  squire:::assert_single_pos(dur_ICase, zero_allowed = FALSE)
  squire:::assert_single_pos(dur_IMild, zero_allowed = FALSE)
  squire:::assert_numeric(prob_hosp)
  squire:::assert_numeric(mixing_matrix)
  squire:::assert_square_matrix(mixing_matrix)
  squire:::assert_same_length(mixing_matrix[,1], prob_hosp)
  
  if(sum(is.na(prob_hosp)) > 0) {
    stop("prob_hosp must not contain NAs")
  }
  
  if(sum(is.na(mixing_matrix)) > 0) {
    stop("mixing_matrix must not contain NAs")
  }
  
  index <- squire:::odin_index(out$model)
  pop <- out$parameters$population
  
  if(is.null(max_date)) {
    max_date <- max(out$pmcmc_results$inputs$data$date)
  }
  t_now <- which(as.Date(rownames(out$output)) == max_date)
  prop_susc <- lapply(seq_len(dim(out$output)[3]), function(x) {
    t(t(out$output[seq_len(t_now), index$S, x])/pop)
  } )
  
  relative_R0_by_age <- prob_hosp*dur_ICase + (1-prob_hosp)*dur_IMild
  
  adjusted_eigens <- lapply(prop_susc, function(x) {
    
    unlist(lapply(seq_len(nrow(x)), function(y) {
      if(any(is.na(x[y,]))) {
        return(NA)
      } else {
        Re(eigen(mixing_matrix*x[y,]*relative_R0_by_age)$values[1])
      }
    }))
    
  })
  
  betas <- lapply(out$replicate_parameters$R0, function(x) {
    squire:::beta_est(squire_model = out$pmcmc_results$inputs$squire_model, 
                      model_params = out$pmcmc_results$inputs$model_params, 
                      R0 = x)
  })
  
  ratios <- lapply(seq_along(betas), function(x) {
    (betas[[x]] * adjusted_eigens[[x]]) / out$replicate_parameters$R0[[x]]
  })
  
  return(ratios)
}


rt_plot_immunity <- function(out) {
  
  iso3c <- squire::get_population(out$parameters$country)$iso3c[1]
  
  if("pmcmc_results" %in% names(out)) {
    wh <- "pmcmc_results"
  } else {
    wh <- "scan_results"
  }
  
  date <- max(as.Date(out$pmcmc_results$inputs$data$date))
  date_0 <- date
  
  # impact of immunity ratios
  ratios <- get_immunity_ratios(out)
  
  # create the Rt data frame
  rts <- lapply(seq_len(length(out$replicate_parameters$R0)), function(y) {
    
    tt <- squire:::intervention_dates_for_odin(dates = out$interventions$date_R0_change, 
                                               change = out$interventions$R0_change, 
                                               start_date = out$replicate_parameters$start_date[y],
                                               steps_per_day = 1/out$parameters$dt)
    
    if(wh == "scan_results") {
      Rt <- c(out$replicate_parameters$R0[y], 
              vapply(tt$change, out[[wh]]$inputs$Rt_func, numeric(1), 
                     R0 = out$replicate_parameters$R0[y], Meff = out$replicate_parameters$Meff[y])) 
    } else {
      Rt <- squire:::evaluate_Rt_pmcmc(
        R0_change = tt$change, 
        date_R0_change = tt$dates, 
        R0 = out$replicate_parameters$R0[y], 
        pars = as.list(out$replicate_parameters[y,]),
        Rt_args = out$pmcmc_results$inputs$Rt_args) 
    }
    
    df <- data.frame(
      "Rt" = Rt,
      "Reff" = Rt*tail(na.omit(ratios[[y]]),length(Rt)),
      "R0" = na.omit(Rt)[1]*tail(na.omit(ratios[[y]]),length(Rt)),
      "date" = tt$dates,
      "iso" = iso3c,
      rep = y,
      stringsAsFactors = FALSE)
    df$pos <- seq_len(nrow(df))
    return(df)
  } )
  
  rt <- do.call(rbind, rts)
  rt$date <- as.Date(rt$date)
  
  rt <- rt[,c(5,4,1,2,3,6,7)]
  
  new_rt_all <- rt %>%
    group_by(iso, rep) %>% 
    arrange(date) %>% 
    complete(date = seq.Date(min(rt$date), date_0, by = "days")) 
  
  column_names <- colnames(new_rt_all)[-c(1,2,3)]
  new_rt_all <- fill(new_rt_all, all_of(column_names), .direction = c("down"))
  new_rt_all <- fill(new_rt_all, all_of(column_names), .direction = c("up"))
  
  suppressMessages(sum_rt <- group_by(new_rt_all, iso, date) %>% 
                     summarise(Rt_min = quantile(Rt, 0.025),
                               Rt_q25 = quantile(Rt, 0.25),
                               Rt_q75 = quantile(Rt, 0.75),
                               Rt_max = quantile(Rt, 0.975),
                               Rt_median = median(Rt),
                               Rt = mean(Rt),
                               R0_min = quantile(R0, 0.025),
                               R0_q25 = quantile(R0, 0.25),
                               R0_q75 = quantile(R0, 0.75),
                               R0_max = quantile(R0, 0.975),
                               R0_median = median(R0),
                               R0 = mean(R0),
                               Reff_min = quantile(Reff, 0.025),
                               Reff_q25 = quantile(Reff, 0.25),
                               Reff_q75 = quantile(Reff, 0.75),
                               Reff_max = quantile(Reff, 0.975),
                               Reff_median = median(Reff),
                               Reff = mean(Reff)))
  
  min_date <- min(as.Date(out$replicate_parameters$start_date))
  
  country_plot <- function(vjust = -1.2) {
    ggplot(sum_rt %>% filter(
      date > min_date & date <= as.Date(as.character(date_0+as.numeric(lubridate::wday(date_0)))))) +
      geom_ribbon(mapping = aes(x=date, ymin=R0_min, ymax = R0_max, group = iso), fill = "#8cbbca") +
      geom_ribbon(mapping = aes(x = date, ymin = R0_q25, ymax = R0_q75, group = iso), fill = "#3f8da7") +
      geom_ribbon(mapping = aes(x=date, ymin=Reff_min, ymax = Reff_max, group = iso), fill = "#96c4aa") +
      geom_ribbon(mapping = aes(x = date, ymin = Reff_q25, ymax = Reff_q75, group = iso), fill = "#48996b") +
      geom_line(mapping = aes(x = date, y = Reff_median), color = "#48996b") +
      geom_hline(yintercept = 1, linetype = "dashed") +
      geom_hline(yintercept = sum_rt$R0_median[1], linetype = "dashed") +
      theme_bw() +
      theme(axis.text = element_text(size=12)) +
      xlab("") +
      ylab("Reff") +
      scale_x_date(breaks = "2 weeks",
                   limits = as.Date(c(as.character(min_date),
                                      as.character(date_0+as.numeric(lubridate::wday(date_0))))), 
                   date_labels = "%d %b",
                   expand = c(0,0)) + 
      theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1, colour = "black"),
            panel.grid.major = element_blank(), 
            panel.grid.minor = element_blank(),
            panel.border = element_blank(),
            panel.background = element_blank(), 
            axis.line = element_line(colour = "black")
      )
  }
  
  
  res <- list("plot" = suppressWarnings(country_plot()), "rts" = sum_rt)
  return(res)  
}



simple_pmcmc_plot <- function(out) {
  
  master <- squire:::create_master_chain(out$pmcmc_results, 0)
  master$chain <-  unlist(lapply(strsplit(rownames(master),".", fixed=TRUE), "[[", 1))
  master$iteration <-  as.numeric(unlist(lapply(strsplit(rownames(master),".", fixed=TRUE), "[[", 2)))
  
  par_pos <- seq_len(which(names(master) == "log_prior")-1)
  pars <- names(master)[par_pos]
  
  hists <- lapply(par_pos, function(i) {
    
    quants <- round(quantile(master[[pars[[i]]]][order(master$log_posterior, decreasing = TRUE)][1:1000], c(0.025, 0.5, 0.975), na.rm = TRUE),2)[c(2,1,3)]
    title <- paste0(pars[i],":\n ", quants[1], " (", quants[2], ", ", quants[3], ")")
    
    ggplot(master, mapping = aes_string(x = pars[i], color = "chain")) + 
      geom_freqpoly(stat = "density") + theme_bw() + theme(legend.position = "none") +
      theme(panel.border = element_blank(), axis.line = element_line()) + 
      ggtitle(title) + scale_color_brewer(type = "qual") +
      theme(axis.title = element_blank(),
            title = element_text(size = 8))
    
  })
  
  chains <- lapply(par_pos, function(i) {
    
    ggplot(master[seq(1,nrow(master),100),], mapping = aes_string(y = pars[i], x = "iteration", color = "chain")) + 
      geom_line() + theme_bw() +
      theme(panel.border = element_blank(), axis.line = element_line()) +
      theme(legend.position = "none") + scale_color_brewer(type = "qual")
    
  }) 
  
  p1 <- cowplot::plot_grid(plotlist = hists)
  p2 <- cowplot::plot_grid(plotlist = chains)
  line <- ggplot() + cowplot::draw_line(x = 1, y=0:10) + 
    theme(panel.background = element_blank(),
          axis.title = element_blank(), 
          axis.text = element_blank(), 
          axis.ticks = element_blank())
  cowplot::plot_grid(p1, line, p2, rel_widths = c(1,0.1,1), ncol = 3)
  
  
}
