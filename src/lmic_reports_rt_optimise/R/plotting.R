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

cumulative_deaths_plot <- function(country, excess = FALSE) {
  if(excess){
    d <- readRDS("excess_deaths.Rds")
    #convert to daily
    d <- d %>%
      mutate(
        deaths = deaths/as.numeric(date_end - date_start),
        date = date_start
      ) %>%
      group_by(iso3c) %>%
      complete(date = seq(min(date), max(date_end), by = 1)) %>%
      arrange(iso3c, date) %>%
      fill(deaths, .direction = "down") %>%
      ungroup() %>%
      transmute(
        iso3c  = iso3c ,
        Region = countrycode::countrycode(iso3c , "iso3c", "country.name"),
        date = date,
        deaths = deaths
      )
    d$Region[d$iso3c == "KSV"] <- "Kosovo"
    d <- d %>%
      left_join(
        readRDS("reported_covid.Rds") %>%
          ungroup() %>%
          select(iso3c, cases, date),
        by = c("iso3c", "date")
      ) %>%
      mutate(cases = if_else(is.na(cases), 0, cases))
  } else {
    d <- readRDS("reported_covid.Rds") %>%
      mutate(Region = countrycode::countrycode(iso3c , "iso3c", "country.name"))
  }
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
    arrange(date) %>%
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

  df_deaths_latest <- df_deaths[df_deaths$date == max(df_deaths$date),]
  continent <- na.omit(unique(df$Continent[df$Region == country]))

  if(excess){
    y_label <- "Cumulative Excess Deaths (Logarithmic Scale)"
  } else {
    y_label <- "Cumulative Deaths (Logarithmic Scale)"
  }

  gg_deaths <- ggplot(df_deaths[which(df_deaths$Continent == continent), ], aes(x=day_since, y=Cum_Deaths, group = Region)) +
    geom_line(data = doubling_lines_deaths, aes(x=x, y=y, linetype = Doubling),alpha=0.7, inherit.aes = FALSE,color = "black") +
    geom_line(show.legend = FALSE, color = "grey", alpha = 0.6) +
    geom_line(data = df_deaths[which(df_deaths$Region == country), ], color = "red", lwd = 1) +
    geom_point(data = df_deaths[which(df_deaths$Region == country), ], color = "red") +
    geom_point(data = df_deaths_latest[which(df_deaths_latest$Continent == continent), ], alpha = 0.5, show.legend = FALSE) +
    ggrepel::geom_text_repel(data =  df_deaths_latest[which(df_deaths_latest$Continent == continent), ],
                             aes(label = Region), show.legend = FALSE, min.segment.length = 2,nudge_x = 1) +
    scale_y_log10(limits=c(start, max(df_deaths$Cum_Deaths[df_deaths$Continent == continent],na.rm=TRUE))) +
    xlim(limits=c(0, max(df_deaths$day_since[df_deaths$Continent == continent],na.rm=TRUE))) +
    theme_bw() +
    ylab(y_label) +
    xlab(paste("Days Since", start, "Deaths"))

  gg_deaths

}


plotly_style <- function(country) {

  #ecdc <- readRDS("ecdc_all.rds")
  ecdc <- readRDS("reported_covid.Rds")

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
    arrange(date) %>%
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

  cases_exist <- nrow(data) > 0

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
                                                col = .data$compartment))
  if(cases_exist){
    gg_cases <- gg_cases + ggplot2::geom_bar(data = data,
                      mapping = ggplot2::aes(x = .data$date, y = .data$cases,
                                             fill = "Reported"),
                      stat = "identity",
                      show.legend = TRUE,
                      inherit.aes = FALSE)
  }
  gg_cases <- gg_cases + ggplot2::geom_ribbon(data = pd_group,
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

  if(cases_exist){
    gg_cases + ggplot2::theme(legend.position = "top",
                              legend.justification = c(0,1),
                              legend.direction = "horizontal") +
      facet_zoom2(ylim = c(0, max((data$cases)*1)), zoom.size = 0.5) +
      ggtitle("Plot on right zoomed in on reported cases") +
      geom_vline(xintercept = date, linetype = "dashed")
  } else {
    gg_cases
  }
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
    ggplot2::ylab("Daily Number of Infections") +
    ggplot2::theme_bw()  +
    ggplot2::scale_y_continuous(expand = c(0,0)) +
    ggplot2::scale_fill_manual(name = "", labels = (c("Estimated")),
                               values = (c("#3f8ea7"))) +
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
    start_date <- as.Date(out$inputs$start_date)
  } else {
    start_date <- as.Date(
      min(out$inputs$data$date_start[out$inputs$data$deaths > 0])
    )
  }

  date <- as.Date(date)
  #set deaths to daily
  data <- out$inputs$data %>% mutate(deaths = deaths/as.numeric(date_end - date_start))
  gg <- plot(out, summarise = TRUE, replicates = FALSE, particle_fit = TRUE)
  ymax <- max(data$deaths, gg$layers[[1]]$data$ymax[gg$layers[[1]]$data$date<=(as.Date(date)+forecast)])

  gg <- gg +
    geom_segment(data = data,
                 mapping = aes(x=date_start, xend = date_end, y=deaths, yend = deaths)) +
    ggplot2::geom_vline(xintercept = date, linetype = "dashed") +
    ggplot2::theme_bw()  +
    ggplot2::scale_y_continuous(limits = c(0, ymax+1)) +
    ggplot2::scale_x_date(date_breaks = "2 month", date_labels = "%b %Y",
                          limits = c(start_date, date + forecast),
                          expand = c(0, 0)) +
    # ggplot2::scale_fill_manual(name = "", labels = rev(c("Estimated")),
    #                            values = (c("#c59e96"))) +
    # ggplot2::scale_color_manual(name = "", labels = rev(c("Estimated")),
    #                             values = (c("#c59e96"))) +
    ggplot2::ylab("") +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1, colour = "black"),
                   axis.title.x = ggplot2::element_blank(),
                   panel.grid.major.x = ggplot2::element_blank(),
                   panel.grid.minor.x = ggplot2::element_blank(),
                   panel.border = ggplot2::element_blank(),
                   panel.background = ggplot2::element_blank(),
                   axis.line = ggplot2::element_line(colour = "black"))

  rg <- gg + ylab("Daily Deaths") +
    #ylim(c(0, max(gg$layers[[1]]$data[gg$layers[[1]]$data$x < date+1,]$ymax,out[[wh]]$inputs$data$deaths+1))) +
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

deaths_plot <-  function(proj, proj_surge, df_excess, df_cases, date_0, date,
                         forecast) {
  # build it from scratch
  if(is.null(proj_surge)){
    r_list <- list(proj)
  } else {
    r_list <- list(proj, proj_surge)
  }

  start_date <- date_0

  pd_list <- lapply(r_list, function(x) {

    df <- x[x$compartment == "deaths" & x$date >= start_date, ]
    df_summ <- group_by(df, date) %>%
      summarise(deaths = median(y, na.rm = TRUE),
                ymin = quantile(y, 0.025, na.rm = TRUE),
                ymax = quantile(y, 0.975, na.rm = TRUE))

    return(df_summ)

  })

  # append scenarios
  scenarios <- c("Current healthcare", "Surge in healthcare")
  for(i in seq_along(pd_list)) {
    pd_list[[i]]$Scenario <- scenarios[i]
    pd_list[[i]]$Scenario <- scenarios[i]
  }

  pds <- do.call(rbind, pd_list) %>% ungroup

  # Plot
  pds <- pds[pds$date <= date + forecast,]
  p <- ggplot(pds, aes(date, deaths, ymin = ymin, ymax = ymax, col = Scenario, fill = Scenario)) +
    geom_line() + geom_ribbon(alpha = 0.25, col = NA)

  # Add remaining formatting
  gg <- p +
    ggplot2::scale_color_discrete(name = "") +
    ggplot2::scale_fill_discrete(guide = "none") +
    ggplot2::theme_bw() +
    theme(legend.position = "top") +
    guides(linetype = "none")

  if(nrow(df_cases) == 0){
    limits <- c(min(df_excess$date_start[which(df_excess$deaths>0)]), date + forecast)
  } else {
    limits <- c(min(df_cases$date[which(df_cases$deaths>0)]), date + forecast)
  }

  if(is.null(df_excess)){
    gg2 <- gg +
      geom_point(data = df_cases, mapping = aes(x=date, y=deaths,shape="Reported Deaths"), inherit.aes = FALSE) +
      ggplot2::theme_bw()  +
      ggplot2::scale_x_date(date_breaks = "1 month", date_labels = "%b %Y",
                            limits = limits,
                            expand = c(0, 0))
  } else {
    gg2 <- gg +
      geom_segment(data = df_excess, mapping = aes(x = date_start, xend = date_end, y = deaths, yend = deaths), inherit.aes = FALSE, size = 1) +
      ggplot2::theme_bw()  +
      ggplot2::scale_x_date(date_breaks = "1 month", date_labels = "%b %Y",
                            limits = limits,
                            expand = c(0, 0))
  }

  gg2 <- gg2 +
    ggplot2::geom_vline(xintercept = date, linetype = "dashed")  +
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
  #don't bother with seperate plots
  gg2 + theme(legend.position = "none") + ylab("Daily Deaths") + ggtitle("Model Fit & 28 Day Projection")

}


healthcare_plot_contrast <- function(projs, what, date, forecast){
  #get correct dates
  df <- projs[projs$date <=  date + forecast + 1 &
                projs$date > date - 7, ] %>%
    #only the needed compartment
    filter(compartment == what) %>%
    group_by(Scenario, date) %>%
    #summarise over dates
    summarise(y = mean(y), .groups = "drop") %>%
    #format scenarios
    mutate(
      Scenario = stringr::str_remove(Scenario, "r_") %>%
        stringr::str_remove("_scenario") %>%
        stringr::str_remove("_leg") %>%
        stringr::str_to_title(),
      Scenario = ordered(Scenario, level = c("Pessimistic", "Central", "Optimistic", "Reverse", "Maintain", "Mitigation"))
    ) %>%
    arrange(date, Scenario)

  # y axis
  if (what == "ICU_demand") {
    title <- "ICU Demand"
  } else if(what == "hospital_demand") {
    title <- "Hospital Bed Demand"
  } else if(what == "infections"){
    title <- "Daily Infections"
  }

  ggplot2::ggplot(df, ggplot2::aes(x = .data$date,
                                   y = .data$y,
                                   fill = .data$Scenario)) +
    geom_bar(
      stat = "identity",
      position = "identity",
      show.legend = TRUE) +
    ggplot2::ylab(title) +
    ggplot2::theme_bw()  +
    ggplot2::scale_y_continuous(expand = c(0,0)) +
    ggplot2::scale_fill_manual(name = "", labels = sort(unique(df$Scenario)),
                               values = c("#9eeccd", "#c59e96", "#3f8ea7")) +
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


rt_creation_vaccine <- function(out, max_date) {

  iso3c <- squire::get_population(out$parameters$country)$iso3c[1]

  #SOME KIND OF CHECK FOR NOT FITTING (might not be neccessary)
  #ret <- rt_creation_simple_out_vaccine(out, date_0, max_date)
  # # simple here change
  # if ("projection_args" %in% names(out)) {
  #   ret[which(ret$date>date_0),-(1:2)] <- ret[which(ret$date>date_0),-(1:2)]*out$projection_args$R0_change
  # }

  # impact of immunity ratios
  ratios <- get_immunity_ratios(out, max_date = max_date, vaccine = TRUE)

  # create the Rt data frame
  rt <- purrr::map_dfr(seq_along(out$samples), function(y) {

    dates <- out$inputs$start_date + seq(0, tail(out$samples[[y]]$tt_R0, 1))
    dates <- c(dates, seq.Date(tail(dates, 1) + 1,
                               as.Date(max_date),
                               1))

    Rt <- as.numeric(squire.page:::block_interpolate(seq_along(dates) - 1, out$samples[[y]]$R0, out$samples[[y]]$tt_R0))

    df <- data.frame(
      "Rt" = Rt,
      "Reff" = Rt*ratios[[y]],
      "date" = dates,
      "iso" = iso3c,
      rep = y,
      stringsAsFactors = FALSE)
    df$pos <- seq_len(nrow(df))
    return(df)
  })
  rt$date <- as.Date(rt$date)

  rt <- rt[,c(3,4,1,2,5,6)]

  sum_rt <- dplyr::group_by(rt, date) %>%
    arrange(date) %>%
    dplyr::summarise(compartment = "Rt",
                     y_025 = quantile(Rt, 0.025),
                     y_25 = quantile(Rt, 0.25),
                     y_median = median(Rt),
                       y_mean = mean(Rt),
                       y_75 = quantile(Rt, 0.75),
                       y_975 = quantile(Rt, 0.975))

  sum_reff <- dplyr::group_by(rt, date) %>%
    arrange(date) %>%
    dplyr::summarise(compartment = "Reff",
                     y_025 = quantile(Reff, 0.025),
                     y_25 = quantile(Reff, 0.25),
                     y_median = median(Reff),
                     y_mean = mean(Reff),
                     y_75 = quantile(Reff, 0.75),
                     y_975 = quantile(Reff, 0.975))

  ret <- rbind(sum_rt, sum_reff)

  return(ret)
}


rt_creation_simple_out_vaccine <- function(out, date_0, max_date) {

  iso3c <- squire::get_population(out$parameters$country)$iso3c[1]

  get_immunity_ratios_simple_vaccine <- function(out, max_date) {

    mixing_matrix <- squire:::process_contact_matrix_scaled_age(
      out$parameters$contact_matrix_set[[1]],
      out$parameters$population
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
      max_date <- max(rownames(out$output))
    }
    t_now <- which(as.Date(rownames(out$output)) == max_date)


    # prop susceptible
    prop_susc <- lapply(seq_len(dim(out$output)[3]), function(x) {

      susceptible <- array(
        out$output[seq_len(t_now),index$S,x],
        dim=c(t_now, dim(index$S))
      )

      # We divide by the total population
      prop_susc <- sweep(susceptible, 2, pop, FUN='/')

      # We multiply by the effect of vaccines on onward infectiousness
      prop_susc <- vapply(
        seq_len(nrow(prop_susc)),
        FUN = function(i){prop_susc[i,,]*out$odin_parameters$vaccine_efficacy_infection[1,,]},
        FUN.VALUE = prop_susc[1,,]
      )

      prop_susc <- aperm(prop_susc, c(3,1,2))

      return(prop_susc)
    } )

    relative_R0_by_age <- prob_hosp*dur_ICase + (1-prob_hosp)*dur_IMild

    adjusted_eigens <- lapply(prop_susc, function(x) {

      unlist(lapply(seq_len(nrow(x)), function(y) {
        if(any(is.na(x[y,,]))) {
          return(NA)
        } else {
          Re(eigen(mixing_matrix*rowSums(x[y,,]*relative_R0_by_age))$values[1])
        }
      }))

    })

    mat <- squire:::process_contact_matrix_scaled_age(out$parameters$contact_matrix_set[[1]],
                                                      out$parameters$population)

    betas <- lapply(rep(out$parameters$R0, dim(out$output)[3]), function(x) {
      squire:::beta_est_explicit(dur_IMild = dur_IMild,
                                 dur_ICase = dur_ICase,
                                 prob_hosp = prob_hosp,
                                 mixing_matrix = mat,
                                 R0 = x)
    })


    ratios <- lapply(seq_along(betas), function(x) {
      (betas[[x]] * adjusted_eigens[[x]]) / out$parameters$R0
    })

    return(ratios)
  }

  # impact of immunity ratios
  ratios <- get_immunity_ratios_simple_vaccine(out, max_date)

  # create the Rt data frame
  rts <- lapply(seq_len(dim(out$output)[3]), function(y) {

    Rt <- rep(out$parameters$R0, length(ratios[[y]]))

    df <- data.frame(
      "Rt" = Rt,
      "Reff" = Rt*tail(na.omit(ratios[[y]]),length(Rt)),
      "date" = head(rownames(out$output),length(Rt)),
      "iso" = iso3c,
      rep = y,
      stringsAsFactors = FALSE)
    df$pos <- seq_len(nrow(df))
    return(df)
  } )

  rt <- do.call(rbind, rts)
  rt$date <- as.Date(rt$date)

  rt <- rt[,c(3,4,1,2,5,6)]

  new_rt_all <- rt %>%
    group_by(iso, rep) %>%
    arrange(date) %>%
    complete(date = seq.Date(min(rt$date), date_0, by = "days"))

  column_names <- colnames(new_rt_all)[-c(1,2,3)]
  new_rt_all <- fill(new_rt_all, all_of(column_names), .direction = c("down"))
  new_rt_all <- fill(new_rt_all, all_of(column_names), .direction = c("up"))

  sum_rt <- dplyr::group_by(new_rt_all, date) %>%
    dplyr::summarise(compartment = "Rt",
                     y_025 = quantile(Rt, 0.025),
                     y_25 = quantile(Rt, 0.25),
                     y_median = median(Rt),
                     y_mean = mean(Rt),
                     y_75 = quantile(Rt, 0.75),
                     y_975 = quantile(Rt, 0.975))

  sum_reff <- dplyr::group_by(new_rt_all, date) %>%
    dplyr::summarise(compartment = "Reff",
                     y_025 = quantile(Reff, 0.025),
                     y_25 = quantile(Reff, 0.25),
                     y_median = median(Reff),
                     y_mean = mean(Reff),
                     y_75 = quantile(Reff, 0.75),
                     y_975 = quantile(Reff, 0.975))

  return(rbind(sum_rt, sum_reff))

}

nim_sq_simulation_plot_prep <- function(x,
                                        var_select,
                                        q = c(0.025, 0.975),
                                        summary_f = mean,
                                        x_var = "t",
                                        ...) {

  pd <- nimue_format(x, var_select = var_select, ...)

  pd <- pd %>%
    dplyr::mutate(x = .data[[x_var]])

  # t sometimes seems to be being rounded weirdly
  if(x_var == "t") {
    pd$x <- round(pd$x, ceiling(log10(1/x$parameters$dt)))
  }

  # remove any NA rows (due to different start dates)
  if(sum(is.na(pd$t) | is.na(pd$y))>0) {
    pd <- pd[-which(is.na(pd$t) | is.na(pd$y)),]
  }

  # Format summary data
  pds <- pd %>%
    dplyr::group_by(.data$x, .data$compartment) %>%
    dplyr::summarise(ymin = stats::quantile(.data$y, q[1]),
                     ymax = stats::quantile(.data$y, q[2]),
                     y = summary_f(.data$y))

  return(list(pd = pd, pds = pds))

}

# variant_timings_plot <- function(variant_timings, date_range){
#   df <- variant_timings %>%
#     mutate(new_start_date = end_date, end_date = lead(start_date, 1, default = date_range[2]),
#            start_date = new_start_date) %>%
#     select(!new_start_date) %>%
#     rbind(
#       variant_timings %>%
#         mutate(variant = "Shift-Period")
#     ) %>%
#     add_row(end_date = min(variant_timings$start_date), start_date = date_range[1],
#             variant = "Wild") %>%
#     arrange(start_date) %>%
#     transmute(
#       t_end = as.numeric(end_date - min(start_date)),
#       t_start = as.numeric(start_date - min(start_date)),
#       centre = (t_end + t_start)/2,
#       date_label = start_date,
#       variant = variant,
#       angle = if_else(variant == "Shift-Period", 90, 0)
#     )
#   height <- 1
#   off_set <- height*0.025
#   date_pos <- c(-off_set, height + off_set)
#   v_pos <- c(1, 0)
#   text_size <- 5
#   ggplot(
#     df,
#     aes(xmin = t_start, xmax = t_end, ymin = 0, ymax = height, fill = variant)
#   ) +
#     geom_rect(show.legend = FALSE) +
#     geom_text(aes(x = centre, y = height/2, label = variant, angle = angle), size = text_size) +
#     geom_text(aes(x = t_start, y = rep(date_pos, nrow(df))[seq_len(nrow(df))], label = date_label,
#                   vjust = rep(v_pos, nrow(df))[seq_len(nrow(df))]),
#               size = text_size*0.75, inherit.aes = FALSE) +
#     geom_segment(aes(
#       x = t_start, xend = t_start, y = rep(date_pos, nrow(df))[seq_len(nrow(df))],
#       yend = rep(c(0, height), nrow(df))[seq_len(nrow(df))]
#     ))+
#     coord_cartesian(xlim = c(-30, as.numeric(diff(date_range))), ylim = date_pos) +
#     theme_void()
# }


variant_timings_plot <- function(variant_timings, date_range){
  df <- variant_timings %>%
    mutate(new_start_date = end_date, end_date = lead(start_date, 1, default = date_range[2]),
           start_date = new_start_date) %>%
    select(!new_start_date) %>%
    rbind(
      variant_timings %>%
        mutate(variant = "Shift-Period")
    ) %>%
    add_row(end_date = min(variant_timings$start_date), start_date = date_range[1],
            variant = "Wild") %>%
    arrange(start_date) %>%
    mutate(
      centre = as.numeric(end_date - start_date)/2 + start_date,
      variant = variant,
      angle = if_else(variant == "Shift-Period", 90, 0),
      end_date = if_else(end_date == date_range[2], date_range[2] + 100, end_date)
    )
  height <- 1
  text_size <- 5
  ggplot(
    df,
    aes(xmin = start_date, xmax = end_date, ymin = 0, ymax = height, fill = variant)
  ) +
    geom_rect(show.legend = FALSE, alpha = 0.75) +
    geom_text(aes(x = centre, y = height/2, label = variant, angle = angle), size = text_size, alpha = 0.75) +
    coord_cartesian(xlim = date_range, ylim = c(0, height)) +
    theme(axis.line.y=element_blank(),
          axis.text.y=element_blank(),
          axis.ticks.y=element_blank(),
          axis.title.y=element_blank(),
          panel.grid.minor.y=element_blank(),
          panel.grid.major.y=element_blank(),
          panel.background = element_blank()) +
    ggplot2::scale_x_date(date_labels = "%b %Y", date_breaks = "1 months") +
    ggplot2::labs(title = "Dominant Variants:") +
    xlab("Date")
}

variant_timings_index_plot <- function(variant_timings, date_range){
  height <- 1
  text_size <- 5
  offset <- 1.05
  df <- variant_timings %>%
    mutate(new_start_date = end_date, end_date = lead(start_date, 1, default = date_range[2]),
           start_date = new_start_date) %>%
    select(!new_start_date) %>%
    rbind(
      variant_timings %>%
        mutate(variant = "")
    ) %>%
    add_row(end_date = min(variant_timings$start_date), start_date = date_range[1],
            variant = "Wild-Type") %>%
    arrange(start_date) %>%
    mutate(
      centre = as.numeric(end_date - start_date)/2 + start_date,
      variant = if_else(variant == "Omicron Sub-Variant", "Sub-Variant", variant),
      y = if_else(as.numeric(end_date - start_date)/as.numeric(max(end_date) - min(start_date)) < 0.1, height*offset, height/2),
      end_date = if_else(end_date == date_range[2], date_range[2] + 100, end_date)
    )
  ggplot(
    df,
    aes(xmin = start_date, xmax = end_date, ymin = 0, ymax = height, fill = variant)
  ) +
    geom_rect(show.legend = FALSE, alpha = 0.75) +
    geom_text(aes(x = centre, y = y, label = variant), size = text_size,
              alpha = 0.75) +
    coord_cartesian(xlim = date_range, ylim = c(0, height*offset)) +
    theme(axis.line.y=element_blank(),
          axis.text.y=element_blank(),
          axis.ticks.y=element_blank(),
          axis.title.y=element_blank(),
          panel.grid.minor.y=element_blank(),
          panel.grid.major.y=element_blank(),
          panel.background = element_blank()) +
    ggplot2::scale_x_date(date_labels = "%b %Y", date_breaks = "3 months") +
    ggplot2::labs(title = "") +
    xlab("Date")
}
