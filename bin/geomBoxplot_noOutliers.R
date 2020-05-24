#modified version of geom_boxplot

require(ggplot2)
geom_boxplot_noOutliers <- function (mapping = NULL, data = NULL, stat = "boxplot",
                                     position = "dodge", outlier.colour = NULL,
                                     outlier.shape = NULL, outlier.size = NULL,
                                     notch = FALSE, notchwidth = .5, varwidth = FALSE,
                                     ...) {
  
  #outlier_defaults <- ggplot2:::Geom$find('point')$default_aes()
  
  #outlier.colour   <- outlier.colour %||% outlier_defaults$colour
  #outlier.shape    <- outlier.shape  %||% outlier_defaults$shape
  #outlier.size     <- outlier.size   %||% outlier_defaults$size
  
  GeomBoxplot_noOutliers$new(mapping = mapping, data = data, stat = stat,
                             position = position, outlier.colour = outlier.colour,
                             outlier.shape = outlier.shape, outlier.size = outlier.size, notch = notch,
                             notchwidth = notchwidth, varwidth = varwidth, ...)
}

GeomBoxplot_noOutliers <- proto(ggplot2:::Geom, {
  objname <- "boxplot_noOutliers"
  
  reparameterise <- function(., df, params) {
    df$width <- df$width %||%
      params$width %||% (resolution(df$x, FALSE) * 0.9)
    
    # if (!is.null(df$outliers)) {
    #    suppressWarnings({
    #      out_min <- vapply(df$outliers, min, numeric(1))
    #      out_max <- vapply(df$outliers, max, numeric(1))
    #    })
    #    
    #    df$ymin_final <- pmin(out_min, df$ymin)
    #    df$ymax_final <- pmax(out_max, df$ymax)
    #   }
    
    # if `varwidth` not requested or not available, don't use it
    if (is.null(params) || is.null(params$varwidth) || !params$varwidth || is.null(df$relvarwidth)) {
      df$xmin <- df$x - df$width / 2
      df$xmax <- df$x + df$width / 2
    } else {
      # make `relvarwidth` relative to the size of the largest group
      df$relvarwidth <- df$relvarwidth / max(df$relvarwidth)
      df$xmin <- df$x - df$relvarwidth * df$width / 2
      df$xmax <- df$x + df$relvarwidth * df$width / 2
    }
    df$width <- NULL
    if (!is.null(df$relvarwidth)) df$relvarwidth <- NULL
    
    df
  }
  
  draw <- function(., data, ..., fatten = 2, outlier.colour = NULL, outlier.shape = NULL, outlier.size = 2,
                   notch = FALSE, notchwidth = .5, varwidth = FALSE) {
    common <- data.frame(
      colour = data$colour,
      size = data$size,
      linetype = data$linetype,
      fill = alpha(data$fill, data$alpha),
      group = data$group,
      stringsAsFactors = FALSE
    )
    
    whiskers <- data.frame(
      x = data$x,
      xend = data$x,
      y = c(data$upper, data$lower),
      yend = c(data$ymax, data$ymin),
      alpha = NA,
      common)
    
    box <- data.frame(
      xmin = data$xmin,
      xmax = data$xmax,
      ymin = data$lower,
      y = data$middle,
      ymax = data$upper,
      ynotchlower = ifelse(notch, data$notchlower, NA),
      ynotchupper = ifelse(notch, data$notchupper, NA),
      notchwidth = notchwidth,
      alpha = data$alpha,
      common)
    
    #  if (!is.null(data$outliers) && length(data$outliers[[1]] >= 1)) {
    #    outliers <- data.frame(
    #      y = data$outliers[[1]],
    #      x = data$x[1],
    #      colour = outlier.colour %||% data$colour[1],
    #      shape = outlier.shape %||% data$shape[1],
    #      size = outlier.size %||% data$size[1],
    #      fill = NA,
    #      alpha = NA,
    #      stringsAsFactors = FALSE)
    #    outliers_grob <- GeomPoint$draw(outliers, ...)
    #  } else {
    outliers_grob <- NULL
    #  }
    
    ggname(.$my_name(), grobTree(
      outliers_grob,
      GeomSegment$draw(whiskers, ...),
      GeomCrossbar$draw(box, fatten = fatten, ...)
    ))
  }
  
  guide_geom <- function(.) "boxplot_noOutliers"
  draw_legend <- function(., data, ...)  {
    data <- aesdefaults(data, .$default_aes(), list(...))
    gp <- with(data, gpar(col=colour, fill=alpha(fill, alpha), lwd=size * .pt, lty = linetype))
    gTree(gp = gp, children = gList(
      linesGrob(0.5, c(0.1, 0.25)),
      linesGrob(0.5, c(0.75, 0.9)),
      rectGrob(height=0.5, width=0.75),
      linesGrob(c(0.125, 0.875), 0.5)
    ))
  }
  
  default_stat <- function(.) StatBoxplot
  default_pos <- function(.) PositionDodge
  default_aes <- function(.) aes(weight=1, colour="grey20", fill="white", size=0.5, alpha = NA, shape = 16, linetype = "solid")
  required_aes <- c("x", "lower", "upper", "middle", "ymin", "ymax")
  
})