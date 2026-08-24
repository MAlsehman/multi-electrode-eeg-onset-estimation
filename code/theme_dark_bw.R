# theme_dark_bw.R
# Strict black-and-white ggplot theme and monochrome plotting scales for the EEG onset project.
# This file affects presentation only and does not alter data generation or statistical analysis.

theme_dark_bw <- function(base_size = 13){
  ggplot2::theme_minimal(base_size = base_size) +
    ggplot2::theme(
      plot.background    = ggplot2::element_rect(fill = "#000000", colour = NA),
      panel.background   = ggplot2::element_rect(fill = "#000000", colour = NA),
      panel.grid.major   = ggplot2::element_blank(),
      panel.grid.minor   = ggplot2::element_blank(),
      panel.border       = ggplot2::element_rect(
        fill = NA, colour = "#ffffff", linewidth = 0.4
      ),
      text               = ggplot2::element_text(colour = "#ffffff"),
      axis.text          = ggplot2::element_text(colour = "#ffffff"),
      axis.title         = ggplot2::element_text(colour = "#ffffff"),
      axis.ticks         = ggplot2::element_line(colour = "#ffffff"),
      plot.title         = ggplot2::element_text(
        colour = "#ffffff", face = "bold"
      ),
      plot.subtitle      = ggplot2::element_text(colour = "#ffffff"),
      legend.background = ggplot2::element_rect(
        fill = "#000000", colour = NA
      ),
      legend.key         = ggplot2::element_rect(
        fill = "#000000", colour = NA
      ),
      legend.text        = ggplot2::element_text(colour = "#ffffff"),
      legend.title       = ggplot2::element_text(colour = "#ffffff"),
      strip.background   = ggplot2::element_rect(
        fill = "#000000", colour = "#ffffff", linewidth = 0.4
      ),
      strip.text         = ggplot2::element_text(
        colour = "#ffffff", face = "bold"
      )
    )
}

mono_linetypes <- c(
  "solid", "dashed", "dotted", "dotdash", "longdash", "twodash"
)

mono_shapes <- c(16, 17, 15, 18, 1, 2, 0, 5, 6)

scale_linetype_bw <- function(...){
  ggplot2::scale_linetype_manual(values = mono_linetypes, ...)
}

scale_shape_bw <- function(...){
  ggplot2::scale_shape_manual(values = mono_shapes, ...)
}

scale_colour_bw <- function(...){
  ggplot2::scale_colour_manual(values = rep("#ffffff", 32), ...)
}

scale_color_bw <- scale_colour_bw

scale_fill_bw <- function(...){
  ggplot2::scale_fill_manual(values = rep("#ffffff", 32), ...)
}