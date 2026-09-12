# theme_dark_bw.R
# Black-background, white-ink ggplot theme for the EEG onset project.
#
# Presentation only: nothing here changes data generation or statistical
# analysis.
#
# Design rationale
# ----------------
# Black panel, thin white panel border, faint grid, white data ink, and legend
# keys long enough to judge a line type. Series are separated by line type,
# line weight, point shape, faceting and direct labelling rather than by
# colour.
#
# Grey is used in exactly one role: context traces are drawn in light grey so
# that a single focal trace can be superimposed in thicker white. Grey never
# encodes a category anywhere else.

# ---- Shared ink constants ---------------------------------------------------
INK        <- "#FFFFFF"    # all data ink and text
PAPER      <- "#000000"    # all backgrounds
GRID_COL   <- "#4D4D4D"    # furniture only, never encodes data
CONTEXT    <- "#BFBFBF"    # context traces behind a focal white trace

# Line weights in mm (ggplot2 linewidth unit). Data lines are kept well above
# the 0.5 pt (~0.18 mm) minimum used for print legibility.
LW_AXIS    <- 0.4
LW_GRID    <- 0.25
LW_DATA    <- 0.6
LW_CONTEXT <- 0.32
LW_FOCAL   <- 1.05

theme_dark_bw <- function(base_size = 13, grid = c("both", "y", "x", "none")){

  grid <- match.arg(grid)

  grid_major <- switch(
    grid,
    both = ggplot2::element_line(colour = GRID_COL, linewidth = LW_GRID),
    none = ggplot2::element_blank(),
    ggplot2::element_blank()
  )

  th <- ggplot2::theme_bw(base_size = base_size) +
    ggplot2::theme(
      # --- paper -------------------------------------------------------------
      plot.background   = ggplot2::element_rect(fill = PAPER, colour = NA),
      panel.background  = ggplot2::element_rect(fill = PAPER, colour = NA),

      # --- furniture ---------------------------------------------------------
      panel.grid.major  = grid_major,
      panel.grid.minor  = ggplot2::element_blank(),
      panel.border      = ggplot2::element_rect(
        fill = NA, colour = INK, linewidth = LW_AXIS
      ),
      axis.ticks        = ggplot2::element_line(colour = INK, linewidth = LW_AXIS),

      # --- text --------------------------------------------------------------
      text              = ggplot2::element_text(colour = INK),
      axis.text         = ggplot2::element_text(colour = INK),
      axis.title        = ggplot2::element_text(colour = INK),
      plot.title        = ggplot2::element_text(colour = INK, face = "bold"),
      plot.subtitle     = ggplot2::element_text(colour = INK),
      plot.caption      = ggplot2::element_text(colour = INK, hjust = 0),

      # --- legend: long keys, so line types can actually be judged ------------
      legend.background = ggplot2::element_rect(fill = PAPER, colour = NA),
      legend.key        = ggplot2::element_rect(fill = PAPER, colour = NA),
      legend.key.width  = grid::unit(1.5, "cm"),
      legend.text       = ggplot2::element_text(colour = INK),
      legend.title      = ggplot2::element_text(colour = INK),

      # --- facets: white strip, thin black rule ------------------------------
      strip.background  = ggplot2::element_rect(
        fill = PAPER, colour = INK, linewidth = LW_AXIS
      ),
      strip.text        = ggplot2::element_text(colour = INK, face = "bold")
    )

  if (grid == "y") {
    th <- th + ggplot2::theme(
      panel.grid.major.y = ggplot2::element_line(colour = GRID_COL, linewidth = LW_GRID),
      panel.grid.major.x = ggplot2::element_blank()
    )
  }
  if (grid == "x") {
    th <- th + ggplot2::theme(
      panel.grid.major.x = ggplot2::element_line(colour = GRID_COL, linewidth = LW_GRID),
      panel.grid.major.y = ggplot2::element_blank()
    )
  }
  th
}

# ---- Global geom ink --------------------------------------------------------
# Set the default drawing colour once so that no layer has to repeat it. This
# keeps drawing colour centralised so layers do not need repeated colour settings.
set_bw_geom_defaults <- function(){
  ggplot2::update_geom_defaults("line",     list(colour = INK))
  ggplot2::update_geom_defaults("path",     list(colour = INK))
  ggplot2::update_geom_defaults("point",    list(colour = INK))
  ggplot2::update_geom_defaults("text",     list(colour = INK))
  ggplot2::update_geom_defaults("label",    list(colour = INK, fill = PAPER))
  ggplot2::update_geom_defaults("segment",  list(colour = INK))
  ggplot2::update_geom_defaults("density",  list(colour = INK, fill = NA))
  ggplot2::update_geom_defaults("hline",    list(colour = INK))
  ggplot2::update_geom_defaults("vline",    list(colour = INK))
  ggplot2::update_geom_defaults("crossbar", list(colour = INK, fill = NA))
  ggplot2::update_geom_defaults("boxplot",  list(colour = INK, fill = NA))
  ggplot2::update_geom_defaults("rect",     list(colour = INK, fill = NA))
  invisible(TRUE)
}

# ---- Monochrome scales ------------------------------------------------------
# Ordered from most to least distinguishable at small print sizes. "dotted" is
# kept last and should be reserved for reference lines, not data series.
mono_linetypes <- c("solid", "longdash", "dashed", "dotdash", "twodash", "dotted")

# Filled and open forms alternate so that neighbouring levels never share a
# silhouette: circle, triangle, square, diamond, open circle, open triangle.
mono_shapes <- c(16, 17, 15, 18, 1, 2, 0, 5, 6)

scale_linetype_bw <- function(...){
  ggplot2::scale_linetype_manual(values = mono_linetypes, na.value = "solid", ...)
}

scale_shape_bw <- function(...){
  ggplot2::scale_shape_manual(values = mono_shapes, ...)
}

scale_colour_bw <- function(...){
  ggplot2::scale_colour_manual(values = rep(INK, 32), na.value = INK, ...)
}

scale_color_bw <- scale_colour_bw

scale_fill_bw <- function(...){
  ggplot2::scale_fill_manual(values = rep(PAPER, 32), na.value = PAPER, ...)
}
