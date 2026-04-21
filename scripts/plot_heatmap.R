#!/usr/bin/env Rscript
# =============================================================================
# plot_heatmap.R  —  Visualisation script for MetalGenie-Evo output
# =============================================================================
# Generates hierarchically clustered heatmaps from MetalGenie-Evo CSV output.
# Automatically detects the input type (gene counts vs coverage) and adapts
# colour scales accordingly.
#
# Usage
# -----
#   Rscript scripts/plot_heatmap.R --input results/MetalGenie-Evo-heatmap-data.csv
#   Rscript scripts/plot_heatmap.R --input results/ --type both --format both
#
# Arguments
# ---------
#   --input   Path to a single CSV file OR a MetalGenie-Evo output directory.
#             If a directory is given, all *heatmap*.csv files are processed.
#   --type    Plot type: "static" | "interactive" | "both"  (default: both)
#   --format  Output format for static: "pdf" | "png" | "both"  (default: both)
#   --out     Output directory for plots (default: same dir as input CSV)
#   --width   Plot width in inches (default: auto)
#   --height  Plot height in inches (default: auto)
#   --min_count Minimum gene count/coverage to include a category row (default: 0)
#   --help    Show this help message
#
# Required R packages
# -------------------
#   ggplot2, pheatmap, plotly, htmlwidgets, optparse, RColorBrewer, scales
#
#   Install with:
#     conda install -c conda-forge r-ggplot2 r-pheatmap r-plotly r-htmlwidgets
#                                  r-optparse r-rcolorbrewer r-scales
#   or:
#     Rscript -e 'install.packages(c("ggplot2","pheatmap","plotly","htmlwidgets",
#                                    "optparse","RColorBrewer","scales"))'
# =============================================================================

suppressPackageStartupMessages({
  if (!require("optparse",     quietly = TRUE)) stop("Package 'optparse' required.")
  if (!require("ggplot2",      quietly = TRUE)) stop("Package 'ggplot2' required.")
  if (!require("pheatmap",     quietly = TRUE)) stop("Package 'pheatmap' required.")
  if (!require("plotly",       quietly = TRUE)) stop("Package 'plotly' required.")
  if (!require("htmlwidgets",  quietly = TRUE)) stop("Package 'htmlwidgets' required.")
  if (!require("RColorBrewer", quietly = TRUE)) stop("Package 'RColorBrewer' required.")
  if (!require("scales",       quietly = TRUE)) stop("Package 'scales' required.")
})

# ─────────────────────────────────────────────────────────────────────────────
# Argument parsing
# ─────────────────────────────────────────────────────────────────────────────

option_list <- list(
  make_option("--input",     type = "character", default = NULL,
              help = "CSV file or MetalGenie-Evo output directory"),
  make_option("--type",      type = "character", default = "both",
              help = "Plot type: static | interactive | both  [default: both]"),
  make_option("--format",    type = "character", default = "both",
              help = "Static format: pdf | png | both  [default: both]"),
  make_option("--out",       type = "character", default = NULL,
              help = "Output directory for plots [default: same as input]"),
  make_option("--width",     type = "double",    default = NULL,
              help = "Plot width in inches [default: auto]"),
  make_option("--height",    type = "double",    default = NULL,
              help = "Plot height in inches [default: auto]"),
  make_option("--min_count", type = "double",    default = 0,
              help = "Min gene count/coverage per category to include [default: 0]"),
  make_option("--no_cluster_rows", action = "store_true", default = FALSE,
              help = "Disable hierarchical clustering of categories"),
  make_option("--no_cluster_cols", action = "store_true", default = FALSE,
              help = "Disable hierarchical clustering of genomes")
)

parser <- OptionParser(
  usage       = "Rscript scripts/plot_heatmap.R --input <file_or_dir> [options]",
  option_list = option_list
)
opt <- parse_args(parser)

if (is.null(opt$input)) {
  print_help(parser)
  stop("--input is required.")
}

# ─────────────────────────────────────────────────────────────────────────────
# Collect input CSV files
# ─────────────────────────────────────────────────────────────────────────────

collect_csvs <- function(path) {
  if (file.info(path)$isdir) {
    files <- list.files(path, pattern = "heatmap.*\\.csv$",
                        full.names = TRUE, ignore.case = TRUE)
    # Also grab the comment-header coverage file
    cov <- list.files(path, pattern = "coverage-heatmap.*\\.csv$",
                      full.names = TRUE, ignore.case = TRUE)
    all_files <- unique(c(files, cov))
    if (length(all_files) == 0)
      stop("No *heatmap*.csv files found in: ", path)
    return(all_files)
  }
  if (!file.exists(path)) stop("File not found: ", path)
  return(path)
}

csv_files <- collect_csvs(opt$input)
cat(sprintf("[INFO] Found %d heatmap file(s)\n", length(csv_files)))

# ─────────────────────────────────────────────────────────────────────────────
# Data loading and preprocessing
# ─────────────────────────────────────────────────────────────────────────────

detect_type <- function(filepath) {
  # Peek at comment lines for coverage metric label
  lines <- readLines(filepath, n = 3)
  if (any(grepl("^#.*coverage metric", lines, ignore.case = TRUE)))
    return("coverage")
  if (any(grepl("^#", lines)))
    return("coverage")   # any comment header → coverage file
  return("count")
}

load_heatmap_csv <- function(filepath) {
  # Skip comment lines starting with #
  lines <- readLines(filepath)
  data_lines <- lines[!grepl("^#", lines)]
  con <- textConnection(data_lines)
  on.exit(close(con))
  df <- read.csv(con, check.names = FALSE, row.names = 1)
  df
}

clean_matrix <- function(df, min_count = 0) {
  # Remove rows where all values <= min_count
  keep <- apply(df, 1, function(x) max(x, na.rm = TRUE) > min_count)
  df <- df[keep, , drop = FALSE]
  if (nrow(df) == 0) stop("No rows remain after filtering (--min_count too high?)")
  # Remove columns that are all zero
  keep_col <- apply(df, 2, function(x) max(x, na.rm = TRUE) > 0)
  df[, keep_col, drop = FALSE]
}

pretty_category <- function(x) {
  x <- gsub("_", " ", x)
  x <- gsub("-", " — ", x)
  x <- tools::toTitleCase(x)
  x
}

# ─────────────────────────────────────────────────────────────────────────────
# Colour palettes
# ─────────────────────────────────────────────────────────────────────────────

palette_for_type <- function(type) {
  if (type == "coverage") {
    # Blue-purple gradient for coverage/TPM
    colorRampPalette(c("#F7FBFF", "#C6DBEF", "#6BAED6", "#2171B5", "#08306B"))(100)
  } else {
    # White → warm orange → dark red for gene counts
    colorRampPalette(c("#FFFFFF", "#FFF3E0", "#FF8F00", "#BF360C"))(100)
  }
}

# ─────────────────────────────────────────────────────────────────────────────
# Static heatmap (pheatmap)
# ─────────────────────────────────────────────────────────────────────────────

make_static <- function(mat, title, type, out_base,
                        formats = "both",
                        width = NULL, height = NULL,
                        cluster_rows = TRUE, cluster_cols = TRUE) {

  pal   <- palette_for_type(type)
  nrow_ <- nrow(mat)
  ncol_ <- ncol(mat)

  # Auto-size
  w <- if (!is.null(width))  width  else max(8,  2 + ncol_ * 0.5)
  h <- if (!is.null(height)) height else max(6,  2 + nrow_ * 0.35)

  # Row labels — pretty-print category names
  rownames(mat) <- pretty_category(rownames(mat))

  # Legend label
  legend_title <- if (type == "coverage") "TPM / Depth" else "Gene count"

  # Value breaks: log1p scaling for counts, linear for coverage
  if (type == "count" && max(mat) > 10) {
    mat_display <- log1p(mat)
    legend_title <- "log(count + 1)"
  } else {
    mat_display <- mat
  }

  save_pheatmap <- function(filename) {
    pheatmap::pheatmap(
      mat_display,
      color            = pal,
      cluster_rows     = cluster_rows && nrow_ > 1,
      cluster_cols     = cluster_cols && ncol_ > 1,
      clustering_method = "ward.D2",
      border_color     = "white",
      cellwidth        = max(12, min(30, 400 / ncol_)),
      cellheight       = max(10, min(20, 300 / nrow_)),
      fontsize_row     = 9,
      fontsize_col     = 8,
      angle_col        = 45,
      main             = title,
      legend           = TRUE,
      filename         = filename,
      width            = w,
      height           = h,
      silent           = TRUE
    )
  }

  if (formats %in% c("pdf", "both")) {
    f <- paste0(out_base, ".pdf")
    save_pheatmap(f)
    cat(sprintf("[INFO]   Static PDF  → %s\n", f))
  }
  if (formats %in% c("png", "both")) {
    f <- paste0(out_base, ".png")
    save_pheatmap(f)
    cat(sprintf("[INFO]   Static PNG  → %s\n", f))
  }
}

# ─────────────────────────────────────────────────────────────────────────────
# Interactive heatmap (plotly)
# ─────────────────────────────────────────────────────────────────────────────

make_interactive <- function(mat, title, type, out_base,
                             cluster_rows = TRUE, cluster_cols = TRUE) {

  # Apply clustering to get dendrogram order
  if (cluster_rows && nrow(mat) > 1) {
    hc_row <- hclust(dist(mat),    method = "ward.D2")
    mat    <- mat[hc_row$order, , drop = FALSE]
  }
  if (cluster_cols && ncol(mat) > 1) {
    hc_col <- hclust(dist(t(mat)), method = "ward.D2")
    mat    <- mat[, hc_col$order,  drop = FALSE]
  }

  # Display values
  if (type == "count" && max(mat) > 10) {
    display_mat <- log1p(mat)
    cbar_title  <- "log(n + 1)"
  } else {
    display_mat <- mat
    cbar_title  <- if (type == "coverage") "TPM / Depth" else "Gene count"
  }

  colorscale <- if (type == "coverage") {
    list(c(0, "#F7FBFF"), c(0.25, "#C6DBEF"), c(0.5, "#6BAED6"),
         c(0.75, "#2171B5"), c(1, "#08306B"))
  } else {
    list(c(0, "#FFFFFF"), c(0.33, "#FFF3E0"), c(0.66, "#FF8F00"),
         c(1, "#BF360C"))
  }

  # Hover text — show original (un-transformed) values
  hover_text <- matrix(
    sprintf("<b>%s</b><br>%s<br>Value: %s",
            rep(pretty_category(rownames(mat)), ncol(mat)),
            rep(colnames(mat), each = nrow(mat)),
            format(as.vector(mat), digits = 3, big.mark = ",")),
    nrow = nrow(mat), ncol = ncol(mat)
  )

  fig <- plotly::plot_ly(
    z           = display_mat,
    x           = colnames(mat),
    y           = pretty_category(rownames(mat)),
    type        = "heatmap",
    colorscale  = colorscale,
    text        = hover_text,
    hoverinfo   = "text",
    colorbar    = list(title = cbar_title)
  ) %>%
    plotly::layout(
      title  = list(text = title, font = list(size = 14)),
      xaxis  = list(title = "", tickangle = -45,
                    tickfont = list(size = 10)),
      yaxis  = list(title = "", autorange = "reversed",
                    tickfont = list(size = 9)),
      margin = list(l = 200, b = 120, t = 60, r = 60)
    ) %>%
    plotly::config(displayModeBar = TRUE,
                   toImageButtonOptions = list(
                     format = "png", filename = basename(out_base),
                     width = 1200, height = 800
                   ))

  out_file <- paste0(out_base, ".html")
  htmlwidgets::saveWidget(fig, file = out_file, selfcontained = TRUE)
  cat(sprintf("[INFO]   Interactive → %s\n", out_file))
}

# ─────────────────────────────────────────────────────────────────────────────
# Main loop
# ─────────────────────────────────────────────────────────────────────────────

for (csv_path in csv_files) {

  cat(sprintf("\n[INFO] Processing: %s\n", basename(csv_path)))

  # Detect data type
  data_type <- detect_type(csv_path)
  cat(sprintf("[INFO]   Detected type: %s\n", data_type))

  # Load and clean
  df  <- tryCatch(load_heatmap_csv(csv_path),
                  error = function(e) { warning("Could not load ", csv_path, ": ", e$message); NULL })
  if (is.null(df)) next

  mat <- tryCatch(clean_matrix(as.matrix(df), min_count = opt$min_count),
                  error = function(e) { warning(e$message); NULL })
  if (is.null(mat)) next

  cat(sprintf("[INFO]   Matrix: %d categories × %d genomes\n", nrow(mat), ncol(mat)))

  # Determine output directory and base name
  out_dir <- if (!is.null(opt$out)) opt$out else dirname(csv_path)
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

  stem     <- tools::file_path_sans_ext(basename(csv_path))
  out_base <- file.path(out_dir, stem)

  # Title
  title <- if (data_type == "coverage") {
    "MetalGenie-Evo — Coverage-based abundance"
  } else {
    "MetalGenie-Evo — Gene presence by category"
  }

  cluster_rows <- !opt$no_cluster_rows
  cluster_cols <- !opt$no_cluster_cols

  # Generate plots
  if (opt$type %in% c("static", "both")) {
    tryCatch(
      make_static(mat, title, data_type, out_base,
                  formats      = opt$format,
                  width        = opt$width,
                  height       = opt$height,
                  cluster_rows = cluster_rows,
                  cluster_cols = cluster_cols),
      error = function(e) warning("Static plot failed: ", e$message)
    )
  }

  if (opt$type %in% c("interactive", "both")) {
    tryCatch(
      make_interactive(mat, title, data_type, out_base,
                       cluster_rows = cluster_rows,
                       cluster_cols = cluster_cols),
      error = function(e) warning("Interactive plot failed: ", e$message)
    )
  }
}

cat("\n[DONE] All plots written.\n")
