# ============================================================
# virtual_stem_manual_linear_las.R
# ============================================================
# Manual/debug script.
# Structure is intentionally simple:
#   function -> call function -> plot result -> next function -> call -> plot
# ============================================================

# ----------------------------
# 0) PARAMETERS / DIRECTORIES
# ----------------------------

get_this_script_dir <- function() {
  # Case 1: script was run with source("path/to/script.R")
  for (i in rev(seq_len(sys.nframe()))) {
    this_file <- sys.frame(i)$ofile
    if (!is.null(this_file)) {
      return(dirname(normalizePath(this_file, winslash = "/", mustWork = TRUE)))
    }
  }
  
  # Case 2: script was run with Rscript
  args <- commandArgs(trailingOnly = FALSE)
  file_arg <- "--file="
  file_path <- sub(file_arg, "", args[startsWith(args, file_arg)])
  
  if (length(file_path) > 0) {
    return(dirname(normalizePath(file_path, winslash = "/", mustWork = TRUE)))
  }
  
  # Case 3: RStudio source editor
  if (requireNamespace("rstudioapi", quietly = TRUE) && rstudioapi::isAvailable()) {
    file_path <- rstudioapi::getSourceEditorContext()$path
    
    if (!is.null(file_path) && nzchar(file_path)) {
      return(dirname(normalizePath(file_path, winslash = "/", mustWork = TRUE)))
    }
  }
  
  stop(
    "Could not determine the script directory. ",
    "Please set SIMULATE_GROWTH_FILE manually.",
    call. = FALSE
  )
}

script_dir <- get_this_script_dir()

las_output_dir <- file.path(script_dir, "virtual_las_files")

if (!dir.exists(las_output_dir)) {
  dir.create(las_output_dir, recursive = TRUE)
}

output_dir <- las_output_dir

output_file <- file.path(output_dir, "manual_virtual_stem.las")
write_las   <- TRUE

seed <- 1

height     <- 0.30
radius     <- 0.20
resolution <- 0.010
center_x   <- 0
center_y   <- 0
z0         <- 0

curvature  <- 20      # cm lateral drift per m height
elasticity <- 10      # larger = smoother centerline

dent_chance <- 0.25
dent_depth  <- 0.030
dent_width  <- 5
min_radius_fraction <- 0.35

imperfection_quantity   <- 0.05
imperfection_depth      <- 0.10
imperfection_elasticity <- 3

noise      <- 0.010
noise_bias <- 0.000

tree_scanned_completely <- FALSE
omitted_sections <- 1
gap_size <- 0.40       # fraction of circumference affected by omission

max_plot_points <- 50000

set.seed(seed)

# ----------------------------
# Tiny plotting helpers
# ----------------------------

plot_stem <- function(stem, title = "stem", color = NULL) {
  if (!requireNamespace("plotly", quietly = TRUE)) {
    stop("Install plotly first: install.packages('plotly')")
  }

  x <- stem
  if (nrow(x) > max_plot_points) x <- x[sample(seq_len(nrow(x)), max_plot_points), ]

  marker <- list(size = 2)
  if (!is.null(color) && color %in% names(x)) {
    marker$color <- x[[color]]
    marker$colorscale <- "Viridis"
    marker$showscale <- TRUE
  }

  print(
    plotly::plot_ly(
      x, x = ~X, y = ~Y, z = ~Z,
      type = "scatter3d", mode = "markers",
      marker = marker
    ) |>
      plotly::layout(title = title, scene = list(aspectmode = "data"))
  )
}

plot_centerline <- function(centerline, title = "centerline") {
  if (!requireNamespace("plotly", quietly = TRUE)) {
    stop("Install plotly first: install.packages('plotly')")
  }

  print(
    plotly::plot_ly(
      centerline, x = ~cx, y = ~cy, z = ~z,
      type = "scatter3d", mode = "lines+markers",
      line = list(width = 6), marker = list(size = 3)
    ) |>
      plotly::layout(title = title, scene = list(aspectmode = "data"))
  )
}

# ============================================================
# 1) Ideal cylinder
# ============================================================

create_ideal_cylinder <- function(height, radius, resolution, center_x = 0, center_y = 0, z0 = 0) {
  n_ring <- max(8, round(2 * pi * radius / resolution))
  n_z <- max(2, round(height / resolution) + 1)

  theta <- seq(0, 2 * pi, length.out = n_ring + 1)[-1]
  z_values <- seq(z0, z0 + height, length.out = n_z)

  stem <- expand.grid(column_id = seq_len(n_ring), row_id = seq_len(n_z))
  stem$theta <- theta[stem$column_id]
  stem$z_base <- z_values[stem$row_id]

  stem$base_radius <- radius
  stem$current_radius <- radius

  stem$centerline_x <- 0
  stem$centerline_y <- 0
  stem$dent_offset <- 0
  stem$imperfection_offset <- 0
  stem$radial_noise <- 0
  stem$tangential_noise <- 0
  stem$vertical_noise <- 0
  stem$omitted <- FALSE

  stem$X <- center_x + stem$current_radius * cos(stem$theta)
  stem$Y <- center_y + stem$current_radius * sin(stem$theta)
  stem$Z <- stem$z_base

  attr(stem, "n_ring") <- n_ring
  attr(stem, "n_z") <- n_z
  attr(stem, "center_x") <- center_x
  attr(stem, "center_y") <- center_y

  stem
}

stem01 <- create_ideal_cylinder(height, radius, resolution, center_x, center_y, z0); plot_stem(stem01, "01 ideal cylinder")

# ============================================================
# 2) Centerline
# ============================================================

generate_centerline <- function(z_values, curvature = 0, elasticity = 1) {
  if (curvature == 0) {
    return(data.frame(z = z_values, cx = 0, cy = 0))
  }

  drift_per_m <- curvature / 100
  dz <- c(0, diff(z_values))

  knot_z <- unique(c(min(z_values), seq(min(z_values), max(z_values), by = elasticity), max(z_values)))
  if (length(knot_z) < 3) knot_z <- c(min(z_values), mean(range(z_values)), max(z_values))

  x_raw <- stats::splinefun(knot_z, stats::rnorm(length(knot_z)), method = "natural")(z_values)
  y_raw <- stats::splinefun(knot_z, stats::rnorm(length(knot_z)), method = "natural")(z_values)

  length_raw <- sqrt(x_raw^2 + y_raw^2)
  length_raw[length_raw == 0] <- 1

  data.frame(
    z = z_values,
    cx = cumsum((x_raw / length_raw) * drift_per_m * dz),
    cy = cumsum((y_raw / length_raw) * drift_per_m * dz)
  )
}

centerline02 <- generate_centerline(sort(unique(stem01$z_base)), curvature, elasticity); plot_centerline(centerline02, "02 centerline")

# ============================================================
# 3) Apply centerline to cylinder
# ============================================================

apply_centerline <- function(stem, centerline) {
  center_x <- attr(stem, "center_x")
  center_y <- attr(stem, "center_y")

  stem$centerline_x <- centerline$cx[stem$row_id]
  stem$centerline_y <- centerline$cy[stem$row_id]

  stem$X <- center_x + stem$centerline_x + stem$current_radius * cos(stem$theta) + stem$tangential_noise * (-sin(stem$theta))
  stem$Y <- center_y + stem$centerline_y + stem$current_radius * sin(stem$theta) + stem$tangential_noise * cos(stem$theta)
  stem$Z <- stem$z_base + stem$vertical_noise

  stem
}

stem03 <- apply_centerline(stem01, centerline02); plot_stem(stem03, "03 cylinder with centerline")

# ============================================================
# 4) Column dents
# ============================================================

add_column_dents <- function(stem, dent_chance, dent_depth, dent_width, min_radius_fraction) {
  n_ring <- attr(stem, "n_ring")
  center_x <- attr(stem, "center_x")
  center_y <- attr(stem, "center_y")

  dent_profile <- rep(0, n_ring)
  dent_centers <- which(stats::runif(n_ring) < dent_chance)

  for (dent_center in dent_centers) {
    dist <- pmin(abs(seq_len(n_ring) - dent_center), n_ring - abs(seq_len(n_ring) - dent_center))
    dent_profile <- dent_profile - abs(stats::rnorm(1, mean = dent_depth, sd = dent_depth / 3)) * exp(-0.5 * (dist / dent_width)^2)
  }

  dent_profile <- pmax(dent_profile, -(1 - min_radius_fraction) * stem$base_radius[1])

  stem$dent_offset <- dent_profile[stem$column_id]
  stem$current_radius <- stem$base_radius + stem$dent_offset + stem$imperfection_offset + stem$radial_noise

  stem$X <- center_x + stem$centerline_x + stem$current_radius * cos(stem$theta) + stem$tangential_noise * (-sin(stem$theta))
  stem$Y <- center_y + stem$centerline_y + stem$current_radius * sin(stem$theta) + stem$tangential_noise * cos(stem$theta)
  stem$Z <- stem$z_base + stem$vertical_noise

  stem
}

stem04 <- add_column_dents(stem03, dent_chance, dent_depth, dent_width, min_radius_fraction); plot_stem(stem04, "04 cylinder with column dents", color = "dent_offset")

# ============================================================
# 5) Point-level imperfections
# ============================================================

add_point_imperfections <- function(stem, imperfection_quantity, imperfection_depth, imperfection_elasticity) {
  n_ring <- attr(stem, "n_ring")
  n_z <- attr(stem, "n_z")
  center_x <- attr(stem, "center_x")
  center_y <- attr(stem, "center_y")

  raw <- stats::rnorm(nrow(stem), 0, imperfection_depth / 2) * (stats::runif(nrow(stem)) < imperfection_quantity)
  imp <- matrix(raw, nrow = n_z, ncol = n_ring, byrow = FALSE)

  if (imperfection_elasticity > 1) {
    smoothed <- imp
    w <- as.integer(imperfection_elasticity)
    for (i in seq_len(n_z)) {
      for (j in seq_len(n_ring)) {
        rows <- max(1, i - w):min(n_z, i + w)
        cols <- ((j - w):(j + w) - 1) %% n_ring + 1
        smoothed[i, j] <- mean(imp[rows, cols])
      }
    }
    imp <- smoothed
  }

  stem$imperfection_offset <- as.vector(imp)
  stem$current_radius <- stem$base_radius + stem$dent_offset + stem$imperfection_offset + stem$radial_noise

  stem$X <- center_x + stem$centerline_x + stem$current_radius * cos(stem$theta) + stem$tangential_noise * (-sin(stem$theta))
  stem$Y <- center_y + stem$centerline_y + stem$current_radius * sin(stem$theta) + stem$tangential_noise * cos(stem$theta)
  stem$Z <- stem$z_base + stem$vertical_noise

  stem
}

stem05 <- add_point_imperfections(stem04, imperfection_quantity, imperfection_depth, imperfection_elasticity); plot_stem(stem05, "05 cylinder with point imperfections", color = "imperfection_offset")

# ============================================================
# 6) Scanner noise
# ============================================================

add_scanner_noise <- function(stem, noise, noise_bias) {
  center_x <- attr(stem, "center_x")
  center_y <- attr(stem, "center_y")

  stem$radial_noise <- noise_bias + stats::rnorm(nrow(stem), 0, noise / 2)
  stem$tangential_noise <- stats::rnorm(nrow(stem), 0, noise / 4)
  stem$vertical_noise <- stats::rnorm(nrow(stem), 0, noise / 2)

  stem$current_radius <- stem$base_radius + stem$dent_offset + stem$imperfection_offset + stem$radial_noise

  stem$X <- center_x + stem$centerline_x + stem$current_radius * cos(stem$theta) + stem$tangential_noise * (-sin(stem$theta))
  stem$Y <- center_y + stem$centerline_y + stem$current_radius * sin(stem$theta) + stem$tangential_noise * cos(stem$theta)
  stem$Z <- stem$z_base + stem$vertical_noise

  stem
}

stem06 <- add_scanner_noise(stem05, noise, noise_bias); plot_stem(stem06, "06 cylinder with scanner noise", color = "radial_noise")

# ============================================================
# 7) Omit scan sections
# ============================================================

omit_scan_sections <- function(stem, omitted_sections, gap_size, tree_scanned_completely) {
  n_ring <- attr(stem, "n_ring")
  width <- max(1, round(gap_size * n_ring))
  remove_column <- rep(FALSE, n_ring)

  if (!tree_scanned_completely && omitted_sections < 1) omitted_sections <- 1

  if (omitted_sections > 0 && gap_size > 0) {
    for (i in seq_len(omitted_sections)) {
      center <- sample(seq_len(n_ring), 1)
      cols <- ((center - floor(width / 2)):(center + floor(width / 2)) - 1) %% n_ring + 1
      remove_column[cols] <- TRUE
    }
  }

  stem$omitted <- remove_column[stem$column_id]
  out <- stem[!stem$omitted, ]

  attr(out, "n_ring") <- attr(stem, "n_ring")
  attr(out, "n_z") <- attr(stem, "n_z")
  attr(out, "center_x") <- attr(stem, "center_x")
  attr(out, "center_y") <- attr(stem, "center_y")

  out
}

stem07 <- omit_scan_sections(stem06, omitted_sections, gap_size, tree_scanned_completely); plot_stem(stem07, "07 cylinder after omitted scan sections", color = "omitted")

# ============================================================
# 8) Convert to LAS object / optionally write LAS file
# ============================================================

stem_to_las <- function(stem, output_file = NULL) {
  if (!requireNamespace("lidR", quietly = TRUE)) {
    stop("Install lidR first: install.packages('lidR')")
  }

  las <- lidR::LAS(stem)

  if (!is.null(output_file)) {
    dir.create(dirname(output_file), recursive = TRUE, showWarnings = FALSE)
    lidR::writeLAS(las, output_file)
  }

  las
}

las08 <- stem_to_las(stem07, if (write_las) output_file else NULL)

# Final objects available in your R environment:
#   stem01, centerline02, stem03, stem04, stem05, stem06, stem07, las08
