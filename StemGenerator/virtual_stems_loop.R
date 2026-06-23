# ============================================================
# virtual_stem_manual_function_loop_las.R
# ============================================================
# Loop script for virtual_stem_stepwise_manual_las_REVISED.R
#
# This script does NOT modify the manual stepwise script.
# It reads only the function definitions from that file, then runs
# the same steps repeatedly with different seeds.
# ============================================================

# ----------------------------
# 0) USER SETTINGS
# ----------------------------

# Use the current R working directory as project/script directory.
# In RStudio, set it to the folder containing this script and
# virtual_stem_stepwise_manual_las_REVISED.R before running.

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

functions_script <- file.path(script_dir, "virtual_stem_stepwise_w_functions.R")

output_dir <- file.path(script_dir, "virtual_las_files")
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

summary_csv <- file.path(output_dir, "virtual_las_summary.csv")
summary_rds <- file.path(output_dir, "virtual_las_summary.rds")

n_iterations <- 100L
seed_start   <- 1L

write_las <- TRUE
overwrite_las <- FALSE
return_stem_objects <- FALSE

# Stem geometry / scanner settings used for every iteration.
# Random shape variation comes from changing the seed.
height     <- 0.30
radius     <- 0.20
resolution <- 0.010
center_x   <- 0
center_y   <- 0
z0         <- 0

curvature  <- 20
elasticity <- 10

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
gap_size <- 0.40

# Diameter / target-slice summary settings.
# This mirrors the old RDS target idea: target radius = configured base radius,
# target center = local centerline around the middle stem slice.
slice_half_height <- 0.15

# ----------------------------
# 1) Load only functions from the manual script
# ----------------------------

source_only_function_definitions <- function(file, envir = new.env(parent = globalenv())) {
  if (!file.exists(file)) {
    stop("Function script not found: ", file)
  }

  exprs <- parse(file = file)

  for (expr in exprs) {
    is_assignment <- is.call(expr) && as.character(expr[[1]]) %in% c("<-", "=")

    if (is_assignment) {
      rhs <- expr[[3]]
      rhs_is_function <- is.call(rhs) && identical(rhs[[1]], as.name("function"))

      if (rhs_is_function) {
        eval(expr, envir = envir)
      }
    }
  }

  envir
}

function_env <- source_only_function_definitions(functions_script)

create_ideal_cylinder    <- function_env$create_ideal_cylinder
generate_centerline      <- function_env$generate_centerline
apply_centerline         <- function_env$apply_centerline
add_column_dents         <- function_env$add_column_dents
add_point_imperfections  <- function_env$add_point_imperfections
add_scanner_noise        <- function_env$add_scanner_noise
omit_scan_sections       <- function_env$omit_scan_sections
stem_to_las              <- function_env$stem_to_las

# ----------------------------
# 2) Small loop helpers
# ----------------------------

fmt_num <- function(x) {
  x <- formatC(as.numeric(x), format = "f", digits = 3)
  x <- sub("0+$", "", x)
  x <- sub("\\.$", "", x)
  gsub("\\.", "p", x)
}

make_las_file_name <- function(iteration, seed, radius, height, resolution, noise, curvature) {
  paste0(
    "virtual_stem_",
    sprintf("%05d", as.integer(iteration)),
    "_seed-", sprintf("%08d", as.integer(seed)),
    "_h-", fmt_num(height),
    "_r-", fmt_num(radius),
    "_res-", fmt_num(resolution),
    "_noise-", fmt_num(noise),
    "_curv-", fmt_num(curvature),
    ".las"
  )
}

summarise_virtual_stem <- function(
    iteration,
    seed,
    stem_before_omission,
    stem_final,
    output_file
) {
  z_center <- z0 + height / 2

  slice <- stem_final[
    stem_final$Z >= (z_center - slice_half_height) &
      stem_final$Z <  (z_center + slice_half_height),
    , drop = FALSE
  ]

  if (nrow(slice) > 0L) {
    target_cx <- center_x + stats::median(slice$centerline_x, na.rm = TRUE)
    target_cy <- center_y + stats::median(slice$centerline_y, na.rm = TRUE)

    observed_radius <- sqrt((slice$X - target_cx)^2 + (slice$Y - target_cy)^2)
  } else {
    target_cx <- NA_real_
    target_cy <- NA_real_
    observed_radius <- numeric(0)
  }

  n_points_before_omission <- nrow(stem_before_omission)
  n_points_remaining <- nrow(stem_final)
  n_points_omitted <- n_points_before_omission - n_points_remaining

  data.frame(
    iteration = as.integer(iteration),
    seed = as.integer(seed),

    # Main annotation / target values
    target_cx_m = target_cx,
    target_cy_m = target_cy,
    target_radius_m = radius,
    target_diameter_m = 2 * radius,

    # Observed diameter estimates from the retained points in the target slice
    observed_radius_mean_m = if (length(observed_radius) > 0L) mean(observed_radius, na.rm = TRUE) else NA_real_,
    observed_radius_median_m = if (length(observed_radius) > 0L) stats::median(observed_radius, na.rm = TRUE) else NA_real_,
    observed_radius_min_m = if (length(observed_radius) > 0L) min(observed_radius, na.rm = TRUE) else NA_real_,
    observed_radius_max_m = if (length(observed_radius) > 0L) max(observed_radius, na.rm = TRUE) else NA_real_,
    observed_diameter_mean_m = if (length(observed_radius) > 0L) 2 * mean(observed_radius, na.rm = TRUE) else NA_real_,
    observed_diameter_median_m = if (length(observed_radius) > 0L) 2 * stats::median(observed_radius, na.rm = TRUE) else NA_real_,

    # Input parameters
    height_m = height,
    base_radius_m = radius,
    base_diameter_m = 2 * radius,
    resolution_m = resolution,
    center_x_m = center_x,
    center_y_m = center_y,
    z0_m = z0,
    curvature_cm_per_m = curvature,
    elasticity_m = elasticity,
    dent_chance = dent_chance,
    dent_depth_m = dent_depth,
    dent_width_columns = dent_width,
    min_radius_fraction = min_radius_fraction,
    imperfection_quantity = imperfection_quantity,
    imperfection_depth_m = imperfection_depth,
    imperfection_elasticity_cells = imperfection_elasticity,
    noise_m = noise,
    noise_bias_m = noise_bias,
    tree_scanned_completely = tree_scanned_completely,
    omitted_sections = omitted_sections,
    gap_size_fraction = gap_size,

    # Shape / point-cloud diagnostics
    n_ring = as.integer(attr(stem_final, "n_ring")),
    n_z = as.integer(attr(stem_final, "n_z")),
    n_points_before_omission = n_points_before_omission,
    n_points_remaining = n_points_remaining,
    n_points_omitted = n_points_omitted,
    omission_fraction = n_points_omitted / n_points_before_omission,
    n_points_target_slice = nrow(slice),
    max_centerline_shift_m = max(sqrt(stem_before_omission$centerline_x^2 + stem_before_omission$centerline_y^2), na.rm = TRUE),
    mean_current_radius_m = mean(stem_before_omission$current_radius, na.rm = TRUE),
    median_current_radius_m = stats::median(stem_before_omission$current_radius, na.rm = TRUE),
    min_current_radius_m = min(stem_before_omission$current_radius, na.rm = TRUE),
    max_current_radius_m = max(stem_before_omission$current_radius, na.rm = TRUE),
    max_abs_dent_offset_m = max(abs(stem_before_omission$dent_offset), na.rm = TRUE),
    max_abs_imperfection_offset_m = max(abs(stem_before_omission$imperfection_offset), na.rm = TRUE),
    mean_radial_noise_m = mean(stem_before_omission$radial_noise, na.rm = TRUE),
    sd_radial_noise_m = stats::sd(stem_before_omission$radial_noise, na.rm = TRUE),

    las_file = output_file,
    generation_time = format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
    stringsAsFactors = FALSE
  )
}

run_one_iteration <- function(iteration, seed) {
  set.seed(seed)

  stem01 <- create_ideal_cylinder(
    height = height,
    radius = radius,
    resolution = resolution,
    center_x = center_x,
    center_y = center_y,
    z0 = z0
  )

  centerline02 <- generate_centerline(
    z_values = sort(unique(stem01$z_base)),
    curvature = curvature,
    elasticity = elasticity
  )

  stem03 <- apply_centerline(stem01, centerline02)

  stem04 <- add_column_dents(
    stem = stem03,
    dent_chance = dent_chance,
    dent_depth = dent_depth,
    dent_width = dent_width,
    min_radius_fraction = min_radius_fraction
  )

  stem05 <- add_point_imperfections(
    stem = stem04,
    imperfection_quantity = imperfection_quantity,
    imperfection_depth = imperfection_depth,
    imperfection_elasticity = imperfection_elasticity
  )

  stem06 <- add_scanner_noise(
    stem = stem05,
    noise = noise,
    noise_bias = noise_bias
  )

  stem07 <- omit_scan_sections(
    stem = stem06,
    omitted_sections = omitted_sections,
    gap_size = gap_size,
    tree_scanned_completely = tree_scanned_completely
  )

  las_file <- file.path(
    output_dir,
    make_las_file_name(iteration, seed, radius, height, resolution, noise, curvature)
  )

  if (file.exists(las_file)) {
    if (isTRUE(overwrite_las)) {
      file.remove(las_file)
    } else {
      stop("LAS file already exists and overwrite_las is FALSE: ", las_file)
    }
  }

  las <- NULL
  if (isTRUE(write_las)) {
    las <- stem_to_las(stem07, output_file = las_file)
  }

  summary_row <- summarise_virtual_stem(
    iteration = iteration,
    seed = seed,
    stem_before_omission = stem06,
    stem_final = stem07,
    output_file = if (isTRUE(write_las)) las_file else NA_character_
  )

  list(
    summary = summary_row,
    stem = if (isTRUE(return_stem_objects)) stem07 else NULL,
    las = if (isTRUE(return_stem_objects)) las else NULL
  )
}

# ----------------------------
# 3) Run loop
# ----------------------------

loop_results <- vector("list", n_iterations)
summary_rows <- vector("list", n_iterations)

for (i in seq_len(n_iterations)) {
  seed_i <- seed_start + i - 1L

  message("[", i, "/", n_iterations, "] seed = ", seed_i)

  loop_results[[i]] <- run_one_iteration(iteration = i, seed = seed_i)
  summary_rows[[i]] <- loop_results[[i]]$summary

  if (!isTRUE(return_stem_objects)) {
    loop_results[[i]] <- NULL
    gc(verbose = FALSE)
  }
}

summary_table <- do.call(rbind, summary_rows)
rownames(summary_table) <- NULL

utils::write.csv(summary_table, summary_csv, row.names = FALSE)
saveRDS(summary_table, summary_rds)

message("Done.")
message("LAS output directory: ", output_dir)
message("Summary CSV: ", summary_csv)
message("Summary RDS: ", summary_rds)

# Final objects in R:
#   summary_table
#   loop_results   only contains objects when return_stem_objects = TRUE
