# ==============================================================================
# CHELSA-TraCE21k processing script
#
# Purpose:
#   1. Read text files containing URLs for 19 CHELSA bioclimatic rasters
#      (bio01 to bio19) for each time slice
#   2. Aggregate rasters to a coarser resolution
#   3. Convert aggregated rasters to data frames
#   4. Combine all time slices into a single matrix
#   5. Run a PCA on bio01-bio19 across all time slices together
#
# Input files:
#   bio0000.txt, bio1000.txt, ..., bio20000.txt
#
# Expected content of each text file:
#   Exactly 19 raster links, corresponding to:
#   bio01, bio02, bio03, bio04, bio05, bio06, bio07, bio08, bio09, bio10,
#   bio11, bio12, bio13, bio14, bio15, bio16, bio17, bio18, bio19
# ==============================================================================

# ------------------------------------------------------------------------------
# Load packages
# ------------------------------------------------------------------------------
library(terra)
library(dplyr)
library(purrr)

options(stringsAsFactors = FALSE)

# ------------------------------------------------------------------------------
# Input/output directories
#
# The script assumes the following directory structure relative to the project root:
#
# project/
# ├── data/
# │   └── chelsa_links/      # text files with raster URLs (bio*.txt)
# ├── output/
# │   └── chelsa_processed/  # output CSV files
#
# If running interactively, make sure your working directory is set to the project root.
# ------------------------------------------------------------------------------

txt_dir <- "data/chelsa_links"
out_dir <- "output/chelsa_processed"


# ------------------------------------------------------------------------------
# User settings
# ------------------------------------------------------------------------------
txt_dir <- "~/Chelsa_2026/CHELSA_new/prub"
out_dir <- "~/Chelsa_2026/CHELSA_new/prub/output_test"

dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# Time slices in years before present
time_slices <- seq(0, 20000, by = 1000)
#test - time_slices <- seq(0, 3000, by = 1000)

# Aggregation factor
agg_factor <- c(100, 100)

# Expected CHELSA bioclimatic variables
expected_vars <- sprintf("bio%02d", 1:19)
#test - expected_vars <- sprintf("bio%02d", 1:3)

# ------------------------------------------------------------------------------
# Helper function: process one time slice
# ------------------------------------------------------------------------------
process_timeslice <- function(time_slice,
                              txt_dir,
                              out_dir,
                              agg_factor = c(100, 100),
                              expected_vars = sprintf("bio%02d", 1:19)) {
  
  message("--------------------------------------------------")
  message("Processing time slice: ", time_slice)
  
  txt_file <- file.path(txt_dir, paste0("bio", time_slice, ".txt"))
  
  if (!file.exists(txt_file)) {
    stop("File not found: ", txt_file)
  }
  
  # Read raster URLs
  rastlist <- readLines(txt_file, warn = FALSE)
  rastlist <- trimws(rastlist)
  rastlist <- rastlist[rastlist != ""]
  
  if (length(rastlist) != length(expected_vars)) {
    stop("Expected exactly ", length(expected_vars), " rasters in ", txt_file,
         ", but found ", length(rastlist))
  }
  
  # Extract variable names from filenames
  file_names <- basename(rastlist)
  var_names <- sub(".*_(bio[0-9]{2})_.*", "\\1", file_names)
  
  if (any(var_names == file_names)) {
    stop("Could not extract variable name from: ",
         paste(file_names, collapse = ", "))
  }
  
  if (anyDuplicated(var_names)) {
    stop("Duplicated variable names in ", txt_file, ": ",
         paste(var_names, collapse = ", "))
  }
  
  if (!setequal(var_names, expected_vars)) {
    stop("Unexpected variables in ", txt_file,
         ". Found: ", paste(sort(var_names), collapse = ", "),
         ". Expected: ", paste(expected_vars, collapse = ", "))
  }
  
  # Reorder rasters to ensure consistent layer order across time slices
  ord <- match(expected_vars, var_names)
  rastlist <- rastlist[ord]
  var_names <- var_names[ord]
  
  message("Variables found: ", paste(var_names, collapse = ", "))
  
  # Read rasters from remote URLs using GDAL's /vsicurl/
  rastlist_vsi <- paste0("/vsicurl/", rastlist)
  
  r <- tryCatch(
    terra::rast(rastlist_vsi),
    error = function(e) {
      stop(
        "Error reading raster(s) from:\n",
        paste(rastlist, collapse = "\n"),
        "\n",
        e$message
      )
    }
  )
  
  names(r) <- var_names
  
  # Aggregate raster layers
  r_agg <- terra::aggregate(r, fact = agg_factor, fun = mean, na.rm = TRUE)
  
  # Convert to data frame
  raster_df <- as.data.frame(r_agg, xy = TRUE, na.rm = TRUE) %>%
    mutate(time_bp = time_slice, .before = 1)
  
  # Save output for this time slice
  out_file <- file.path(out_dir, paste0("matrix_", time_slice, ".csv"))
  write.csv(raster_df, out_file, row.names = FALSE)
  
  message("Saved: ", out_file)
  message("Rows: ", nrow(raster_df), " | Columns: ", ncol(raster_df))
  
  return(raster_df)
}

# ------------------------------------------------------------------------------
# Process all time slices
# ------------------------------------------------------------------------------
all_results <- map(
  time_slices,
  process_timeslice,
  txt_dir = txt_dir,
  out_dir = out_dir,
  agg_factor = agg_factor,
  expected_vars = expected_vars
)

names(all_results) <- paste0("bio_", time_slices)

# ------------------------------------------------------------------------------
# Combine results across all time slices
# ------------------------------------------------------------------------------
combined_climate <- bind_rows(all_results)

combined_file <- file.path(out_dir, "matrix_all_timeslices.csv")
write.csv(combined_climate, combined_file, row.names = FALSE)

message("Combined matrix saved: ", combined_file)
message("Combined rows: ", nrow(combined_climate),
        " | Columns: ", ncol(combined_climate))

# ------------------------------------------------------------------------------
# PCA of bio01-bio19 using all time slices together
# ------------------------------------------------------------------------------
bio_vars <- expected_vars

# Keep only complete cases for PCA
pca_input <- combined_climate %>%
  select(time_bp, x, y, all_of(bio_vars)) %>%
  na.omit()

if (nrow(pca_input) == 0) {
  stop("PCA input contains no complete cases after removing NA values.")
}

# Run PCA on scaled variables
pca_res <- prcomp(pca_input[, bio_vars], center = TRUE, scale. = TRUE)

# PCA scores for each grid cell and time slice
pca_scores <- as.data.frame(pca_res$x) %>%
  bind_cols(pca_input %>% select(time_bp, x, y), .) %>%
  select(time_bp, x, y, everything())

# Variable loadings
pca_loadings <- as.data.frame(pca_res$rotation)
pca_loadings$variable <- rownames(pca_loadings)
rownames(pca_loadings) <- NULL
pca_loadings <- pca_loadings %>%
  select(variable, everything())

# Explained variance
eigenvalues <- pca_res$sdev^2
pca_variance <- data.frame(
  PC = paste0("PC", seq_along(eigenvalues)),
  eigenvalue = eigenvalues,
  proportion_variance = eigenvalues / sum(eigenvalues),
  cumulative_variance = cumsum(eigenvalues / sum(eigenvalues))
)

# Save PCA outputs
write.csv(
  pca_scores,
  file.path(out_dir, "matrix_all_timeslices_PCA_scores.csv"),
  row.names = FALSE
)

write.csv(
  pca_loadings,
  file.path(out_dir, "matrix_all_timeslices_PCA_loadings.csv"),
  row.names = FALSE
)

write.csv(
  pca_variance,
  file.path(out_dir, "matrix_all_timeslices_PCA_variance.csv"),
  row.names = FALSE
)

message("PCA completed successfully.")
message("Saved PCA scores, loadings, and variance tables.")
message("==================================================")
message("SCRIPT COMPLETED SUCCESSFULLY")
