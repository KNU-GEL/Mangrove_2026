library(raster)
library(sf)
library(rjson)

############################################################
# Parameters
############################################################

params <- fromJSON(file = "/home/hjlee/project/mangrove/code/PLCM_Parameters.json")

outDir <- params$setup$outDir
numChunks <- params$setup$numChunks

set.seed(20260515)

n_sample <- 1500

# Minimum distance between sampled points in pixel units
# 3 means center-to-center distance should be at least 3 pixels.
min_dist_pixels <- 3

site_periods <- list(
  # sanur     = c("1721", "2024", "2326"),
  # gilimanuk = c("1722", "2126"),
  # jembrana  = c("1722", "2126"),
  # lembongan = c("1722", "2126")
  # hungyen = c("1721", "2024", "2326"),
  # taklong = c("1722", "2126"),
  aceh = c("1722", "2126")
)

out_csv <- file.path(outDir, "training_random_points.csv")
out_gpkg <- file.path(outDir, "training_random_points.gpkg")

############################################################
# Helper functions
############################################################

cell_to_chunk <- function(cell, ncell_total, numChunks) {
  chunk_size <- ncell_total %/% numChunks
  
  chunk_id <- ((cell - 1) %/% chunk_size) + 1
  chunk_id[chunk_id > numChunks] <- numChunks
  
  chunk_start <- chunk_size * (chunk_id - 1) + 1
  row_in_chunk <- cell - chunk_start + 1
  
  data.frame(
    chunk_id = chunk_id,
    row_in_chunk = row_in_chunk
  )
}


allocate_samples_by_area <- function(n_total, weights) {
  raw_n <- n_total * weights / sum(weights)
  base_n <- floor(raw_n)
  remain <- n_total - sum(base_n)
  
  if (remain > 0) {
    add_idx <- order(raw_n - base_n, decreasing = TRUE)[seq_len(remain)]
    base_n[add_idx] <- base_n[add_idx] + 1
  }
  
  base_n
}


make_distance_offsets <- function(min_dist_pixels) {
  max_offset <- ceiling(min_dist_pixels)
  
  offsets <- expand.grid(
    dr = -max_offset:max_offset,
    dc = -max_offset:max_offset
  )
  
  offsets$dist <- sqrt(offsets$dr^2 + offsets$dc^2)
  
  # Reject points that are closer than min_dist_pixels.
  # If you want to reject points exactly 3 pixels away too,
  # change "< min_dist_pixels" to "<= min_dist_pixels".
  offsets <- offsets[offsets$dist <= min_dist_pixels, ]
  
  offsets
}


sample_cells_min_distance <- function(
    r,
    n_target,
    min_dist_pixels = 3,
    candidate_factor = 30,
    max_rounds = 200
) {
  
  nr <- nrow(r)
  nc <- ncol(r)
  nt <- ncell(r)
  
  offsets <- make_distance_offsets(min_dist_pixels)
  
  occupied <- matrix(FALSE, nrow = nr, ncol = nc)
  
  accepted_cells <- integer(0)
  accepted_rows <- integer(0)
  accepted_cols <- integer(0)
  
  batch_size <- min(
    nt,
    max(50000, n_target * candidate_factor)
  )
  
  round_id <- 1
  
  while (length(accepted_cells) < n_target && round_id <= max_rounds) {
    
    candidates <- sample.int(
      nt,
      size = batch_size,
      replace = FALSE
    )
    
    rc <- rowColFromCell(r, candidates)
    
    for (j in seq_along(candidates)) {
      
      if (length(accepted_cells) >= n_target) break
      
      rr <- rc[j, 1]
      cc <- rc[j, 2]
      
      if (occupied[rr, cc]) next
      
      rr_check <- rr + offsets$dr
      cc_check <- cc + offsets$dc
      
      valid <- rr_check >= 1 &
        rr_check <= nr &
        cc_check >= 1 &
        cc_check <= nc
      
      rr_check <- rr_check[valid]
      cc_check <- cc_check[valid]
      
      has_near_point <- any(occupied[cbind(rr_check, cc_check)])
      
      if (!has_near_point) {
        accepted_cells <- c(accepted_cells, candidates[j])
        accepted_rows <- c(accepted_rows, rr)
        accepted_cols <- c(accepted_cols, cc)
        occupied[rr, cc] <- TRUE
      }
    }
    
    cat(
      "  round", round_id,
      "| accepted", length(accepted_cells),
      "of", n_target, "\n"
    )
    
    round_id <- round_id + 1
  }
  
  if (length(accepted_cells) < n_target) {
    warning(
      "Only sampled ", length(accepted_cells),
      " cells out of requested ", n_target,
      ". Try lowering min_dist_pixels or increasing max_rounds."
    )
  }
  
  data.frame(
    cell = accepted_cells,
    row = accepted_rows,
    col = accepted_cols
  )
}

############################################################
# Build site information
############################################################

site_info <- data.frame(
  site = names(site_periods),
  ncell = NA_integer_,
  n_target = NA_integer_,
  stringsAsFactors = FALSE
)

for (i in seq_len(nrow(site_info))) {
  r <- raster(file.path(outDir, site_info$site[i], "base_image.tif"))
  site_info$ncell[i] <- ncell(r)
}

site_info$n_target <- allocate_samples_by_area(
  n_total = n_sample,
  weights = site_info$ncell
)

print(site_info)
cat("Total target samples:", sum(site_info$n_target), "\n")

############################################################
# Sample points site by site with minimum distance
############################################################

samples_list <- list()
sf_list <- list()

global_id <- 1

for (i in seq_len(nrow(site_info))) {
  
  site <- site_info$site[i]
  n_site <- site_info$n_target[i]
  
  cat("\nSampling site:", site, "\n")
  cat("Target points:", n_site, "\n")
  
  r <- raster(file.path(outDir, site, "base_image.tif"))
  
  sampled_cells <- sample_cells_min_distance(
    r = r,
    n_target = n_site,
    min_dist_pixels = min_dist_pixels
  )
  
  if (nrow(sampled_cells) == 0) next
  
  ck <- cell_to_chunk(
    cell = sampled_cells$cell,
    ncell_total = ncell(r),
    numChunks = numChunks
  )
  
  xy <- xyFromCell(r, sampled_cells$cell)
  
  tmp <- data.frame(
    sample_id = global_id:(global_id + nrow(sampled_cells) - 1),
    site = site,
    cell = sampled_cells$cell,
    raster_row = sampled_cells$row,
    raster_col = sampled_cells$col,
    chunk_id = ck$chunk_id,
    row_in_chunk = ck$row_in_chunk,
    x = xy[, 1],
    y = xy[, 2],
    lon = NA_real_,
    lat = NA_real_,
    periods = paste(site_periods[[site]], collapse = ","),
    label = NA_character_,
    memo = NA_character_,
    stringsAsFactors = FALSE
  )
  
  tmp_sf <- st_as_sf(
    tmp,
    coords = c("x", "y"),
    crs = st_crs(projection(r)),
    remove = FALSE
  )
  
  tmp_ll <- st_transform(tmp_sf, 4326)
  ll <- st_coordinates(tmp_ll)
  
  tmp$lon <- ll[, 1]
  tmp$lat <- ll[, 2]
  tmp_ll$lon <- ll[, 1]
  tmp_ll$lat <- ll[, 2]
  
  samples_list[[site]] <- tmp
  sf_list[[site]] <- tmp_ll
  
  global_id <- global_id + nrow(sampled_cells)
  
  cat("Finished site:", site, "| sampled:", nrow(sampled_cells), "\n")
}

############################################################
# Combine and save
############################################################

samples <- do.call(rbind, samples_list)
samples_sf <- do.call(rbind, sf_list)

# Reassign sample_id just in case any site sampled fewer than requested
samples$sample_id <- seq_len(nrow(samples))
samples_sf$sample_id <- seq_len(nrow(samples_sf))

write.csv(samples, out_csv, row.names = FALSE)

if (file.exists(out_gpkg)) file.remove(out_gpkg)

st_write(
  samples_sf,
  out_gpkg,
  layer = "random_points_6000_min3pix",
  delete_layer = TRUE
)

cat("\nSaved CSV :", out_csv, "\n")
cat("Saved GPKG:", out_gpkg, "\n")
cat("Total sampled points:", nrow(samples), "\n")

cat("\nSamples by site:\n")
print(table(samples$site))
