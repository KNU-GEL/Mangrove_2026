############################################################
# Test labeling Excel and PDF sheets for selected training points
#
# Input:
#   /home/hjlee/project/mangrove/training_data/training_points.gpkg
#
# Outputs:
#   /home/hjlee/project/mangrove/training_data/
#     training_points_labeling_test_1_10_n10.xlsx
#     training_points_label_sheets_test_1_10_n10.pdf
#
# PDF layout:
#   Left  : SWIR annual series in gray + EVI2 monthly series in blue
#   Right : ID, site, ESA WorldCover, Hansen LCLU, DEM
#
# Notes:
#   - EVI2 y-axis is fixed from 0 to 0.8.
#   - SWIR rasters are 10000-scaled reflectance.
############################################################

library(raster)
library(dplyr)
library(lubridate)
library(rjson)

############################################################
# 0. Paths and settings
############################################################

params <- fromJSON(file = "/home/hjlee/project/mangrove/code/PLCM_Parameters.json")

outDir <- params$setup$outDir
numChunks <- params$setup$numChunks

training_dir <- "/home/hjlee/project/mangrove/training_data"

in_points_gpkg <- file.path(training_dir, "training_points.gpkg")

analysis_start <- as.Date("2017-07-01")
analysis_end   <- as.Date("2026-06-30")

############################################################
# 0-1. Test / subset settings
############################################################

# selected_final_ids <- 4990:5090
selected_final_ids <- integer(0)

run_subset <- length(selected_final_ids) > 0

if (run_subset) {
  run_tag <- paste0(
    "test_",
    min(selected_final_ids),
    "_",
    max(selected_final_ids),
    "_n",
    length(selected_final_ids)
  )
} else {
  run_tag <- "all"
}

out_xlsx <- file.path(training_dir, paste0("training_points_labeling_", run_tag, ".xlsx"))
out_pdf  <- file.path(training_dir, paste0("training_points_label_sheets_", run_tag, ".pdf"))

out_evi2_rds <- file.path(training_dir, paste0("training_points_evi2_timeseries_", run_tag, ".rds"))
out_swir_rds <- file.path(training_dir, paste0("training_points_swir_timeseries_", run_tag, ".rds"))
out_monthly_rds <- file.path(training_dir, paste0("training_points_monthly_timeseries_", run_tag, ".rds"))

############################################################
# 0-2. Input products
############################################################

swir_years <- 2018:2025

# SWIR tif가 10000배 스케일된 반사도이므로 10000으로 나눔
swir_scale_factor <- 10000

site_periods <- list(
  # sanur     = c("1721", "2024", "2326"),
  # gilimanuk = c("1722", "2126"),
  # jembrana  = c("1722", "2126"),
  # lembongan = c("1722", "2126")
  # hungyen   = c("1721", "2024", "2326"),
  aceh   = c("1722", "2126")
)

lct_dir <- "/home/hjlee/project/mangrove/lct_product"

class_rasters <- list(
  ESA_WorldCover_2021 = file.path(lct_dir, "ESA_WorldCover_2021.tif"),
  Hansen_LCLU = file.path(lct_dir, "Hansen_LCLU.tif")
)

############################################################
# 1. Helper functions
############################################################

safe_median <- function(x) {
  if (all(is.na(x))) return(NA_real_)
  median(x, na.rm = TRUE)
}


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


esa_name <- function(x) {
  out <- rep(NA_character_, length(x))
  
  out[x == 10]  <- "Tree cover"
  out[x == 20]  <- "Shrubland"
  out[x == 30]  <- "Grassland"
  out[x == 40]  <- "Cropland"
  out[x == 50]  <- "Built-up"
  out[x == 60]  <- "Bare/sparse vegetation"
  out[x == 70]  <- "Snow and ice"
  out[x == 80]  <- "Permanent water"
  out[x == 90]  <- "Herbaceous wetland"
  out[x == 95]  <- "Mangroves"
  out[x == 100] <- "Moss and lichen"
  
  out
}


hansen_annual_name <- function(x) {
  
  out <- rep(NA_character_, length(x))
  
  out[is.na(x)] <- NA_character_
  
  # 0-24: Terra Firma short vegetation / desert / semi-arid / dense short vegetation
  out[x >= 0 & x <= 1] <- "Terra Firma / True desert / 3-7% short vegetation"
  out[x >= 2 & x <= 18] <- "Terra Firma / Semi-arid / 11-75% short vegetation"
  out[x >= 19 & x <= 24] <- "Terra Firma / Dense short vegetation / 79-100% short vegetation"
  
  # 25-48: Terra Firma tree cover
  out[x >= 25 & x <= 48] <- "Terra Firma / Tree cover"
  
  # 49-96: Terra Firma not used
  out[x >= 49 & x <= 96] <- "Terra Firma / Not used"
  
  # 97-99: Not used
  out[x >= 97 & x <= 99] <- "Not used"
  
  # 100-124: Wetland short vegetation
  out[x >= 100 & x <= 101] <- "Wetland / Salt pan / 3-7% short vegetation"
  out[x >= 102 & x <= 118] <- "Wetland / Sparse vegetation / 11-75% short vegetation"
  out[x >= 119 & x <= 124] <- "Wetland / Dense short vegetation / 79-100% short vegetation"
  
  # 125-148: Wetland tree cover
  out[x >= 125 & x <= 148] <- "Wetland / Tree cover"
  
  # 149-196: Wetland not used
  out[x >= 149 & x <= 196] <- "Wetland / Not used"
  
  # 197-199: Not used
  out[x >= 197 & x <= 199] <- "Not used"
  
  # 200-207: Open surface water
  out[x == 200] <- "Open surface water / 20-29% of year"
  out[x == 201] <- "Open surface water / 30-39% of year"
  out[x == 202] <- "Open surface water / 40-49% of year"
  out[x == 203] <- "Open surface water / 50-59% of year"
  out[x == 204] <- "Open surface water / 60-69% of year"
  out[x == 205] <- "Open surface water / 70-79% of year"
  out[x == 206] <- "Open surface water / 80-89% of year"
  out[x == 207] <- "Open surface water / 90-100% of year"
  
  # 208-240: Not used
  out[x >= 208 & x <= 240] <- "Not used"
  
  # Other classes
  out[x == 241] <- "Snow/ice"
  out[x >= 242 & x <= 243] <- "Not used"
  out[x == 244] <- "Cropland"
  out[x >= 245 & x <= 249] <- "Not used"
  out[x == 250] <- "Built-up"
  out[x >= 251 & x <= 253] <- "Not used"
  out[x == 254] <- "Ocean"
  out[x == 255] <- "No data"
  
  out
}


extract_values_from_raster <- function(r, pts_sf) {
  pts_sp <- as(pts_sf, "Spatial")
  raster::extract(r, pts_sp)
}


ensure_point_index_columns <- function(pts) {
  
  for (site in unique(pts$site)) {
    
    idx <- which(pts$site == site)
    if (length(idx) == 0) next
    
    base_img_path <- file.path(outDir, site, "base_image.tif")
    
    if (!file.exists(base_img_path)) {
      stop("Missing base image: ", base_img_path)
    }
    
    r <- raster(base_img_path)
    
    need_cell <- !("cell" %in% names(pts)) || any(is.na(pts$cell[idx]))
    
    need_chunk <- !all(c("chunk_id", "row_in_chunk") %in% names(pts)) ||
      any(is.na(pts$chunk_id[idx])) ||
      any(is.na(pts$row_in_chunk[idx]))
    
    if (need_cell || need_chunk) {
      
      pts_site <- st_transform(pts[idx, ], st_crs(projection(r)))
      xy <- st_coordinates(pts_site)
      
      cell <- cellFromXY(r, xy)
      
      ck <- cell_to_chunk(
        cell = cell,
        ncell_total = ncell(r),
        numChunks = numChunks
      )
      
      pts$cell[idx] <- cell
      pts$chunk_id[idx] <- ck$chunk_id
      pts$row_in_chunk[idx] <- ck$row_in_chunk
    }
  }
  
  pts
}


extract_raster_values <- function(pts, raster_path, out_col) {
  
  if (is.na(raster_path) || !file.exists(raster_path)) {
    warning("Raster not found: ", raster_path)
    pts[[out_col]] <- NA
    return(pts)
  }
  
  r <- raster(raster_path)
  pts_r <- st_transform(pts, st_crs(projection(r)))
  
  vals <- extract_values_from_raster(r, pts_r)
  
  pts[[out_col]] <- vals
  
  pts
}


extract_site_dem_values <- function(pts) {
  
  pts$DEM <- NA_real_
  
  for (site in unique(pts$site)) {
    
    idx <- which(pts$site == site)
    if (length(idx) == 0) next
    
    dem_dir <- file.path(outDir, site, "dem")
    
    if (!dir.exists(dem_dir)) {
      warning("DEM directory not found: ", dem_dir)
      next
    }
    
    dem_files <- list.files(
      dem_dir,
      pattern = "\\.tif$",
      full.names = TRUE
    )
    
    if (length(dem_files) == 0) {
      warning("No DEM tif found in: ", dem_dir)
      next
    }
    
    if (length(dem_files) > 1) {
      warning(
        "Multiple DEM tifs found in ",
        dem_dir,
        ". Using first one: ",
        dem_files[1]
      )
    }
    
    dem_r <- raster(dem_files[1])
    pts_site <- st_transform(pts[idx, ], st_crs(projection(dem_r)))
    
    pts$DEM[idx] <- extract_values_from_raster(dem_r, pts_site)
    
    cat("Extracted DEM:", site, "|", basename(dem_files[1]), "\n")
  }
  
  pts
}


extract_evi2_timeseries_from_chunks <- function(pts) {
  
  pts_df <- sf::st_drop_geometry(pts)
  
  required_cols <- c(
    "final_id",
    "sample_id",
    "site",
    "chunk_id",
    "row_in_chunk"
  )
  
  if (!all(required_cols %in% names(pts_df))) {
    stop(
      "Required columns are missing: ",
      paste(setdiff(required_cols, names(pts_df)), collapse = ", ")
    )
  }
  
  out_list <- list()
  counter <- 1
  
  for (site in unique(pts_df$site)) {
    
    cat("\nExtracting EVI2 for site:", site, "\n")
    
    if (!(site %in% names(site_periods))) {
      warning("Unknown site in points: ", site)
      next
    }
    
    site_df <- pts_df[pts_df$site == site, ]
    periods <- site_periods[[site]]
    site_dir <- file.path(outDir, site)
    
    for (period in periods) {
      
      chunk_dir <- file.path(site_dir, paste0("chunk_", period))
      
      if (!dir.exists(chunk_dir)) {
        warning("Missing chunk directory: ", chunk_dir)
        next
      }
      
      cat("  Period:", period, "\n")
      
      for (ck in sort(unique(site_df$chunk_id))) {
        
        idx <- which(site_df$chunk_id == ck)
        if (length(idx) == 0) next
        
        ck_num <- sprintf("%03d", ck)
        chunk_file <- file.path(chunk_dir, paste0("chunk_", ck_num, ".rda"))
        
        if (!file.exists(chunk_file)) {
          warning("Missing chunk file: ", chunk_file)
          next
        }
        
        e <- new.env()
        load(chunk_file, envir = e)
        
        rows <- site_df$row_in_chunk[idx]
        
        valid_row <- rows >= 1 & rows <= nrow(e$band1)
        
        if (!all(valid_row)) {
          warning("Invalid row_in_chunk found in: ", chunk_file)
          idx <- idx[valid_row]
          rows <- rows[valid_row]
        }
        
        if (length(rows) == 0) {
          rm(e)
          gc()
          next
        }
        
        dates <- as.Date(e$dates)
        n_dates <- length(dates)
        n_pts <- length(rows)
        
        red <- e$band3[rows, , drop = FALSE] / 10000
        nir <- e$band4[rows, , drop = FALSE] / 10000
        
        evi2 <- 2.5 * (nir - red) / (nir + 2.4 * red + 1)
        
        tmp <- data.frame(
          final_id = rep(site_df$final_id[idx], each = n_dates),
          id = rep(site_df$id[idx], each = n_dates),
          sample_id = rep(site_df$sample_id[idx], each = n_dates),
          site = rep(site, length(idx) * n_dates),
          period = period,
          chunk_id = rep(site_df$chunk_id[idx], each = n_dates),
          row_in_chunk = rep(site_df$row_in_chunk[idx], each = n_dates),
          date = rep(dates, times = n_pts),
          evi2 = as.vector(t(evi2))
        )
        
        out_list[[counter]] <- tmp
        counter <- counter + 1
        
        rm(e, red, nir, evi2, tmp)
        gc()
      }
    }
  }
  
  dplyr::bind_rows(out_list)
}


extract_swir_timeseries_from_annual_rasters <- function(pts) {
  
  out_list <- list()
  counter <- 1
  
  for (site in unique(pts$site)) {
    
    idx <- which(pts$site == site)
    if (length(idx) == 0) next
    
    swir_dir <- file.path(outDir, site, "swir")
    
    if (!dir.exists(swir_dir)) {
      warning("SWIR directory not found: ", swir_dir)
      next
    }
    
    cat("\nExtracting SWIR for site:", site, "\n")
    
    for (yy in swir_years) {
      
      swir_file <- file.path(swir_dir, paste0("swir1_", yy, ".tif"))
      
      if (!file.exists(swir_file)) {
        warning("Missing SWIR file: ", swir_file)
        next
      }
      
      r <- raster(swir_file)
      pts_site <- st_transform(pts[idx, ], st_crs(projection(r)))
      
      vals <- extract_values_from_raster(r, pts_site)
      vals <- vals / swir_scale_factor
      
      tmp <- data.frame(
        final_id = pts$final_id[idx],
        id = pts$id[idx],
        sample_id = pts$sample_id[idx],
        site = site,
        swir_year = yy,
        date = as.Date(paste0(yy, "-07-01")),
        swir = vals
      )
      
      out_list[[counter]] <- tmp
      counter <- counter + 1
      
      cat("  SWIR:", yy, "|", basename(swir_file), "\n")
    }
  }
  
  dplyr::bind_rows(out_list)
}


make_monthly_timeseries <- function(evi2_ts, swir_ts) {
  
  evi2_monthly <- evi2_ts %>%
    dplyr::filter(date >= analysis_start, date <= analysis_end) %>%
    dplyr::mutate(
      evi2 = ifelse(is.finite(evi2), evi2, NA_real_),
      month = floor_date(date, "month")
    ) %>%
    dplyr::group_by(final_id, id, sample_id, site, month) %>%
    dplyr::summarise(
      evi2 = safe_median(evi2),
      .groups = "drop"
    )
  
  if (!is.null(swir_ts) && nrow(swir_ts) > 0) {
    
    swir_monthly <- swir_ts %>%
      dplyr::mutate(month = floor_date(date, "month")) %>%
      dplyr::group_by(final_id, month) %>%
      dplyr::summarise(
        swir = safe_median(swir),
        .groups = "drop"
      )
    
    out <- evi2_monthly %>%
      dplyr::left_join(swir_monthly, by = c("final_id", "month"))
    
  } else {
    out <- evi2_monthly
    out$swir <- NA_real_
  }
  
  out %>%
    dplyr::arrange(site, id, month)
}


make_labeling_excel <- function(pts, out_xlsx) {
  
  if (!requireNamespace("openxlsx", quietly = TRUE)) {
    stop("Package 'openxlsx' is required. Install it with install.packages('openxlsx')")
  }
  
  df <- sf::st_drop_geometry(pts)
  
  if (!("id" %in% names(df))) {
    df$id <- df$final_id
  }
  
  wanted_cols <- c(
    "id",
    "site",
    "label",
    "confidence",
    "memo"
  )
  
  for (cc in wanted_cols) {
    if (!(cc %in% names(df))) {
      df[[cc]] <- NA
    }
  }
  
  df_out <- df[, wanted_cols]
  
  wb <- openxlsx::createWorkbook()
  openxlsx::addWorksheet(wb, "labeling")
  
  openxlsx::writeData(wb, "labeling", df_out)
  
  header_style <- openxlsx::createStyle(
    textDecoration = "bold",
    fgFill = "#D9EAF7",
    border = "Bottom"
  )
  
  openxlsx::addStyle(
    wb,
    "labeling",
    header_style,
    rows = 1,
    cols = seq_len(ncol(df_out)),
    gridExpand = TRUE
  )
  
  openxlsx::freezePane(wb, "labeling", firstRow = TRUE)
  openxlsx::setColWidths(wb, "labeling", cols = 1:ncol(df_out), widths = "auto")
  
  label_cols <- which(names(df_out) %in% c("label", "confidence", "memo"))
  
  if (length(label_cols) > 0) {
    openxlsx::setColWidths(wb, "labeling", cols = label_cols, widths = 22)
  }
  
  openxlsx::saveWorkbook(wb, out_xlsx, overwrite = TRUE)
}


plot_label_page <- function(point_row, ts_monthly) {
  
  fid <- point_row$final_id[1]
  pid <- point_row$id[1]
  
  d <- ts_monthly %>%
    dplyr::filter(final_id == fid) %>%
    dplyr::arrange(month)
  
  site_txt <- as.character(point_row$site[1])
  
  esa_val <- if ("ESA_WorldCover_2021" %in% names(point_row)) {
    point_row$ESA_WorldCover_2021[1]
  } else {
    NA
  }
  
  esa_txt <- if ("ESA_WorldCover_name" %in% names(point_row)) {
    as.character(point_row$ESA_WorldCover_name[1])
  } else {
    "NA"
  }
  
  hansen_val <- if ("Hansen_LCLU" %in% names(point_row)) {
    point_row$Hansen_LCLU[1]
  } else {
    NA
  }
  
  hansen_txt <- if ("Hansen_LCLU_name" %in% names(point_row)) {
    as.character(point_row$Hansen_LCLU_name[1])
  } else {
    "NA"
  }
  
  dem_txt <- if ("DEM" %in% names(point_row)) {
    ifelse(
      is.na(point_row$DEM[1]),
      "NA",
      as.character(round(as.numeric(point_row$DEM[1]), 2))
    )
  } else {
    "NA"
  }
  
  old_par <- par(no.readonly = TRUE)
  on.exit(par(old_par))
  
  layout(matrix(c(1, 2), nrow = 1), widths = c(0.78, 0.22))
  
  ##########################################################
  # Left: time series
  ##########################################################
  
  par(mar = c(4.6, 4.8, 2.2, 1.0))
  
  plot(
    NA,
    xlim = c(analysis_start, analysis_end),
    ylim = c(0, 0.8),
    xlab = "Date",
    ylab = "Reflectance / Index value",
    main = paste0("Point ", pid, " | ", site_txt),
    type = "n"
  )
  
  abline(
    v = as.Date(paste0(2018:2026, "-01-01")),
    col = "gray85",
    lty = 3
  )
  
  abline(
    h = seq(0, 0.8, by = 0.2),
    col = "gray90",
    lty = 3
  )
  
  if (nrow(d) > 0) {
    
    # SWIR first, behind EVI2
    if ("swir" %in% names(d) && any(!is.na(d$swir))) {
      
      d_swir <- d[!is.na(d$swir), ]
      
      lines(
        d_swir$month,
        d_swir$swir,
        lwd = 1.6,
        col = "gray60"
      )
      
      points(
        d_swir$month,
        d_swir$swir,
        pch = 16,
        cex = 0.55,
        col = "gray60"
      )
    }
    
    # EVI2 after SWIR, in front
    points(
      d$month,
      d$evi2,
      pch = 16,
      cex = 0.45,
      col = "royalblue3"
    )
    
    lines(
      d$month,
      d$evi2,
      lwd = 2.4,
      col = "royalblue3"
    )
  }
  
  legend(
    "topright",
    legend = c("EVI2 monthly median", "SWIR annual"),
    col = c("royalblue3", "gray60"),
    lty = c(1, 1),
    lwd = c(2.4, 1.6),
    pch = c(16, 16),
    bty = "n",
    cex = 0.85
  )
  
  box(lwd = 1.5)
  
  ##########################################################
  # Right: text panel
  ##########################################################
  
  par(mar = c(2, 1.2, 2, 1.2))
  plot.new()
  plot.window(xlim = c(0, 1), ylim = c(0, 1))
  
  text(
    x = 0,
    y = 0.96,
    labels = paste0("ID: ", pid),
    adj = c(0, 1),
    cex = 1.30,
    font = 2
  )
  
  text(
    x = 0,
    y = 0.84,
    labels = paste0("Site: ", site_txt),
    adj = c(0, 1),
    cex = 1.08
  )
  
  text(
    x = 0,
    y = 0.68,
    labels = "ESA WorldCover:",
    adj = c(0, 1),
    cex = 0.96,
    font = 2
  )
  
  text(
    x = 0,
    y = 0.60,
    labels = paste0(esa_val, " / ", esa_txt),
    adj = c(0, 1),
    cex = 0.82
  )
  
  text(
    x = 0,
    y = 0.45,
    labels = "Hansen LCLU:",
    adj = c(0, 1),
    cex = 0.96,
    font = 2
  )
  
  hansen_line <- paste0(hansen_val, " / ", hansen_txt)
  hansen_wrapped <- strwrap(hansen_line, width = 28)
  
  y0 <- 0.37
  for (ii in seq_along(hansen_wrapped)) {
    text(
      x = 0,
      y = y0 - (ii - 1) * 0.055,
      labels = hansen_wrapped[ii],
      adj = c(0, 1),
      cex = 0.78
    )
  }
  
  text(
    x = 0,
    y = 0.15,
    labels = paste0("DEM: ", dem_txt, " m"),
    adj = c(0, 1),
    cex = 0.96,
    font = 2
  )
}


make_labeling_pdf <- function(pts, ts_monthly, out_pdf) {
  
  pdf(out_pdf, width = 11.5, height = 4.8)
  
  for (i in seq_len(nrow(pts))) {
    
    plot_label_page(pts[i, ], ts_monthly)
    
    if (i %% 100 == 0 || i == nrow(pts)) {
      cat("PDF pages written:", i, "of", nrow(pts), "\n")
    }
  }
  
  dev.off()
}

############################################################
# 2. Read training points
############################################################

if (!file.exists(in_points_gpkg)) {
  stop("Input training points not found: ", in_points_gpkg)
}

layers <- st_layers(in_points_gpkg)

pts <- st_read(
  in_points_gpkg,
  layer = layers$name[1],
  quiet = TRUE
)

if (!("final_id" %in% names(pts))) {
  pts$final_id <- seq_len(nrow(pts))
}

if (!("sample_id" %in% names(pts))) {
  pts$sample_id <- pts$final_id
}

if (!("site" %in% names(pts))) {
  stop("site column is missing in training points.")
}

# 라벨링용 대표 ID
pts$id <- pts$final_id

############################################################
# 2-1. Subset points for test run
############################################################

if (run_subset) {
  
  cat("Running subset mode.\n")
  cat("Requested final_id values:\n")
  print(selected_final_ids)
  
  pts <- pts[pts$final_id %in% selected_final_ids, ]
  
  if (nrow(pts) == 0) {
    stop("No points matched selected_final_ids.")
  }
  
  cat("Number of selected points:", nrow(pts), "\n")
  cat("Selected final_id values actually found:\n")
  print(sort(pts$final_id))
}

############################################################
# 3. Ensure cell/chunk/row columns
############################################################

pts <- ensure_point_index_columns(pts)

############################################################
# 4. Extract class products and DEM
############################################################

for (nm in names(class_rasters)) {
  pts <- extract_raster_values(
    pts = pts,
    raster_path = class_rasters[[nm]],
    out_col = nm
  )
}

pts$ESA_WorldCover_name <- esa_name(pts$ESA_WorldCover_2021)
pts$Hansen_LCLU_name <- hansen_annual_name(pts$Hansen_LCLU)

pts <- extract_site_dem_values(pts)

############################################################
# 5. Add empty labeling columns
############################################################

pts$label <- NA_character_
pts$confidence <- NA_character_
pts$memo <- NA_character_

############################################################
# 6. Extract EVI2 time series
############################################################

if (file.exists(out_evi2_rds)) {
  cat("Loading cached EVI2 time series:\n")
  cat(out_evi2_rds, "\n")
  evi2_ts <- readRDS(out_evi2_rds)
} else {
  evi2_ts <- extract_evi2_timeseries_from_chunks(pts)
  
  evi2_ts <- evi2_ts %>%
    dplyr::filter(date >= analysis_start, date <= analysis_end)
  
  saveRDS(evi2_ts, out_evi2_rds)
}

############################################################
# 7. Extract SWIR time series
############################################################

if (file.exists(out_swir_rds)) {
  cat("Loading cached SWIR time series:\n")
  cat(out_swir_rds, "\n")
  swir_ts <- readRDS(out_swir_rds)
} else {
  swir_ts <- extract_swir_timeseries_from_annual_rasters(pts)
  saveRDS(swir_ts, out_swir_rds)
}

############################################################
# 8. Create monthly time series
############################################################

ts_monthly <- make_monthly_timeseries(
  evi2_ts = evi2_ts,
  swir_ts = swir_ts
)

saveRDS(ts_monthly, out_monthly_rds)

############################################################
# 9. Save Excel
############################################################

make_labeling_excel(
  pts = pts,
  out_xlsx = out_xlsx
)

cat("Saved Excel:\n")
cat(out_xlsx, "\n")

############################################################
# 10. Save PDF
############################################################

make_labeling_pdf(
  pts = pts,
  ts_monthly = ts_monthly,
  out_pdf = out_pdf
)

cat("Saved PDF:\n")
cat(out_pdf, "\n")

cat("\nDone.\n")