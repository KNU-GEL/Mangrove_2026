############################################################
# Adjusted area estimation + confusion matrices
#
# Input:
#   1. Labeled Excel files:
#      /home/hjlee/project/mangrove/Img_cliped/{site}/reference/YYYY/reference_labeling_YYYY.xlsx
#
#      sheet: labeling
#      columns:
#        reference_id
#        year
#        predicted_class
#        reference_class
#
#   2. Final HMM maps:
#      /home/hjlee/project/mangrove/Img_cliped/{site}/output_hmm/02_hmm_map_YYYY.tif
#
# Output:
#   /home/hjlee/project/mangrove/Img_cliped/{site}/reference/evaluation/
#     {site}_adjusted_area_2018_2025.xlsx
#     {site}_confusion_matrices_2018_2025.xlsx
#
# Notes:
#   - GPKG is not used.
#   - Labeled Excel is the reference-label input.
#   - Confusion matrix:
#       row    = predicted_class
#       column = reference_class
#   - Adjusted area:
#       p_hat_.k = sum_i W_i * n_ik / n_i.
############################################################

library(terra)
library(dplyr)
library(openxlsx)

############################################################
# 0. Settings
############################################################

site <- "hungyen"
years <- 2018:2025

rootDir <- "/home/hjlee/project/mangrove/Img_cliped"
BaseDir <- file.path(rootDir, site)

hmm_dir <- file.path(BaseDir, "output_hmm")
reference_root <- file.path(BaseDir, "reference")
evaluation_dir <- file.path(reference_root, "evaluation")

if (!dir.exists(evaluation_dir)) {
  dir.create(evaluation_dir, recursive = TRUE)
}

class_table <- data.frame(
  class_code = 1:6,
  class_name = c(
    "mangrove",
    "terrestrial",
    "crop",
    "built_up",
    "deep_water",
    "shallow_water"
  ),
  stringsAsFactors = FALSE
)

class_codes <- class_table$class_code
class_labels <- paste0(class_table$class_code, "_", class_table$class_name)

############################################################
# 1. Helper functions
############################################################

normalize_class <- function(x) {
  
  x_chr <- trimws(as.character(x))
  x_chr[x_chr == ""] <- NA_character_
  
  x_num <- suppressWarnings(as.integer(x_chr))
  
  out <- x_num
  out[!(out %in% class_codes)] <- NA_integer_
  
  unresolved <- is.na(out) & !is.na(x_chr)
  
  if (any(unresolved)) {
    
    x_low <- tolower(x_chr[unresolved])
    x_low <- gsub("-", "_", x_low)
    x_low <- gsub(" ", "_", x_low)
    
    name_map <- setNames(class_table$class_code, class_table$class_name)
    matched <- name_map[x_low]
    
    out[unresolved] <- as.integer(matched)
  }
  
  out
}


read_labeling_excel <- function(yy) {
  
  fp <- file.path(
    reference_root,
    as.character(yy),
    sprintf("reference_labeling_%d.xlsx", yy)
  )
  
  if (!file.exists(fp)) {
    stop("Labeling Excel not found: ", fp)
  }
  
  cat("Reading labeling Excel:", fp, "\n")
  
  df <- openxlsx::read.xlsx(
    xlsxFile = fp,
    sheet = "labeling"
  )
  
  required_cols <- c(
    "reference_id",
    "year",
    "predicted_class",
    "reference_class"
  )
  
  if (!all(required_cols %in% names(df))) {
    stop(
      "Missing required columns in ",
      fp,
      ": ",
      paste(setdiff(required_cols, names(df)), collapse = ", ")
    )
  }
  
  df <- df %>%
    mutate(
      reference_id = as.integer(reference_id),
      year = as.integer(year),
      predicted_class = normalize_class(predicted_class),
      reference_class = normalize_class(reference_class)
    )
  
  df$year[is.na(df$year)] <- yy
  
  df
}


get_hmm_area_weights <- function(yy) {
  
  hmm_path <- file.path(
    hmm_dir,
    sprintf("02_hmm_map_%d.tif", yy)
  )
  
  if (!file.exists(hmm_path)) {
    stop("HMM map not found: ", hmm_path)
  }
  
  cat("Reading HMM map:", hmm_path, "\n")
  
  r <- terra::rast(hmm_path)
  
  vals <- terra::values(r, mat = FALSE)
  vals <- as.integer(vals)
  vals[!(vals %in% class_codes)] <- NA_integer_
  
  if (terra::is.lonlat(r)) {
    
    warning("Raster is lon/lat. Using terra::cellSize() for area calculation.")
    
    area_r <- terra::cellSize(r, unit = "m")
    area_vals <- terra::values(area_r, mat = FALSE)
    
    pixel_count <- sapply(class_codes, function(cls) {
      sum(vals == cls, na.rm = TRUE)
    })
    
    area_m2 <- sapply(class_codes, function(cls) {
      sum(area_vals[vals == cls], na.rm = TRUE)
    })
    
  } else {
    
    pixel_area_m2 <- abs(prod(terra::res(r)))
    
    pixel_count <- sapply(class_codes, function(cls) {
      sum(vals == cls, na.rm = TRUE)
    })
    
    area_m2 <- pixel_count * pixel_area_m2
  }
  
  total_area_m2 <- sum(area_m2, na.rm = TRUE)
  
  if (!is.finite(total_area_m2) || total_area_m2 <= 0) {
    stop("Total valid map area is invalid for year: ", yy)
  }
  
  data.frame(
    year = yy,
    predicted_class = class_codes,
    predicted_class_name = class_table$class_name,
    pixel_count = as.numeric(pixel_count),
    area_ha = as.numeric(area_m2) / 10000,
    W = as.numeric(area_m2) / total_area_m2,
    total_area_ha = total_area_m2 / 10000,
    source_map = hmm_path,
    stringsAsFactors = FALSE
  )
}


make_count_matrix <- function(df_valid) {
  
  mat <- table(
    factor(df_valid$predicted_class, levels = class_codes),
    factor(df_valid$reference_class, levels = class_codes)
  )
  
  mat <- as.matrix(mat)
  storage.mode(mat) <- "integer"
  
  rownames(mat) <- class_labels
  colnames(mat) <- class_labels
  
  mat
}


make_area_adjusted_matrix <- function(count_mat, area_weights_df) {
  
  W <- area_weights_df$W[match(class_codes, area_weights_df$predicted_class)]
  row_n <- rowSums(count_mat)
  
  if (any(row_n == 0)) {
    stop(
      "Zero valid reference samples in predicted class: ",
      paste(class_labels[row_n == 0], collapse = ", ")
    )
  }
  
  sample_prop <- count_mat
  
  for (i in seq_along(class_codes)) {
    sample_prop[i, ] <- count_mat[i, ] / row_n[i]
  }
  
  area_prop <- sample_prop
  
  for (i in seq_along(class_codes)) {
    area_prop[i, ] <- W[i] * sample_prop[i, ]
  }
  
  rownames(sample_prop) <- class_labels
  colnames(sample_prop) <- class_labels
  
  rownames(area_prop) <- class_labels
  colnames(area_prop) <- class_labels
  
  list(
    sample_prop = sample_prop,
    area_prop = area_prop,
    row_n = row_n,
    W = W
  )
}


calculate_adjusted_area <- function(yy, count_mat, area_weights_df) {
  
  adj_obj <- make_area_adjusted_matrix(
    count_mat = count_mat,
    area_weights_df = area_weights_df
  )
  
  area_prop <- adj_obj$area_prop
  sample_prop <- adj_obj$sample_prop
  row_n <- adj_obj$row_n
  W <- adj_obj$W
  
  total_area_ha <- unique(area_weights_df$total_area_ha)
  
  if (length(total_area_ha) != 1) {
    stop("Invalid total_area_ha for year: ", yy)
  }
  
  adjusted_prop <- colSums(area_prop, na.rm = TRUE)
  adjusted_area_ha <- total_area_ha * adjusted_prop
  
  se_prop <- rep(NA_real_, length(class_codes))
  
  for (k in seq_along(class_codes)) {
    
    var_terms <- rep(NA_real_, length(class_codes))
    
    for (i in seq_along(class_codes)) {
      
      p_ik <- sample_prop[i, k]
      n_i <- row_n[i]
      
      if (n_i > 1) {
        var_terms[i] <- W[i]^2 * p_ik * (1 - p_ik) / (n_i - 1)
      } else {
        var_terms[i] <- NA_real_
      }
    }
    
    se_prop[k] <- sqrt(sum(var_terms, na.rm = TRUE))
  }
  
  SE_ha <- total_area_ha * se_prop
  CI95_margin_ha <- 1.96 * SE_ha
  
  mapped_area_ha <- area_weights_df$area_ha[
    match(class_codes, area_weights_df$predicted_class)
  ]
  
  data.frame(
    year = yy,
    reference_class = class_codes,
    reference_class_name = class_table$class_name,
    mapped_area_ha = mapped_area_ha,
    adjusted_area_ha = adjusted_area_ha,
    SE_ha = SE_ha,
    CI95_margin_ha = CI95_margin_ha,
    stringsAsFactors = FALSE
  )
}


############################################################
# 2. Main calculation
############################################################

all_area_weights <- list()
all_adjusted_area <- list()
count_matrix_list <- list()

for (yy in years) {
  
  cat("\n====================================\n")
  cat("Processing year:", yy, "\n")
  cat("====================================\n")
  
  label_df <- read_labeling_excel(yy)
  
  label_valid <- label_df %>%
    filter(
      predicted_class %in% class_codes,
      reference_class %in% class_codes
    )
  
  if (nrow(label_valid) == 0) {
    stop("No valid reference labels for year: ", yy)
  }
  
  area_weights_df <- get_hmm_area_weights(yy)
  count_mat <- make_count_matrix(label_valid)
  
  adjusted_area_df <- calculate_adjusted_area(
    yy = yy,
    count_mat = count_mat,
    area_weights_df = area_weights_df
  )
  
  all_area_weights[[as.character(yy)]] <- area_weights_df
  all_adjusted_area[[as.character(yy)]] <- adjusted_area_df
  count_matrix_list[[as.character(yy)]] <- count_mat
}

area_weights_all <- bind_rows(all_area_weights)
adjusted_area_all <- bind_rows(all_adjusted_area)

############################################################
# 3. Write adjusted area workbook
############################################################

adjusted_area_xlsx <- file.path(
  evaluation_dir,
  sprintf("%s_adjusted_area_%d_%d.xlsx", site, min(years), max(years))
)

wb_area <- createWorkbook()

addWorksheet(wb_area, "adjusted_area")
writeData(wb_area, "adjusted_area", adjusted_area_all)

addWorksheet(wb_area, "area_weights")
writeData(wb_area, "area_weights", area_weights_all)

header_style <- createStyle(
  textDecoration = "bold",
  fgFill = "#D9EAF7",
  border = "Bottom"
)

for (sh in names(wb_area)) {
  
  df_tmp <- readWorkbook(wb_area, sheet = sh)
  
  addStyle(
    wb_area,
    sh,
    header_style,
    rows = 1,
    cols = 1:ncol(df_tmp),
    gridExpand = TRUE
  )
  
  freezePane(wb_area, sh, firstRow = TRUE)
  setColWidths(wb_area, sh, cols = 1:ncol(df_tmp), widths = "auto")
}

saveWorkbook(wb_area, adjusted_area_xlsx, overwrite = TRUE)

cat("\nSaved adjusted area workbook:\n")
cat(adjusted_area_xlsx, "\n")

############################################################
# 4. Write confusion matrix workbook
############################################################

confusion_xlsx <- file.path(
  evaluation_dir,
  sprintf("%s_confusion_matrices_%d_%d.xlsx", site, min(years), max(years))
)

wb_conf <- createWorkbook()
addWorksheet(wb_conf, "annual_confusion_matrices")

title_style <- createStyle(
  textDecoration = "bold",
  fontSize = 13,
  fgFill = "#E2F0D9"
)

matrix_header_style <- createStyle(
  textDecoration = "bold",
  fgFill = "#D9EAF7",
  border = "Bottom"
)

number_style <- createStyle(numFmt = "0")
percent_style <- createStyle(numFmt = "0.0%")

write_confusion_block <- function(wb, sheet, yy, count_mat, start_row) {
  
  row_total <- rowSums(count_mat)
  col_total <- colSums(count_mat)
  
  UA <- diag(count_mat) / row_total
  PA <- diag(count_mat) / col_total
  OA <- sum(diag(count_mat)) / sum(count_mat)
  
  block <- data.frame(
    predicted_class = class_labels,
    as.data.frame.matrix(count_mat),
    row_total = as.integer(row_total),
    UA = as.numeric(UA),
    check.names = FALSE
  )
  
  column_total_row <- data.frame(
    predicted_class = "column_total",
    as.data.frame.matrix(t(as.matrix(col_total))),
    row_total = sum(count_mat),
    UA = NA_real_,
    check.names = FALSE
  )
  
  colnames(column_total_row) <- colnames(block)
  
  PA_row <- data.frame(
    predicted_class = "PA",
    as.data.frame.matrix(t(as.matrix(PA))),
    row_total = NA_real_,
    UA = OA,
    check.names = FALSE
  )
  
  colnames(PA_row) <- colnames(block)
  
  out_block <- bind_rows(
    block,
    column_total_row,
    PA_row
  )
  
  writeData(
    wb,
    sheet,
    x = paste0("Year ", yy, " confusion matrix"),
    startRow = start_row,
    startCol = 1
  )
  
  addStyle(
    wb,
    sheet,
    title_style,
    rows = start_row,
    cols = 1,
    gridExpand = TRUE
  )
  
  writeData(
    wb,
    sheet,
    x = out_block,
    startRow = start_row + 1,
    startCol = 1
  )
  
  addStyle(
    wb,
    sheet,
    matrix_header_style,
    rows = start_row + 1,
    cols = 1:ncol(out_block),
    gridExpand = TRUE
  )
  
  count_rows <- (start_row + 2):(start_row + 1 + nrow(out_block))
  count_cols <- 2:(1 + length(class_codes) + 1)
  
  addStyle(
    wb,
    sheet,
    number_style,
    rows = count_rows,
    cols = count_cols,
    gridExpand = TRUE
  )
  
  percent_col <- ncol(out_block)
  
  addStyle(
    wb,
    sheet,
    percent_style,
    rows = count_rows,
    cols = percent_col,
    gridExpand = TRUE
  )
  
  PA_row_excel <- start_row + 1 + nrow(out_block)
  
  addStyle(
    wb,
    sheet,
    percent_style,
    rows = PA_row_excel,
    cols = 2:(1 + length(class_codes)),
    gridExpand = TRUE
  )
  
  start_row + nrow(out_block) + 4
}

current_row <- 1

for (yy in years) {
  
  current_row <- write_confusion_block(
    wb = wb_conf,
    sheet = "annual_confusion_matrices",
    yy = yy,
    count_mat = count_matrix_list[[as.character(yy)]],
    start_row = current_row
  )
}

setColWidths(
  wb_conf,
  "annual_confusion_matrices",
  cols = 1:12,
  widths = "auto"
)

saveWorkbook(wb_conf, confusion_xlsx, overwrite = TRUE)

cat("\nSaved confusion matrix workbook:\n")
cat(confusion_xlsx, "\n")

cat("\nDone.\n")