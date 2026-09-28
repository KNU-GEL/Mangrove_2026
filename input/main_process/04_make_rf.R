# =========================
# chunk_feat_1721 : 2018, 2019, 2020
# chunk_feat_2024 : 2021, 2022, 2023
# chunk_feat_2326 : 2024, 2025
# =========================

library(terra)
library(sf)
library(randomForest)
library(rjson)

## ========== Environment and Path Settings ==========
params <- fromJSON(file = "/home/hjlee/project/mangrove/code/PLCM_Parameters.json")
source(params$setup$rFunctions)

sites    <- c("taklong")
rootDir  <- "/home/hjlee/project/mangrove/Img_cliped"
shpRoot  <- "/home/hjlee/project/mangrove/shp"

## Model Output Path
outRoot  <- file.path(rootDir, "taklong", "model")
if (!dir.exists(outRoot)) dir.create(outRoot, recursive = TRUE)

## SWIR / DEM Paths
swirRoot <- "/home/hjlee/project/mangrove/Img_cliped/taklong/swir"
demRoot  <- "/home/hjlee/project/mangrove/Img_cliped/taklong/dem"

dem_files <- list.files(demRoot, pattern = "\\.tif$", full.names = FALSE)

year_to_train <- 2018 
numFeatUse <- 40

set.seed(456789)

## ========== Utilities ==========

## Select the chunk_feat directory and calculate yToDo for year_to_train
# get_chunk_dir_and_yToDo <- function(BaseDir, target_year) {
#   if (target_year >= 2018 && target_year <= 2020) {
#     return(list(
#       ckPheDir = file.path(BaseDir, "chunk_feat_1721"),
#       yToDo = target_year - 2018 + 1L
#     ))
#   }
#   
#   if (target_year >= 2021 && target_year <= 2023) {
#     return(list(
#       ckPheDir = file.path(BaseDir, "chunk_feat_2024"),
#       yToDo = target_year - 2021 + 1L
#     ))
#   }
#   
#   if (target_year >= 2024 && target_year <= 2025) {
#     return(list(
#       ckPheDir = file.path(BaseDir, "chunk_feat_2326"),
#       yToDo = target_year - 2024 + 1L
#     ))
#   }

  
get_chunk_dir_and_yToDo <- function(BaseDir, target_year) {
  if (target_year >= 2018 && target_year <= 2021) {
    return(list(
      ckPheDir = file.path(BaseDir, "chunk_feat_1722"),
      yToDo = target_year - 2018 + 1L
    ))
  }
  
  if (target_year >= 2022 && target_year <= 2025) {
    return(list(
      ckPheDir = file.path(BaseDir, "chunk_feat_2126"),
      yToDo = target_year - 2021 + 1L
    ))
  }
  
    
  stop("Unsupported year_to_train: ", target_year)
}

prep_vec <- function(v, imgBase, site_name) {
  nms <- names(v)
  hit <- which(tolower(nms) == "site")
  
  if (length(hit) > 0) {
    key <- hit[1]
    v <- v[tolower(trimws(v[[key]])) == tolower(site_name), ]
  }
  
  v <- project(v, crs(imgBase))
  vv <- try(suppressWarnings(crop(v, imgBase)), silent = TRUE)
  if (!inherits(vv, "try-error")) v <- vv
  
  v
}

safe_extract <- function(r_idx, v) {
  if (is.null(v)) return(NULL)
  if (inherits(v, "SpatVector") && nrow(v) == 0) return(NULL)
  extract(r_idx, v)
}

split_80_10_10 <- function(pt_df) {
  if (is.null(pt_df) || nrow(pt_df) == 0) {
    return(list(train = NULL, val = NULL, test = NULL))
  }
  
  list(
    train = pt_df,
    val = NULL,
    test = NULL
  )
}

get_train_idx <- function(sp) {
  if (is.null(sp$train)) integer(0) else as.integer(sp$train[, 2])
}

sorted_chunk_files <- function(dir_path) {
  files <- list.files(dir_path, pattern = "^chunk_feat_\\d+\\.rda$", full.names = TRUE)
  
  if (length(files) == 0) return(character(0))
  
  idx <- suppressWarnings(as.integer(sub(".*chunk_feat_(\\d+)\\.rda$", "\\1", files)))
  files[order(idx, na.last = NA)]
}

build_Fmat_site <- function(ckPheDir, numPix, yToDo, numFeatUse) {
  files <- sorted_chunk_files(ckPheDir)
  
  if (length(files) == 0) {
    stop(sprintf("No chunk_feat files in %s", ckPheDir))
  }
  
  sel_cols <- (58 * (yToDo - 1) + 1):(58 * (yToDo - 1) + numFeatUse)
  
  Fmat <- matrix(NA_real_, nrow = numPix, ncol = numFeatUse)
  colnames(Fmat) <- paste0("F", seq_len(numFeatUse))
  
  pos <- 1L
  
  for (cf in files) {
    log <- try(load(cf), silent = TRUE)  # -> f_mat
    
    if (inherits(log, "try-error") || !exists("f_mat")) next
    
    n <- nrow(f_mat)
    
    if (n <= 0) {
      rm(f_mat)
      next
    }
    
    if (max(sel_cols) > ncol(f_mat)) {
      stop(sprintf(
        "sel_cols exceed ncol(f_mat)=%d in %s",
        ncol(f_mat),
        basename(cf)
      ))
    }
    
    end_pos <- min(pos + n - 1L, numPix)
    block <- f_mat[1:(end_pos - pos + 1L), sel_cols, drop = FALSE]
    
    if (nrow(block) != (end_pos - pos + 1L)) {
      stop(sprintf("row mismatch at %s", basename(cf)))
    }
    
    Fmat[pos:end_pos, ] <- block
    pos <- end_pos + 1L
    
    rm(f_mat)
    
    if ((pos - 1L) %% 1e6 == 0) {
      cat("    filled rows:", pos - 1L, "\n")
    }
  }
  
  if (pos - 1L != numPix) {
    warning(sprintf("Filled %d of %d rows (remaining rows are NA).", pos - 1L, numPix))
  }
  
  Fmat
}

## ========== Prepare DEM Once Globally ==========
dem_paths <- file.path(demRoot, dem_files)
dem_paths <- dem_paths[file.exists(dem_paths)]

if (length(dem_paths) == 0) {
  warning("DEM file not found:", demRoot)
  dem_mosaic <- NULL
  
} else if (length(dem_paths) == 1) {
  dem_mosaic <- rast(dem_paths[1])
  dem_mosaic[dem_mosaic <= -1000] <- NA
  
} else {
  dem_list <- lapply(dem_paths, rast)
  
  dem_mosaic <- do.call(mosaic, c(dem_list, fun = "mean"))
  dem_mosaic[dem_mosaic <= -1000] <- NA
  
  rm(dem_list)
}

## ========== Data Collection and Combination Loop ==========
train_list <- vector("list", 6)
names(train_list) <- paste0("c", 1:6)

errors <- list()

for (site in sites) {
  message("\n===============================\nProcessing site: ", site)
  
  BaseDir <- file.path(rootDir, site)
  
  ## Select the chunk_feat directory for year_to_train
  chunk_info <- get_chunk_dir_and_yToDo(BaseDir, year_to_train)
  ckPheDir <- chunk_info$ckPheDir
  yToDo <- chunk_info$yToDo
  
  message("  chunk feature dir: ", ckPheDir)
  message("  yToDo within selected chunk folder: ", yToDo)
  
  tryCatch({
    
    ## 1) Base Raster / Pixel Index
    base_fp <- file.path(BaseDir, "base_image.tif")
    
    if (!file.exists(base_fp)) {
      stop("base_image not found: ", base_fp)
    }
    
    imgBase <- rast(base_fp)
    numPix  <- ncol(imgBase) * nrow(imgBase)
    imgNum  <- setValues(imgBase, 1:numPix)
    
    ## 2) Convert GPKG Training Data to Pixel Indices
    ## Class Order:
    ## 1 = mangrove
    ## 2 = terrestrial
    ## 3 = crop
    ## 4 = built-up
    ## 5 = deep_water
    ## 6 = shallow_water
    
    t1 <- prep_vec(vect(file.path(shpRoot, "mangrove.gpkg")), imgBase, site)
    t2 <- prep_vec(vect(file.path(shpRoot, "terrestrial.gpkg")), imgBase, site)
    t3 <- prep_vec(vect(file.path(shpRoot, "crop.gpkg")), imgBase, site)
    t4 <- prep_vec(vect(file.path(shpRoot, "built-up.gpkg")), imgBase, site)
    t5 <- prep_vec(vect(file.path(shpRoot, "deep_water.gpkg")), imgBase, site)
    t6 <- prep_vec(vect(file.path(shpRoot, "shallow_water.gpkg")), imgBase, site)
    
    pt1 <- safe_extract(imgNum, t1)
    pt2 <- safe_extract(imgNum, t2)
    pt3 <- safe_extract(imgNum, t3)
    pt4 <- safe_extract(imgNum, t4)
    pt5 <- safe_extract(imgNum, t5)
    pt6 <- safe_extract(imgNum, t6)
    
    count_pt <- function(pt) if (is.null(pt)) 0L else nrow(pt)
    
    message(
      "  pts per class: ",
      paste(sapply(list(pt1, pt2, pt3, pt4, pt5, pt6), count_pt), collapse = ", ")
    )
    
    sp1 <- split_80_10_10(pt1)
    sp2 <- split_80_10_10(pt2)
    sp3 <- split_80_10_10(pt3)
    sp4 <- split_80_10_10(pt4)
    sp5 <- split_80_10_10(pt5)
    sp6 <- split_80_10_10(pt6)
    
    ## 3) Build Fmat (40 features) using yToDo within the selected period folder
    files <- sorted_chunk_files(ckPheDir)
    
    if (length(files) == 0) {
      stop("chunk_feat files not found: ", ckPheDir)
    }
    
    load(files[1])  # -> f_mat
    ncol_f <- ncol(f_mat)
    rm(f_mat)
    
    years_avail <- floor(ncol_f / 58)
    
    if (yToDo > years_avail) {
      stop(sprintf(
        "Insufficient columns for requested year (yToDo=%d; available years=%d; ncol=%d).",
        yToDo,
        years_avail,
        ncol_f
      ))
    }
    
    Fmat <- build_Fmat_site(ckPheDir, numPix, yToDo, numFeatUse)
    
    ## 4) Read SWIR1
    ## File Naming Convention:
    ## /home/hjlee/project/mangrove/Img_cliped/hungyen/swir/swir1_yyyy.tif
    swir_fp <- file.path(swirRoot, sprintf("swir1_%d.tif", year_to_train))
    
    if (!file.exists(swir_fp)) {
      warning("SWIR file not found; replacing with NA: ", swir_fp)
      swir_vec <- rep(NA_real_, numPix)
      
    } else {
      swir <- rast(swir_fp)
      
      if (as.character(crs(swir)) != as.character(crs(imgBase))) {
        swir <- project(swir, imgBase, method = "bilinear")
      }
      
      swir3m <- resample(swir, imgBase, method = "bilinear")
      swir_vec <- as.numeric(values(swir3m))
      
      if (length(swir_vec) != numPix) {
        stop("SWIR length mismatch: ", length(swir_vec), " vs ", numPix)
      }
    }
    
    ## 5) Read DEM and Align to 3 m Resolution
    if (is.null(dem_mosaic)) {
      warning("DEM mosaic not available; replacing with NA.")
      dem_vec <- rep(NA_real_, numPix)
      
    } else {
      dem_use <- dem_mosaic
      
      if (as.character(crs(dem_use)) != as.character(crs(imgBase))) {
        img_in_dem <- project(imgBase, crs(dem_use))
        dem_use <- crop(dem_use, ext(img_in_dem))
        
      } else {
        dem_use <- crop(dem_use, ext(imgBase))
      }
      
      dem3m <- project(dem_use, imgBase, method = "bilinear")
      dem_vec <- as.numeric(values(dem3m))
      
      if (length(dem_vec) != numPix) {
        stop("DEM 길이 불일치: ", length(dem_vec), " vs ", numPix)
      }
    }
    
    ## 6) Combine Feature Stack: F1..F40 + SWIR1 + DEM
    Fmat_ext <- cbind(
      Fmat,
      SWIR1_20m = swir_vec,
      DEM30m = dem_vec
    )
    
    rm(Fmat, swir_vec, dem_vec)
    gc()
    
    ## 7) Accumulate All Training Samples by Class
    idx1 <- get_train_idx(sp1)
    idx2 <- get_train_idx(sp2)
    idx3 <- get_train_idx(sp3)
    idx4 <- get_train_idx(sp4)
    idx5 <- get_train_idx(sp5)
    idx6 <- get_train_idx(sp6)
    
    if (length(idx1) > 0) {
      train_list[[1]] <- rbind(train_list[[1]], Fmat_ext[idx1, , drop = FALSE])
    }
    if (length(idx2) > 0) {
      train_list[[2]] <- rbind(train_list[[2]], Fmat_ext[idx2, , drop = FALSE])
    }
    if (length(idx3) > 0) {
      train_list[[3]] <- rbind(train_list[[3]], Fmat_ext[idx3, , drop = FALSE])
    }
    if (length(idx4) > 0) {
      train_list[[4]] <- rbind(train_list[[4]], Fmat_ext[idx4, , drop = FALSE])
    }
    if (length(idx5) > 0) {
      train_list[[5]] <- rbind(train_list[[5]], Fmat_ext[idx5, , drop = FALSE])
    }
    if (length(idx6) > 0) {
      train_list[[6]] <- rbind(train_list[[6]], Fmat_ext[idx6, , drop = FALSE])
    }
    
    message("  Appended training samples (+SWIR + DEM): ", site)
    
    rm(Fmat_ext)
    gc()
    
  }, error = function(e) {
    message("  [SKIP] ", site, " -> ", conditionMessage(e))
    errors[[site]] <<- conditionMessage(e)
  })
}

if (length(errors)) {
  message("\n--- Summary of Skipped Sites ---")
  print(errors)
}

## ========== Combine Training Inputs ==========
rows_list <- list()
labels <- integer(0)

for (k in 1:6) {
  if (!is.null(train_list[[k]]) && nrow(train_list[[k]]) > 0) {
    rows_list[[length(rows_list) + 1L]] <- train_list[[k]]
    labels <- c(labels, rep(k, nrow(train_list[[k]])))
  }
}

if (!length(rows_list)) {
  stop("No training samples were collected from any site. Check the error summary above.")
}

mat_input <- do.call(rbind, rows_list)

input <- data.frame(
  lct = factor(labels),
  mat_input,
  check.names = FALSE
)

cat("\nFinal training samples:", nrow(input), "\n")
print(table(input$lct))

## ========== Robust Model Training and Saving ==========
y <- factor(input$lct)
X <- input[, -1, drop = FALSE]

Xmat <- as.matrix(X)
Xmat[!is.finite(Xmat)] <- NA
X <- as.data.frame(Xmat)

thr_seq <- c(0.30, 0.50, 0.70, 0.90)
keep_cols <- integer(0)

for (thr in thr_seq) {
  cols <- which(
    colSums(!is.na(X)) > 0 &
      colMeans(is.na(X)) < thr
  )
  
  if (length(cols) > 0) {
    keep_cols <- cols
    cat("Number of features used:", length(cols), "(NA <", thr, ")\n")
    break
  }
}

if (length(keep_cols) == 0L) {
  cat("Warning: All features have high NA proportions. Attempting to use all features after roughfix.\n")
  X2_imp <- randomForest::na.roughfix(X)
  keep_names <- colnames(X)
  
} else {
  X2 <- X[, keep_cols, drop = FALSE]
  X2_imp <- randomForest::na.roughfix(X2)
  keep_names <- colnames(X2)
}

cat("Final training dimensions:", nrow(X2_imp), "x", ncol(X2_imp), "\n")

set.seed(456789)

rf_all <- randomForest(
  X2_imp,
  y,
  importance = TRUE,
  ntree = 500
)

print(rf_all)

imp <- try(importance(rf_all, type = 1), silent = TRUE)

if (!inherits(imp, "try-error")) {
  write.csv(
    imp,
    file = file.path(outRoot, sprintf("rf_%d_varimp.csv", year_to_train)),
    row.names = TRUE
  )
}

preproc <- list(
  keep_cols = keep_names,
  medians = vapply(
    X,
    function(v) if (is.numeric(v)) median(v, na.rm = TRUE) else NA_real_,
    numeric(1)
  )
)

objs <- c("rf_all", "preproc")
missing <- objs[!sapply(objs, exists, inherits = TRUE)]

if (length(missing)) {
  stop("Missing objects: ", paste(missing, collapse = ", "))
}

save(
  list = objs,
  file = file.path(outRoot, sprintf("rf_%d.rda", year_to_train))
)

cat(
  "Saved model + preproc to:",
  file.path(outRoot, sprintf("rf_%d.rda", year_to_train)),
  "\n"
)