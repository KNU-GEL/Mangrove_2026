library(raster)
library(terra)
library(randomForest)

## ===== 사용자 설정 =====
site <- "hungyen"             # 예측 지역
year_to_predict <- 2025       # Planet 시계열 열 선택 + SWIR 연도
use_nearest_swir <- TRUE      # 해당 연도 SWIR 없으면 가장 가까운 연도 자동 선택

numFeatUse <- 40

## ===== 경로 설정 =====
rootDir   <- "/home/hjlee/project/mangrove/Img_cliped"
BaseDir   <- file.path(rootDir, site)

## 예측 결과 저장 경로
outputDir <- file.path(BaseDir, "output_balimodel")

## 학습된 모델 경로
# modelDir  <- file.path(BaseDir, "model")
modelDir  <- "/home/hjlee/project/mangrove/model_bali"

## SWIR / DEM 경로
swirRoot  <- file.path(BaseDir, "swir")
demRoot   <- file.path(BaseDir, "dem")

if (!dir.exists(outputDir)) dir.create(outputDir, recursive = TRUE)

## ===== 연도별 chunk_feat 폴더 선택 함수 =====
## chunk_feat_1721 : 2018, 2019, 2020
## chunk_feat_2024 : 2021, 2022, 2023
## chunk_feat_2326 : 2024, 2025

get_chunk_dir_and_yToDo <- function(BaseDir, target_year) {
  if (target_year >= 2018 && target_year <= 2020) {
    return(list(
      FeatDir = file.path(BaseDir, "chunk_feat_1721"),
      yToDo = target_year - 2018 + 1L
    ))
  }
  
  if (target_year >= 2021 && target_year <= 2023) {
    return(list(
      FeatDir = file.path(BaseDir, "chunk_feat_2024"),
      yToDo = target_year - 2021 + 1L
    ))
  }
  
  if (target_year >= 2024 && target_year <= 2025) {
    return(list(
      FeatDir = file.path(BaseDir, "chunk_feat_2326"),
      yToDo = target_year - 2024 + 1L
    ))
  }
  
  stop("지원하지 않는 year_to_predict 입니다: ", target_year)
}

chunk_info <- get_chunk_dir_and_yToDo(BaseDir, year_to_predict)
FeatDir <- chunk_info$FeatDir
yToDo   <- chunk_info$yToDo

message("예측 연도: ", year_to_predict)
message("사용 chunk feature dir: ", FeatDir)
message("선택된 chunk 폴더 내부 yToDo: ", yToDo)

## ===== chunk 파일 정렬 함수 =====
sorted_chunk_files <- function(dir_path) {
  files <- list.files(
    dir_path,
    pattern = "^chunk_feat_\\d+\\.rda$",
    full.names = TRUE
  )
  
  if (length(files) == 0) return(character(0))
  
  idx <- suppressWarnings(
    as.integer(sub(".*chunk_feat_(\\d+)\\.rda$", "\\1", files))
  )
  
  files[order(idx, na.last = NA)]
}

## ===== 기본 이미지 =====
base_fp <- file.path(BaseDir, "base_image.tif")

if (!file.exists(base_fp)) {
  stop("base_image 없음: ", base_fp)
}

imgBase <- rast(base_fp)
NumPix  <- ncell(imgBase)

## ===== Planet 시계열 피처(F1..F40) =====
## 기존 numChunks <- 200 방식 대신,
## FeatDir 안의 chunk_feat_*.rda 파일 개수와 row 수를 기준으로 자동 결합

build_Fmat_prediction <- function(FeatDir, NumPix, yToDo, numFeatUse) {
  files <- sorted_chunk_files(FeatDir)
  
  if (length(files) == 0) {
    stop("chunk_feat 파일이 없습니다: ", FeatDir)
  }
  
  sel_cols <- (58L * (yToDo - 1L) + 1L):(58L * (yToDo - 1L) + numFeatUse)
  
  Fmat <- matrix(NA_real_, nrow = NumPix, ncol = numFeatUse)
  colnames(Fmat) <- paste0("F", seq_len(numFeatUse))
  
  pos <- 1L
  
  for (cf in files) {
    log <- try(load(cf), silent = TRUE)  # 기대 객체: f_mat
    
    if (inherits(log, "try-error") || !exists("f_mat")) {
      warning("f_mat 로드 실패, 건너뜀: ", cf)
      next
    }
    
    n <- nrow(f_mat)
    
    if (n <= 0) {
      rm(f_mat)
      next
    }
    
    if (max(sel_cols) > ncol(f_mat)) {
      stop(sprintf(
        "요청 연도열 부족: yToDo=%d, 필요 최대열=%d, f_mat ncol=%d, file=%s",
        yToDo,
        max(sel_cols),
        ncol(f_mat),
        basename(cf)
      ))
    }
    
    end_pos <- min(pos + n - 1L, NumPix)
    n_use <- end_pos - pos + 1L
    
    Fmat[pos:end_pos, ] <- f_mat[1:n_use, sel_cols, drop = FALSE]
    
    pos <- end_pos + 1L
    
    rm(f_mat)
    
    cat("Loaded:", basename(cf), " / filled rows:", pos - 1L, "\n")
    
    if (pos > NumPix) break
  }
  
  if (pos - 1L != NumPix) {
    warning(sprintf("Filled %d of %d pixels. 남은 픽셀은 NA입니다.", pos - 1L, NumPix))
  }
  
  Fmat
}

Fmat <- build_Fmat_prediction(
  FeatDir = FeatDir,
  NumPix = NumPix,
  yToDo = yToDo,
  numFeatUse = numFeatUse
)

## ===== SWIR1 → 3 m 정렬 =====
## 현재 파일명 규칙:
## /home/hjlee/project/mangrove/Img_cliped/hungyen/swir/swir1_yyyy.tif

get_swir_fp <- function(target_year, swirRoot, use_nearest = TRUE) {
  fs <- list.files(
    swirRoot,
    pattern = "^swir1_\\d{4}\\.tif$",
    full.names = TRUE
  )
  
  if (!length(fs)) {
    return(list(fp = NA_character_, year = NA_integer_))
  }
  
  yrs <- as.integer(sub("^swir1_(\\d{4})\\.tif$", "\\1", basename(fs)))
  
  hit <- which(yrs == target_year)
  
  if (length(hit)) {
    return(list(fp = fs[hit[1]], year = yrs[hit[1]]))
  }
  
  if (!use_nearest) {
    return(list(fp = NA_character_, year = NA_integer_))
  }
  
  le <- yrs[yrs <= target_year]
  
  if (length(le)) {
    y <- max(le)
    return(list(fp = fs[which(yrs == y)[1]], year = y))
  }
  
  y <- yrs[which.min(abs(yrs - target_year))]
  list(fp = fs[which(yrs == y)[1]], year = y)
}

swir_info <- get_swir_fp(
  target_year = year_to_predict,
  swirRoot = swirRoot,
  use_nearest = use_nearest_swir
)

if (is.na(swir_info$fp)) {
  warning("SWIR 파일을 찾지 못했습니다: ", swirRoot)
  swir_vec <- rep(NA_real_, NumPix)
  
} else {
  message(sprintf(
    "SWIR 사용 파일: %s (연도=%d)",
    basename(swir_info$fp),
    swir_info$year
  ))
  
  swir <- rast(swir_info$fp)
  
  if (as.character(crs(swir)) != as.character(crs(imgBase))) {
    swir <- project(swir, imgBase, method = "bilinear")
  }
  
  swir3m <- resample(swir, imgBase, method = "bilinear")
  swir_vec <- as.numeric(values(swir3m))
  
  if (length(swir_vec) != NumPix) {
    stop("SWIR 길이 불일치: ", length(swir_vec), " vs ", NumPix)
  }
}

## ===== DEM → 3 m 정렬 =====
## 현재 구조:
## /home/hjlee/project/mangrove/Img_cliped/hungyen/dem/
## 안에 tif 파일 1개

dem_paths <- list.files(
  demRoot,
  pattern = "\\.tif$",
  full.names = TRUE
)

dem_paths <- dem_paths[file.exists(dem_paths)]

if (!length(dem_paths)) {
  warning("DEM 파일을 찾지 못했습니다: ", demRoot)
  dem_vec <- rep(NA_real_, NumPix)
  
} else if (length(dem_paths) == 1) {
  message("DEM 사용 파일: ", basename(dem_paths[1]))
  
  dem_mosaic <- rast(dem_paths[1])
  dem_mosaic[dem_mosaic <= -1000] <- NA
  
  dem_use <- dem_mosaic
  
  if (as.character(crs(dem_use)) != as.character(crs(imgBase))) {
    img_in_dem <- project(imgBase, crs(dem_use))
    dem_use <- crop(dem_use, ext(img_in_dem))
  } else {
    dem_use <- crop(dem_use, ext(imgBase))
  }
  
  dem3m <- project(dem_use, imgBase, method = "bilinear")
  dem_vec <- as.numeric(values(dem3m))
  
  if (length(dem_vec) != NumPix) {
    stop("DEM 길이 불일치: ", length(dem_vec), " vs ", NumPix)
  }
  
} else {
  message("DEM tif가 여러 개 있어 mosaic 처리합니다.")
  
  dem_list <- lapply(dem_paths, rast)
  dem_mosaic <- do.call(mosaic, c(dem_list, fun = "mean"))
  dem_mosaic[dem_mosaic <= -1000] <- NA
  
  dem_use <- dem_mosaic
  
  if (as.character(crs(dem_use)) != as.character(crs(imgBase))) {
    img_in_dem <- project(imgBase, crs(dem_use))
    dem_use <- crop(dem_use, ext(img_in_dem))
  } else {
    dem_use <- crop(dem_use, ext(imgBase))
  }
  
  dem3m <- project(dem_use, imgBase, method = "bilinear")
  dem_vec <- as.numeric(values(dem3m))
  
  if (length(dem_vec) != NumPix) {
    stop("DEM 길이 불일치: ", length(dem_vec), " vs ", NumPix)
  }
  
  rm(dem_list)
}

## ===== 예측 입력 스택 구성: F1..F40 + SWIR1_20m + DEM30m =====
Xpred <- as.data.frame(Fmat, check.names = FALSE)
Xpred$SWIR1_20m <- swir_vec
Xpred$DEM30m    <- dem_vec

rm(Fmat, swir_vec, dem_vec)
gc()

## ===== 모델/전처리 로드 =====
modelPath <- file.path(modelDir, sprintf("rf_%d.rda", year_to_predict))

if (!file.exists(modelPath)) {
  stop("모델 파일 없음: ", modelPath)
}

load(modelPath)  # 기대 객체: rf_all, preproc

model_obj <- if (exists("rf_all")) {
  rf_all
} else {
  stop("모델 객체 rf_all을 찾지 못했습니다.")
}

if (!exists("preproc")) {
  stop("preproc 전처리 정보 객체가 없습니다.")
}

## ===== 예측 전처리: 학습과 동일한 컬럼만 사용 =====
want <- preproc$keep_cols

if (is.null(want)) {
  stop("preproc$keep_cols 가 없습니다.")
}

miss_cols <- setdiff(want, colnames(Xpred))

if (length(miss_cols) > 0) {
  stop("예측 입력에 누락된 컬럼: ", paste(miss_cols, collapse = ", "))
}

Xpred <- Xpred[, want, drop = FALSE]

## 유한성 / 결측 처리
Xmat <- as.matrix(Xpred)
Xmat[!is.finite(Xmat)] <- NA
Xpred <- as.data.frame(Xmat)

## 학습 때 저장한 중앙값으로 NA 대체
if (!is.null(preproc$medians)) {
  for (nm in intersect(names(preproc$medians), colnames(Xpred))) {
    if (is.numeric(Xpred[[nm]])) {
      m <- preproc$medians[[nm]]
      
      if (is.finite(m)) {
        nas <- is.na(Xpred[[nm]])
        
        if (any(nas)) {
          Xpred[[nm]][nas] <- m
        }
      }
    }
  }
}

## 그래도 남은 NA는 randomForest roughfix로 처리
if (anyNA(Xpred)) {
  Xpred <- randomForest::na.roughfix(Xpred)
}

## ===== 예측 & 저장 =====
pred_rf  <- predict(model_obj, Xpred)
pred_map <- setValues(imgBase, as.integer(pred_rf))

outPath <- file.path(outputDir, sprintf("01_map_%d.tif", year_to_predict))

writeRaster(
  pred_map,
  filename = outPath,
  overwrite = TRUE
)

cat("Saved prediction map:", outPath, "\n")