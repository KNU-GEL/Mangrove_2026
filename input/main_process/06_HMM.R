library(terra)

## ============================================================
## HMM Forward-Backward Postprocessing
## 별도 실행용 코드
##
## 입력:
##   output/01_prob_all_classes_YYYY.tif
##
## 출력:
##   output_hmm/02_hmm_map_YYYY.tif
##
## 특징:
##   - Viterbi 사용 안 함
##   - Forward-Backward만 사용
##   - HMM 보정 후 확률 raster는 저장하지 않음
##   - 2018~2025를 한 해씩 처리
##   - block 처리 없음
## ============================================================


## ===== 사용자 설정 =====
site <- "hungyen"

years <- 2018:2025

rootDir <- "/home/hjlee/project/mangrove/Img_cliped"
BaseDir <- file.path(rootDir, site)

## 기존 RF 확률 지도 위치
inputDir <- file.path(BaseDir, "output_balimodel")

## HMM 결과 저장 위치
outputDir <- file.path(BaseDir, "output_hmm_balimodel")

if (!dir.exists(outputDir)) {
  dir.create(outputDir, recursive = TRUE)
}

## 입력 확률 파일명 규칙
prob_file_pattern <- "01_prob_all_classes_%d.tif"

## 출력 HMM map 파일명 규칙
hmm_map_pattern <- "02_hmm_map_%d.tif"


## ===== 클래스 정의 =====
## RF probability tif의 band 순서와 반드시 일치해야 함
## band 1 = mangrove
## band 2 = terrestrial
## band 3 = crop
## band 4 = built-up
## band 5 = deep_water
## band 6 = shallow_water

class_codes <- 1:6

class_names <- c(
  "mangrove",
  "terrestrial",
  "crop",
  "built_up",
  "deep_water",
  "shallow_water"
)

K <- length(class_codes)


## ===== Transition matrix 설정 =====
## p_change = 0.1
## 같은 클래스 유지 확률 = 0.9
## 다른 클래스 전이 확률 = 0.1 / 5 = 0.02

p_change <- 0.1

trans_mat <- matrix(
  p_change / (K - 1),
  nrow = K,
  ncol = K
)

diag(trans_mat) <- 1 - p_change

rownames(trans_mat) <- class_names
colnames(trans_mat) <- class_names

cat("Transition matrix:\n")
print(trans_mat)

write.csv(
  trans_mat,
  file = file.path(outputDir, "hmm_transition_matrix.csv"),
  row.names = TRUE
)

write.csv(
  data.frame(
    class_code = class_codes,
    class_name = class_names,
    probability_band = 1:K
  ),
  file = file.path(outputDir, "hmm_class_order.csv"),
  row.names = FALSE
)


## ===== 입력 확률 파일 확인 =====
prob_paths <- file.path(
  inputDir,
  sprintf(prob_file_pattern, years)
)

names(prob_paths) <- as.character(years)

missing_files <- prob_paths[!file.exists(prob_paths)]

if (length(missing_files) > 0) {
  stop(
    "다음 RF 확률 파일이 없습니다:\n",
    paste(missing_files, collapse = "\n")
  )
}


## ===== raster 객체 열기 =====
prob_rasters <- lapply(prob_paths, function(fp) {
  r <- terra::rast(fp)
  
  if (terra::nlyr(r) < K) {
    stop("확률 raster의 band 수가 6보다 작습니다: ", fp)
  }
  
  r[[1:K]]
})

names(prob_rasters) <- as.character(years)

template <- prob_rasters[[1]][[1]]

## 모든 연도 raster geometry 확인
for (yy in years) {
  ok <- terra::compareGeom(
    template,
    prob_rasters[[as.character(yy)]][[1]],
    stopOnError = FALSE
  )
  
  if (!ok) {
    stop("Raster geometry가 서로 다릅니다. 확인 필요 연도: ", yy)
  }
}

n_pix_total <- terra::ncell(template)

cat("Total pixels:", n_pix_total, "\n")


## ============================================================
## 유틸 함수
## ============================================================

normalize_rows <- function(M, uniform_if_bad = TRUE) {
  M <- as.matrix(M)
  
  M[!is.finite(M)] <- NA_real_
  M[M < 0] <- NA_real_
  
  rs <- rowSums(M, na.rm = TRUE)
  bad <- !is.finite(rs) | rs <= 0
  
  M[is.na(M)] <- 0
  
  if (any(!bad)) {
    M[!bad, ] <- M[!bad, , drop = FALSE] / rs[!bad]
  }
  
  if (any(bad)) {
    if (uniform_if_bad) {
      M[bad, ] <- 1 / ncol(M)
    } else {
      M[bad, ] <- NA_real_
    }
  }
  
  M
}


read_obs_year <- function(r, K) {
  ## 전체 raster 값을 읽습니다.
  ## 반환 구조:
  ##   prob     : ncell x 6 normalized probability matrix
  ##   has_info : 해당 픽셀에 유효한 확률 정보가 있었는지
  
  v <- terra::values(r, mat = TRUE)
  v <- as.matrix(v)
  
  if (ncol(v) < K) {
    stop("읽은 확률 raster의 column 수가 K보다 작습니다.")
  }
  
  v <- v[, 1:K, drop = FALSE]
  
  has_info <- rowSums(is.finite(v) & v >= 0, na.rm = TRUE) > 0 &
    rowSums(v, na.rm = TRUE) > 0
  
  v <- normalize_rows(v, uniform_if_bad = TRUE)
  
  list(
    prob = v,
    has_info = has_info
  )
}


calc_alpha_to_target <- function(
    prob_rasters,
    years,
    target_year,
    trans_mat,
    K
) {
  target_idx <- which(years == target_year)
  
  alpha <- NULL
  has_any_info <- NULL
  
  for (ii in seq_len(target_idx)) {
    yy <- years[ii]
    
    cat("    Forward reading year:", yy, "\n")
    
    obs_obj <- read_obs_year(
      r = prob_rasters[[as.character(yy)]],
      K = K
    )
    
    obs <- obs_obj$prob
    
    if (is.null(has_any_info)) {
      has_any_info <- obs_obj$has_info
    } else {
      has_any_info <- has_any_info | obs_obj$has_info
    }
    
    if (ii == 1) {
      ## 첫 해는 RF observation probability 그대로 시작
      ## uniform initial prior와 같은 효과
      alpha <- normalize_rows(obs, uniform_if_bad = TRUE)
    } else {
      ## alpha_t(k)
      ## = obs_t(k) * sum_i alpha_{t-1}(i) * A(i,k)
      pred <- alpha %*% trans_mat
      alpha <- pred * obs
      alpha <- normalize_rows(alpha, uniform_if_bad = TRUE)
    }
    
    rm(obs_obj, obs)
    gc(verbose = FALSE)
  }
  
  list(
    alpha = alpha,
    has_any_info = has_any_info
  )
}


calc_beta_from_future <- function(
    prob_rasters,
    years,
    target_year,
    trans_mat,
    K,
    n_pix_total
) {
  target_idx <- which(years == target_year)
  n_years <- length(years)
  
  ## 마지막 연도에서는 beta = 1
  beta <- matrix(
    1,
    nrow = n_pix_total,
    ncol = K
  )
  
  has_future_info <- rep(FALSE, n_pix_total)
  
  if (target_idx < n_years) {
    for (ii in seq(from = n_years, to = target_idx + 1, by = -1)) {
      yy <- years[ii]
      
      cat("    Backward reading year:", yy, "\n")
      
      obs_obj <- read_obs_year(
        r = prob_rasters[[as.character(yy)]],
        K = K
      )
      
      obs <- obs_obj$prob
      has_future_info <- has_future_info | obs_obj$has_info
      
      ## beta_t(i)
      ## = sum_j A(i,j) * obs_{t+1}(j) * beta_{t+1}(j)
      beta <- (obs * beta) %*% t(trans_mat)
      beta <- normalize_rows(beta, uniform_if_bad = TRUE)
      
      rm(obs_obj, obs)
      gc(verbose = FALSE)
    }
  }
  
  list(
    beta = beta,
    has_future_info = has_future_info
  )
}


## ============================================================
## Main loop
## 2018년부터 2025년까지 한 해씩 처리
## ============================================================

for (target_year in years) {
  
  cat("\n======================================\n")
  cat("HMM target year:", target_year, "\n")
  cat("======================================\n")
  
  out_fp <- file.path(
    outputDir,
    sprintf(hmm_map_pattern, target_year)
  )
  
  ## ------------------------------------------------------------
  ## 1. Forward 계산
  ## 2018년부터 target_year까지 사용
  ## ------------------------------------------------------------
  alpha_obj <- calc_alpha_to_target(
    prob_rasters = prob_rasters,
    years = years,
    target_year = target_year,
    trans_mat = trans_mat,
    K = K
  )
  
  alpha_target <- alpha_obj$alpha
  has_any_info <- alpha_obj$has_any_info
  
  
  ## ------------------------------------------------------------
  ## 2. Backward 계산
  ## 2025년부터 target_year + 1까지 사용
  ## target_year가 2025이면 beta는 모두 1입니다.
  ## ------------------------------------------------------------
  beta_obj <- calc_beta_from_future(
    prob_rasters = prob_rasters,
    years = years,
    target_year = target_year,
    trans_mat = trans_mat,
    K = K,
    n_pix_total = n_pix_total
  )
  
  beta_target <- beta_obj$beta
  has_any_info <- has_any_info | beta_obj$has_future_info
  
  
  ## ------------------------------------------------------------
  ## 3. Forward와 Backward 결합
  ## gamma_t(k) = alpha_t(k) * beta_t(k)
  ## ------------------------------------------------------------
  gamma <- alpha_target * beta_target
  gamma <- normalize_rows(gamma, uniform_if_bad = TRUE)
  
  
  ## ------------------------------------------------------------
  ## 4. 해당 연도의 최종 label 선택
  ## Viterbi가 아니라 해당 연도 marginal probability의 argmax
  ## ------------------------------------------------------------
  label_out <- max.col(gamma, ties.method = "first")
  
  ## 모든 연도에서 확률 정보가 전혀 없던 픽셀은 NA 처리
  label_out[!has_any_info] <- NA_integer_
  
  if (length(label_out) != n_pix_total) {
    stop(
      "label_out 길이 불일치: ",
      length(label_out),
      " vs ",
      n_pix_total
    )
  }
  
  
  ## ------------------------------------------------------------
  ## 5. raster로 변환 후 저장
  ## ------------------------------------------------------------
  out_r <- terra::setValues(template, label_out)
  names(out_r) <- paste0("hmm_class_", target_year)
  
  terra::writeRaster(
    out_r,
    filename = out_fp,
    overwrite = TRUE,
    datatype = "INT1U",
    NAflag = 0,
    gdal = c("COMPRESS=LZW")
  )
  
  cat("Saved HMM map:", out_fp, "\n")
  
  rm(
    alpha_obj,
    beta_obj,
    alpha_target,
    beta_target,
    gamma,
    label_out,
    out_r
  )
  
  gc(verbose = FALSE)
}

cat("\nDone. HMM-corrected annual maps saved to:\n")
cat(outputDir, "\n")


## ============================================================
## HMM 후처리 평가 지표 저장
## - HMM map 생성이 모두 끝난 뒤 실행
## - 최종 산출물: Excel 파일 1개
##
## 평가 지표:
## A. 평균 라벨 변화 횟수
## B. Stable pixel proportion
## C. % Year-to-Year change
## D. Pairwise label transition frequency
## E. A-B-A flip
## G. Illogical transition
##    - deep_water -> mangrove
##    - deep_water -> terrestrial
##    - deep_water -> crop
##
## 필요 패키지:
##   openxlsx
## ============================================================

cat("\n======================================\n")
cat("Start HMM evaluation metrics\n")
cat("======================================\n")

if (!requireNamespace("openxlsx", quietly = TRUE)) {
  stop(
    "openxlsx 패키지가 필요합니다. 먼저 R에서 install.packages('openxlsx')를 실행해 주세요."
  )
}

## ===== 평가 입력 파일 =====
raw_map_pattern <- "01_map_%d.tif"
hmm_map_pattern <- "02_hmm_map_%d.tif"

raw_paths <- file.path(inputDir, sprintf(raw_map_pattern, years))
hmm_paths <- file.path(outputDir, sprintf(hmm_map_pattern, years))

names(raw_paths) <- as.character(years)
names(hmm_paths) <- as.character(years)

missing_eval_files <- c(
  raw_paths[!file.exists(raw_paths)],
  hmm_paths[!file.exists(hmm_paths)]
)

if (length(missing_eval_files) > 0) {
  stop(
    "평가에 필요한 map 파일이 없습니다:\n",
    paste(missing_eval_files, collapse = "\n")
  )
}

## ===== 평가용 함수 =====

read_map_vec_eval <- function(fp, class_codes) {
  v <- terra::values(terra::rast(fp), mat = FALSE)
  v <- as.integer(v)
  v[!(v %in% class_codes)] <- NA_integer_
  v
}

calc_label_change_count <- function(map_mat) {
  ## map_mat: n_pixels x n_years
  
  n_years <- ncol(map_mat)
  change_count <- rep(NA_integer_, nrow(map_mat))
  
  valid_rows <- rowSums(!is.na(map_mat)) >= 2
  
  if (any(valid_rows)) {
    x <- map_mat[valid_rows, , drop = FALSE]
    
    ch <- matrix(FALSE, nrow = nrow(x), ncol = n_years - 1)
    
    for (i in seq_len(n_years - 1)) {
      ok <- !is.na(x[, i]) & !is.na(x[, i + 1])
      ch[, i] <- ok & x[, i] != x[, i + 1]
    }
    
    change_count[valid_rows] <- rowSums(ch, na.rm = TRUE)
  }
  
  change_count
}

label_change_summary <- function(change_count, product_name) {
  valid <- !is.na(change_count)
  
  data.frame(
    product = product_name,
    valid_pixels = sum(valid),
    mean_label_changes = mean(change_count[valid], na.rm = TRUE),
    median_label_changes = median(change_count[valid], na.rm = TRUE),
    max_label_changes = max(change_count[valid], na.rm = TRUE),
    stable_pixels = sum(change_count[valid] == 0),
    stable_pixel_proportion = sum(change_count[valid] == 0) / sum(valid)
  )
}

label_change_histogram <- function(change_count, product_name) {
  valid <- !is.na(change_count)
  max_ch <- max(change_count[valid], na.rm = TRUE)
  
  tab <- table(
    factor(change_count[valid], levels = 0:max_ch),
    useNA = "no"
  )
  
  data.frame(
    product = product_name,
    label_change_count = as.integer(names(tab)),
    pixel_count = as.integer(tab),
    proportion = as.integer(tab) / sum(tab)
  )
}

annual_change_summary_eval <- function(v1, v2, product_name, year_from, year_to) {
  ok <- !is.na(v1) & !is.na(v2)
  
  n_valid <- sum(ok)
  n_change <- sum(v1[ok] != v2[ok])
  
  data.frame(
    product = product_name,
    year_from = year_from,
    year_to = year_to,
    valid_pixels = n_valid,
    changed_pixels = n_change,
    unchanged_pixels = n_valid - n_change,
    change_rate = ifelse(n_valid > 0, n_change / n_valid, NA_real_)
  )
}

transition_table_eval <- function(v1, v2, product_name, year_from, year_to,
                                  class_codes, class_names) {
  ok <- !is.na(v1) & !is.na(v2)
  
  from_fac <- factor(v1[ok], levels = class_codes)
  to_fac   <- factor(v2[ok], levels = class_codes)
  
  tab <- as.data.frame(table(from_fac, to_fac))
  colnames(tab) <- c("from_code", "to_code", "count")
  
  tab$from_code <- as.integer(as.character(tab$from_code))
  tab$to_code   <- as.integer(as.character(tab$to_code))
  
  tab$product <- product_name
  tab$year_from <- year_from
  tab$year_to <- year_to
  
  tab$from_class <- class_names[match(tab$from_code, class_codes)]
  tab$to_class   <- class_names[match(tab$to_code, class_codes)]
  
  tab <- tab[, c(
    "product",
    "year_from",
    "year_to",
    "from_code",
    "from_class",
    "to_code",
    "to_class",
    "count"
  )]
  
  tab
}

aba_flip_count_eval <- function(v_prev, v_mid, v_next,
                                product_name,
                                year_prev, year_mid, year_next) {
  ok <- !is.na(v_prev) & !is.na(v_mid) & !is.na(v_next)
  
  aba <- ok &
    v_prev == v_next &
    v_prev != v_mid
  
  data.frame(
    product = product_name,
    year_prev = year_prev,
    year_mid = year_mid,
    year_next = year_next,
    valid_pixels = sum(ok),
    aba_flip_pixels = sum(aba),
    aba_flip_rate = ifelse(sum(ok) > 0, sum(aba) / sum(ok), NA_real_)
  )
}

count_illogical_eval <- function(v1, v2, product_name, year_from, year_to,
                                 illogical_pairs) {
  ok <- !is.na(v1) & !is.na(v2)
  
  n_valid <- sum(ok)
  n_changed <- sum(v1[ok] != v2[ok])
  
  if (n_valid == 0) {
    return(data.frame(
      product = product_name,
      year_from = year_from,
      year_to = year_to,
      valid_pixels = 0,
      changed_pixels = 0,
      illogical_pixels = 0,
      illogical_rate_all = NA_real_,
      illogical_share_of_changes = NA_real_
    ))
  }
  
  pair_key <- paste(v1[ok], v2[ok], sep = "_")
  ill_key <- paste(illogical_pairs$from_code, illogical_pairs$to_code, sep = "_")
  
  n_ill <- sum(pair_key %in% ill_key)
  
  data.frame(
    product = product_name,
    year_from = year_from,
    year_to = year_to,
    valid_pixels = n_valid,
    changed_pixels = n_changed,
    illogical_pixels = n_ill,
    illogical_rate_all = n_ill / n_valid,
    illogical_share_of_changes = ifelse(n_changed > 0, n_ill / n_changed, NA_real_)
  )
}

illogical_by_pair_eval <- function(v1, v2, product_name, year_from, year_to,
                                   illogical_pairs, class_codes, class_names) {
  ok <- !is.na(v1) & !is.na(v2)
  
  out <- illogical_pairs
  
  out$product <- product_name
  out$year_from <- year_from
  out$year_to <- year_to
  
  out$from_class <- class_names[match(out$from_code, class_codes)]
  out$to_class   <- class_names[match(out$to_code, class_codes)]
  
  if (sum(ok) == 0) {
    out$count <- 0L
  } else {
    pair_key <- paste(v1[ok], v2[ok], sep = "_")
    
    out$count <- vapply(
      seq_len(nrow(out)),
      function(i) {
        sum(pair_key == paste(out$from_code[i], out$to_code[i], sep = "_"))
      },
      integer(1)
    )
  }
  
  out <- out[, c(
    "product",
    "year_from",
    "year_to",
    "from_code",
    "from_class",
    "to_code",
    "to_class",
    "count"
  )]
  
  out
}

## ============================================================
## 불가능 전이 정의
## class code:
## 1 = mangrove
## 2 = terrestrial
## 3 = crop
## 4 = built_up / bare land 포함
## 5 = deep_water
## 6 = shallow_water
##
## 선생님 설정:
##   deep_water -> mangrove
##   deep_water -> terrestrial
##   deep_water -> crop
## 만 불가능 전이로 계산
## ============================================================

illogical_pairs <- data.frame(
  from_code = c(5, 5, 5),
  to_code   = c(1, 2, 3)
)

illogical_pairs$from_class <- class_names[match(illogical_pairs$from_code, class_codes)]
illogical_pairs$to_class   <- class_names[match(illogical_pairs$to_code, class_codes)]

## ============================================================
## 1. 전체 연도 map을 matrix로 읽기
## ============================================================

cat("\nReading annual maps...\n")

raw_mat <- NULL
hmm_mat <- NULL

for (yy in years) {
  cat("  reading year:", yy, "\n")
  
  raw_v <- read_map_vec_eval(raw_paths[as.character(yy)], class_codes)
  hmm_v <- read_map_vec_eval(hmm_paths[as.character(yy)], class_codes)
  
  if (is.null(raw_mat)) {
    raw_mat <- matrix(NA_integer_, nrow = length(raw_v), ncol = length(years))
    hmm_mat <- matrix(NA_integer_, nrow = length(hmm_v), ncol = length(years))
    colnames(raw_mat) <- as.character(years)
    colnames(hmm_mat) <- as.character(years)
  }
  
  raw_mat[, as.character(yy)] <- raw_v
  hmm_mat[, as.character(yy)] <- hmm_v
  
  rm(raw_v, hmm_v)
  gc(verbose = FALSE)
}

## ============================================================
## A. Number of label changes per pixel
## B. Stable pixel proportion
## ============================================================

cat("\nCalculating label change count...\n")

raw_change_count <- calc_label_change_count(raw_mat)
hmm_change_count <- calc_label_change_count(hmm_mat)

label_change_summary_df <- rbind(
  label_change_summary(raw_change_count, "before_rf"),
  label_change_summary(hmm_change_count, "after_hmm")
)

label_change_histogram_df <- rbind(
  label_change_histogram(raw_change_count, "before_rf"),
  label_change_histogram(hmm_change_count, "after_hmm")
)

stable_pixel_summary_df <- label_change_summary_df[, c(
  "product",
  "valid_pixels",
  "stable_pixels",
  "stable_pixel_proportion"
)]

## ============================================================
## C. % Year-to-Year change
## D. Pairwise label transition frequency
## G. Illogical transition
## ============================================================

cat("\nCalculating annual changes, transitions, and illogical transitions...\n")

annual_change_list <- list()
transition_list <- list()
illogical_summary_list <- list()
illogical_by_pair_list <- list()

for (i in seq_len(length(years) - 1)) {
  y1 <- years[i]
  y2 <- years[i + 1]
  
  cat("  transition:", y1, "->", y2, "\n")
  
  raw_1 <- raw_mat[, as.character(y1)]
  raw_2 <- raw_mat[, as.character(y2)]
  
  hmm_1 <- hmm_mat[, as.character(y1)]
  hmm_2 <- hmm_mat[, as.character(y2)]
  
  annual_change_list[[length(annual_change_list) + 1]] <- annual_change_summary_eval(
    raw_1,
    raw_2,
    "before_rf",
    y1,
    y2
  )
  
  annual_change_list[[length(annual_change_list) + 1]] <- annual_change_summary_eval(
    hmm_1,
    hmm_2,
    "after_hmm",
    y1,
    y2
  )
  
  transition_list[[length(transition_list) + 1]] <- transition_table_eval(
    raw_1,
    raw_2,
    "before_rf",
    y1,
    y2,
    class_codes,
    class_names
  )
  
  transition_list[[length(transition_list) + 1]] <- transition_table_eval(
    hmm_1,
    hmm_2,
    "after_hmm",
    y1,
    y2,
    class_codes,
    class_names
  )
  
  illogical_summary_list[[length(illogical_summary_list) + 1]] <- count_illogical_eval(
    raw_1,
    raw_2,
    "before_rf",
    y1,
    y2,
    illogical_pairs
  )
  
  illogical_summary_list[[length(illogical_summary_list) + 1]] <- count_illogical_eval(
    hmm_1,
    hmm_2,
    "after_hmm",
    y1,
    y2,
    illogical_pairs
  )
  
  illogical_by_pair_list[[length(illogical_by_pair_list) + 1]] <- illogical_by_pair_eval(
    raw_1,
    raw_2,
    "before_rf",
    y1,
    y2,
    illogical_pairs,
    class_codes,
    class_names
  )
  
  illogical_by_pair_list[[length(illogical_by_pair_list) + 1]] <- illogical_by_pair_eval(
    hmm_1,
    hmm_2,
    "after_hmm",
    y1,
    y2,
    illogical_pairs,
    class_codes,
    class_names
  )
}

annual_change_df <- do.call(rbind, annual_change_list)
transition_counts_df <- do.call(rbind, transition_list)
illogical_summary_df <- do.call(rbind, illogical_summary_list)
illogical_by_pair_df <- do.call(rbind, illogical_by_pair_list)

## ============================================================
## E. A-B-A flip
## 예: 2018=A, 2019=B, 2020=A
## ============================================================

cat("\nCalculating A-B-A flips...\n")

aba_list <- list()

if (length(years) >= 3) {
  for (i in 2:(length(years) - 1)) {
    y_prev <- years[i - 1]
    y_mid  <- years[i]
    y_next <- years[i + 1]
    
    cat("  A-B-A:", y_prev, y_mid, y_next, "\n")
    
    raw_prev <- raw_mat[, as.character(y_prev)]
    raw_mid  <- raw_mat[, as.character(y_mid)]
    raw_next <- raw_mat[, as.character(y_next)]
    
    hmm_prev <- hmm_mat[, as.character(y_prev)]
    hmm_mid  <- hmm_mat[, as.character(y_mid)]
    hmm_next <- hmm_mat[, as.character(y_next)]
    
    aba_list[[length(aba_list) + 1]] <- aba_flip_count_eval(
      raw_prev,
      raw_mid,
      raw_next,
      "before_rf",
      y_prev,
      y_mid,
      y_next
    )
    
    aba_list[[length(aba_list) + 1]] <- aba_flip_count_eval(
      hmm_prev,
      hmm_mid,
      hmm_next,
      "after_hmm",
      y_prev,
      y_mid,
      y_next
    )
  }
}

aba_flip_df <- do.call(rbind, aba_list)

## ============================================================
## 전체 요약 생성
## ============================================================

cat("\nCreating total summary...\n")

mean_annual_change <- aggregate(
  change_rate ~ product,
  data = annual_change_df,
  FUN = mean,
  na.rm = TRUE
)

total_changed <- aggregate(
  changed_pixels ~ product,
  data = annual_change_df,
  FUN = sum,
  na.rm = TRUE
)

mean_aba <- aggregate(
  aba_flip_rate ~ product,
  data = aba_flip_df,
  FUN = mean,
  na.rm = TRUE
)

total_aba <- aggregate(
  aba_flip_pixels ~ product,
  data = aba_flip_df,
  FUN = sum,
  na.rm = TRUE
)

mean_illogical_rate <- aggregate(
  illogical_rate_all ~ product,
  data = illogical_summary_df,
  FUN = mean,
  na.rm = TRUE
)

mean_illogical_share <- aggregate(
  illogical_share_of_changes ~ product,
  data = illogical_summary_df,
  FUN = mean,
  na.rm = TRUE
)

total_illogical <- aggregate(
  illogical_pixels ~ product,
  data = illogical_summary_df,
  FUN = sum,
  na.rm = TRUE
)

eval_total_summary_df <- Reduce(
  function(x, y) merge(x, y, by = "product", all = TRUE),
  list(
    label_change_summary_df,
    mean_annual_change,
    total_changed,
    mean_aba,
    total_aba,
    mean_illogical_rate,
    mean_illogical_share,
    total_illogical
  )
)

## 이름 조금 명확하게 변경
names(eval_total_summary_df)[names(eval_total_summary_df) == "change_rate"] <- "mean_year_to_year_change_rate"
names(eval_total_summary_df)[names(eval_total_summary_df) == "changed_pixels"] <- "total_changed_pixels"
names(eval_total_summary_df)[names(eval_total_summary_df) == "aba_flip_rate"] <- "mean_aba_flip_rate"
names(eval_total_summary_df)[names(eval_total_summary_df) == "aba_flip_pixels"] <- "total_aba_flip_pixels"
names(eval_total_summary_df)[names(eval_total_summary_df) == "illogical_rate_all"] <- "mean_illogical_rate_all"
names(eval_total_summary_df)[names(eval_total_summary_df) == "illogical_share_of_changes"] <- "mean_illogical_share_of_changes"
names(eval_total_summary_df)[names(eval_total_summary_df) == "illogical_pixels"] <- "total_illogical_pixels"

## ============================================================
## Excel 파일 하나로 저장
## ============================================================

cat("\nWriting Excel workbook...\n")

xlsx_path <- file.path(outputDir, "eval_hmm_summary.xlsx")

wb <- openxlsx::createWorkbook()

add_sheet <- function(wb, sheet_name, df) {
  openxlsx::addWorksheet(wb, sheet_name)
  openxlsx::writeData(wb, sheet_name, df)
  openxlsx::freezePane(wb, sheet_name, firstRow = TRUE)
  openxlsx::setColWidths(wb, sheet_name, cols = 1:ncol(df), widths = "auto")
}

add_sheet(wb, "total_summary", eval_total_summary_df)
add_sheet(wb, "label_change_summary", label_change_summary_df)
add_sheet(wb, "label_change_histogram", label_change_histogram_df)
add_sheet(wb, "stable_pixel_summary", stable_pixel_summary_df)
add_sheet(wb, "annual_change_summary", annual_change_df)
add_sheet(wb, "transition_counts", transition_counts_df)
add_sheet(wb, "aba_flip_summary", aba_flip_df)
add_sheet(wb, "illogical_definition", illogical_pairs)
add_sheet(wb, "illogical_summary", illogical_summary_df)
add_sheet(wb, "illogical_by_pair", illogical_by_pair_df)

openxlsx::saveWorkbook(
  wb,
  file = xlsx_path,
  overwrite = TRUE
)

cat("\nDone. Evaluation Excel saved to:\n")
cat(xlsx_path, "\n")

## 메모리 정리
rm(
  raw_mat,
  hmm_mat,
  raw_change_count,
  hmm_change_count,
  wb
)

gc(verbose = FALSE)

