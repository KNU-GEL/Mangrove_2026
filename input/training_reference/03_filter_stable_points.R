############################################################
# Site random points binary filtering using
# CCDC-style harmonic seasonal break detection
#
# Binary classes:
#   changed
#   no_change
#
# Main idea:
#   - Fit one harmonic seasonal model for the full time series
#   - Search for the best breakpoint
#   - Fit separate harmonic models before and after the break
#   - If the two-segment model is much better AND the seasonal/level
#     structure changes persistently, classify as changed
#
# Usage:
#   Rscript this_script.R gilimanuk
#   Rscript this_script.R sanur
#   Rscript this_script.R jembrana
#   Rscript this_script.R lembongan
############################################################

library(sf)
library(raster)
library(dplyr)
library(lubridate)
library(rjson)

############################################################
# 0. Settings
############################################################

args <- commandArgs(trailingOnly = TRUE)

if (length(args) >= 1) {
  site <- args[1]
} else {
  site <- "aceh"
}

params <- fromJSON(file = "/home/hjlee/project/mangrove/code/PLCM_Parameters.json")

base_out_dir <- params$setup$outDir
numChunks <- params$setup$numChunks

site_periods <- list(
  # sanur     = c("1721", "2024", "2326"),
  # gilimanuk = c("1722", "2126"),
  # jembrana  = c("1722", "2126"),
  # lembongan = c("1722", "2126"),
  # hungyen     = c("1721", "2024", "2326"),
  # taklong = c("1722", "2126"),
  aceh = c("1722", "2126")
)

if (!(site %in% names(site_periods))) {
  stop("Unknown site: ", site)
}

site_dir <- file.path(base_out_dir, site)

gpkg_path <- file.path(site_dir, paste0(site, "_random.gpkg"))
base_img_path <- file.path(site_dir, "base_image.tif")

periods <- site_periods[[site]]

analysis_start <- as.Date("2017-07-01")
analysis_end   <- as.Date("2026-06-30")

out_dir <- file.path(site_dir, paste0(site, "_change_filter"))

if (!dir.exists(out_dir)) {
  dir.create(out_dir, recursive = TRUE)
}

out_all_gpkg <- file.path(out_dir, paste0(site, "_harmonic_change_status.gpkg"))
out_all_csv  <- file.path(out_dir, paste0(site, "_harmonic_change_status.csv"))

out_changed_gpkg <- file.path(out_dir, paste0(site, "_changed_points.gpkg"))
out_changed_csv  <- file.path(out_dir, paste0(site, "_changed_points.csv"))

out_nochange_gpkg <- file.path(out_dir, paste0(site, "_stable_points.gpkg"))
out_nochange_csv  <- file.path(out_dir, paste0(site, "_stable_points.csv"))

out_ts_rds <- file.path(out_dir, paste0(site, "_daily_timeseries.rds"))
out_monthly_rds <- file.path(out_dir, paste0(site, "_monthly_evi2.rds"))

out_pdf_changed <- file.path(out_dir, paste0(site, "_changed_points_timeseries_all.pdf"))
out_pdf_nochange <- file.path(out_dir, paste0(site, "_stable_points_timeseries_all.pdf"))

cat("Site:", site, "\n")
cat("Site directory:", site_dir, "\n")
cat("Input GPKG:", gpkg_path, "\n")
cat("Output directory:", out_dir, "\n")
cat("Periods:", paste(periods, collapse = ", "), "\n")

if (!file.exists(gpkg_path)) {
  stop("Input GPKG does not exist: ", gpkg_path)
}

if (!file.exists(base_img_path)) {
  stop("Base image does not exist: ", base_img_path)
}

############################################################
# 1. Sensitivity settings
############################################################

# QC only. Low data quality will be recorded but not automatically changed.
min_valid_months <- 60
max_missing_gap_months <- 10

# Harmonic model settings
n_harmonics <- 2

# Break search settings
min_segment_months <- 24
min_segment_valid_months <- 16
candidate_step_months <- 3

# Model improvement threshold
# Larger = less sensitive
min_bic_gain <- 10

# Moderate persistent level shift.
# BIC gain이 충분히 크고, 전후 평균 EVI2 level이 이 값 이상 바뀌면 changed로 판단
moderate_level_shift_thr <- 0.12

# Persistent structural change thresholds
# Level shift: difference between before/after harmonic curves
level_shift_thr <- 0.12

# Shape RMSE: difference in seasonal shape after removing mean level
shape_rmse_thr <- 0.05

# Very strong level shift can be changed even if shape change is weak
strong_level_shift_thr <- 0.15

# Pattern present/absent thresholds from harmonic curve amplitude
# Hysteresis prevents unstable toggling around one threshold.
pattern_absent_amp_thr <- 0.08
pattern_present_amp_thr <- 0.14

# EVI2 valid range
evi2_min_allowed <- -0.2
evi2_max_allowed <- 1.2

############################################################
# 2. Helper functions
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


safe_median <- function(x) {
  if (all(is.na(x))) return(NA_real_)
  median(x, na.rm = TRUE)
}


safe_quantile <- function(x, p) {
  if (all(is.na(x))) return(NA_real_)
  as.numeric(quantile(x, probs = p, na.rm = TRUE, names = FALSE))
}


max_consecutive_na <- function(x) {
  z <- is.na(x)
  if (!any(z)) return(0)
  r <- rle(z)
  max(r$lengths[r$values])
}


make_harmonic_df <- function(month_vec, y = NULL, n_harmonics = 2) {
  
  t_year <- as.numeric(month_vec - analysis_start) / 365.25
  
  out <- data.frame(
    month = month_vec,
    t_year = t_year
  )
  
  for (k in seq_len(n_harmonics)) {
    out[[paste0("sin", k)]] <- sin(2 * pi * k * t_year)
    out[[paste0("cos", k)]] <- cos(2 * pi * k * t_year)
  }
  
  if (!is.null(y)) {
    out$evi2 <- y
  }
  
  out
}


harmonic_formula <- function(n_harmonics = 2) {
  rhs <- c("t_year")
  
  for (k in seq_len(n_harmonics)) {
    rhs <- c(rhs, paste0("sin", k), paste0("cos", k))
  }
  
  as.formula(paste("evi2 ~", paste(rhs, collapse = " + ")))
}


fit_harmonic_model <- function(d, n_harmonics = 2, min_obs = 16) {
  
  d <- d[!is.na(d$evi2), ]
  
  min_coef <- 2 + 2 * n_harmonics
  
  if (nrow(d) < max(min_obs, min_coef + 2)) {
    return(NULL)
  }
  
  fit <- try(
    lm(harmonic_formula(n_harmonics), data = d),
    silent = TRUE
  )
  
  if (inherits(fit, "try-error")) {
    return(NULL)
  }
  
  fit
}


model_rss <- function(fit, d) {
  d <- d[!is.na(d$evi2), ]
  
  if (is.null(fit) || nrow(d) == 0) return(NA_real_)
  
  pred <- try(predict(fit, newdata = d), silent = TRUE)
  
  if (inherits(pred, "try-error")) return(NA_real_)
  
  sum((d$evi2 - pred)^2, na.rm = TRUE)
}


calc_bic <- function(rss, n, k) {
  if (is.na(rss) || !is.finite(rss) || n <= 0) return(NA_real_)
  
  rss <- max(rss, 1e-12)
  
  n * log(rss / n) + k * log(n)
}


curve_metrics <- function(fit_before, fit_after, break_month, n_harmonics = 2) {
  
  pred_months <- seq(
    break_month %m-% months(5),
    by = "month",
    length.out = 12
  )
  
  newd <- make_harmonic_df(pred_months, n_harmonics = n_harmonics)
  
  pred_before <- try(predict(fit_before, newdata = newd), silent = TRUE)
  pred_after  <- try(predict(fit_after, newdata = newd), silent = TRUE)
  
  if (inherits(pred_before, "try-error") || inherits(pred_after, "try-error")) {
    return(list(
      level_shift = NA_real_,
      shape_rmse = NA_real_,
      amp_before = NA_real_,
      amp_after = NA_real_,
      amp_change = NA_real_,
      pattern_transition = FALSE,
      pattern_transition_type = NA_character_
    ))
  }
  
  level_shift <- abs(mean(pred_after, na.rm = TRUE) - mean(pred_before, na.rm = TRUE))
  
  cb <- pred_before - mean(pred_before, na.rm = TRUE)
  ca <- pred_after  - mean(pred_after, na.rm = TRUE)
  
  shape_rmse <- sqrt(mean((ca - cb)^2, na.rm = TRUE))
  
  amp_before <- safe_quantile(pred_before, 0.90) - safe_quantile(pred_before, 0.10)
  amp_after  <- safe_quantile(pred_after, 0.90) - safe_quantile(pred_after, 0.10)
  amp_change <- abs(amp_after - amp_before)
  
  before_absent <- amp_before <= pattern_absent_amp_thr
  before_present <- amp_before >= pattern_present_amp_thr
  
  after_absent <- amp_after <= pattern_absent_amp_thr
  after_present <- amp_after >= pattern_present_amp_thr
  
  pattern_transition <- FALSE
  pattern_transition_type <- NA_character_
  
  if (before_absent && after_present) {
    pattern_transition <- TRUE
    pattern_transition_type <- "pattern_appeared"
  }
  
  if (before_present && after_absent) {
    pattern_transition <- TRUE
    pattern_transition_type <- "pattern_disappeared"
  }
  
  list(
    level_shift = level_shift,
    shape_rmse = shape_rmse,
    amp_before = amp_before,
    amp_after = amp_after,
    amp_change = amp_change,
    pattern_transition = pattern_transition,
    pattern_transition_type = pattern_transition_type
  )
}


find_best_harmonic_break <- function(x) {
  
  x_model <- make_harmonic_df(
    month_vec = x$month,
    y = x$evi2,
    n_harmonics = n_harmonics
  )
  
  valid_months <- sum(!is.na(x_model$evi2))
  
  if (valid_months < min_valid_months) {
    return(NULL)
  }
  
  fit_full <- fit_harmonic_model(
    x_model,
    n_harmonics = n_harmonics,
    min_obs = min_valid_months
  )
  
  if (is.null(fit_full)) {
    return(NULL)
  }
  
  n_valid <- sum(!is.na(x_model$evi2))
  k_full <- length(coef(fit_full))
  rss_full <- model_rss(fit_full, x_model)
  bic_full <- calc_bic(rss_full, n_valid, k_full)
  
  if (is.na(bic_full)) {
    return(NULL)
  }
  
  n_months <- nrow(x_model)
  
  candidate_idx <- seq(
    min_segment_months,
    n_months - min_segment_months,
    by = candidate_step_months
  )
  
  if (length(candidate_idx) == 0) {
    return(NULL)
  }
  
  best <- NULL
  
  for (idx in candidate_idx) {
    
    d_before <- x_model[seq_len(idx), ]
    d_after  <- x_model[(idx + 1):n_months, ]
    
    if (sum(!is.na(d_before$evi2)) < min_segment_valid_months) next
    if (sum(!is.na(d_after$evi2)) < min_segment_valid_months) next
    
    fit_before <- fit_harmonic_model(
      d_before,
      n_harmonics = n_harmonics,
      min_obs = min_segment_valid_months
    )
    
    fit_after <- fit_harmonic_model(
      d_after,
      n_harmonics = n_harmonics,
      min_obs = min_segment_valid_months
    )
    
    if (is.null(fit_before) || is.null(fit_after)) next
    
    rss_before <- model_rss(fit_before, d_before)
    rss_after  <- model_rss(fit_after, d_after)
    
    rss_two <- rss_before + rss_after
    k_two <- length(coef(fit_before)) + length(coef(fit_after)) + 1
    bic_two <- calc_bic(rss_two, n_valid, k_two)
    
    if (is.na(bic_two)) next
    
    bic_gain <- bic_full - bic_two
    
    cm <- curve_metrics(
      fit_before = fit_before,
      fit_after = fit_after,
      break_month = x_model$month[idx],
      n_harmonics = n_harmonics
    )
    
    candidate <- list(
      break_index = idx,
      break_month = x_model$month[idx],
      bic_full = bic_full,
      bic_two = bic_two,
      bic_gain = bic_gain,
      rss_full = rss_full,
      rss_two = rss_two,
      fit_full = fit_full,
      fit_before = fit_before,
      fit_after = fit_after,
      level_shift = cm$level_shift,
      shape_rmse = cm$shape_rmse,
      amp_before = cm$amp_before,
      amp_after = cm$amp_after,
      amp_change = cm$amp_change,
      pattern_transition = cm$pattern_transition,
      pattern_transition_type = cm$pattern_transition_type
    )
    
    if (is.null(best) || candidate$bic_gain > best$bic_gain) {
      best <- candidate
    }
  }
  
  best
}


extract_timeseries_from_chunks <- function(points_df) {
  
  out_list <- list()
  counter <- 1
  
  for (period in periods) {
    
    chunk_dir <- file.path(site_dir, paste0("chunk_", period))
    cat("Reading period:", period, "\n")
    
    for (ck in sort(unique(points_df$chunk_id))) {
      
      idx <- which(points_df$chunk_id == ck)
      if (length(idx) == 0) next
      
      ck_num <- sprintf("%03d", ck)
      chunk_file <- file.path(chunk_dir, paste0("chunk_", ck_num, ".rda"))
      
      if (!file.exists(chunk_file)) {
        warning("Missing chunk file: ", chunk_file)
        next
      }
      
      e <- new.env()
      load(chunk_file, envir = e)
      
      rows <- points_df$row_in_chunk[idx]
      
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
      
      blue  <- e$band1[rows, , drop = FALSE] / 10000
      green <- e$band2[rows, , drop = FALSE] / 10000
      red   <- e$band3[rows, , drop = FALSE] / 10000
      nir   <- e$band4[rows, , drop = FALSE] / 10000
      
      ndvi <- (nir - red) / (nir + red)
      gndvi <- (nir - green) / (nir + green)
      nirv <- nir * ndvi
      evi2 <- 2.5 * (nir - red) / (nir + 2.4 * red + 1)
      
      tmp <- data.frame(
        sample_id = rep(points_df$sample_id[idx], each = n_dates),
        period = period,
        chunk_id = rep(points_df$chunk_id[idx], each = n_dates),
        row_in_chunk = rep(points_df$row_in_chunk[idx], each = n_dates),
        date = rep(dates, times = n_pts),
        blue = as.vector(t(blue)),
        green = as.vector(t(green)),
        red = as.vector(t(red)),
        nir = as.vector(t(nir)),
        ndvi = as.vector(t(ndvi)),
        gndvi = as.vector(t(gndvi)),
        nirv = as.vector(t(nirv)),
        evi2 = as.vector(t(evi2))
      )
      
      out_list[[counter]] <- tmp
      counter <- counter + 1
      
      rm(e, blue, green, red, nir, ndvi, gndvi, nirv, evi2, tmp)
      gc()
    }
  }
  
  bind_rows(out_list)
}


classify_one_point_harmonic <- function(d, sid) {
  
  full_months <- seq(
    floor_date(analysis_start, "month"),
    floor_date(analysis_end, "month"),
    by = "month"
  )
  
  x <- data.frame(month = full_months) %>%
    left_join(d, by = "month") %>%
    arrange(month)
  
  valid_months <- sum(!is.na(x$evi2))
  max_gap <- max_consecutive_na(x$evi2)
  
  qc_reasons <- c()
  
  if (valid_months < min_valid_months) {
    qc_reasons <- c(qc_reasons, "too_few_valid_months")
  }
  
  if (max_gap > max_missing_gap_months) {
    qc_reasons <- c(qc_reasons, "long_missing_gap")
  }
  
  qc_flag <- length(qc_reasons) > 0
  qc_reason <- ifelse(qc_flag, paste(qc_reasons, collapse = ";"), "ok")
  
  best <- find_best_harmonic_break(x)
  
  if (is.null(best)) {
    return(data.frame(
      sample_id = sid,
      change_binary = "no_change",
      reason = "insufficient_model_evidence",
      qc_flag = qc_flag,
      qc_reason = qc_reason,
      valid_months = valid_months,
      max_missing_gap_months = max_gap,
      break_month = as.Date(NA),
      bic_gain = NA_real_,
      level_shift = NA_real_,
      shape_rmse = NA_real_,
      amp_before = NA_real_,
      amp_after = NA_real_,
      amp_change = NA_real_,
      pattern_transition = FALSE,
      pattern_transition_type = NA_character_,
      model_improved = FALSE,
      changed_by_level_shape = FALSE,
      changed_by_persistent_level = FALSE,
      changed_by_strong_level = FALSE,
      changed_by_pattern = FALSE
    ))
  }
  
  model_improved <- !is.na(best$bic_gain) &&
    best$bic_gain >= min_bic_gain
  
  changed_by_level_shape <- model_improved &&
    !is.na(best$level_shift) &&
    !is.na(best$shape_rmse) &&
    best$level_shift >= level_shift_thr &&
    best$shape_rmse >= shape_rmse_thr
  
  changed_by_persistent_level <- model_improved &&
    !is.na(best$level_shift) &&
    best$level_shift >= moderate_level_shift_thr
  
  changed_by_strong_level <- model_improved &&
    !is.na(best$level_shift) &&
    best$level_shift >= strong_level_shift_thr
  
  changed_by_pattern <- model_improved &&
    isTRUE(best$pattern_transition)
  
  is_changed <- changed_by_level_shape ||
    changed_by_persistent_level ||
    changed_by_strong_level ||
    changed_by_pattern
  
  reason_vec <- c()
  
  if (changed_by_level_shape) {
    reason_vec <- c(reason_vec, "persistent_level_and_seasonal_shape_change")
  }
  
  if (changed_by_persistent_level) {
    reason_vec <- c(reason_vec, "persistent_level_shift")
  }
  
  if (changed_by_strong_level) {
    reason_vec <- c(reason_vec, "strong_persistent_level_shift")
  }
  
  if (changed_by_pattern) {
    reason_vec <- c(reason_vec, best$pattern_transition_type)
  }
  
  if (!is_changed) {
    if (!model_improved) {
      reason_vec <- "no_significant_harmonic_break"
    } else {
      reason_vec <- "harmonic_break_but_pattern_consistent"
    }
  }
  
  data.frame(
    sample_id = sid,
    change_binary = ifelse(is_changed, "changed", "no_change"),
    reason = paste(unique(reason_vec), collapse = ";"),
    qc_flag = qc_flag,
    qc_reason = qc_reason,
    valid_months = valid_months,
    max_missing_gap_months = max_gap,
    break_month = best$break_month,
    bic_gain = best$bic_gain,
    level_shift = best$level_shift,
    shape_rmse = best$shape_rmse,
    amp_before = best$amp_before,
    amp_after = best$amp_after,
    amp_change = best$amp_change,
    pattern_transition = best$pattern_transition,
    pattern_transition_type = best$pattern_transition_type,
    model_improved = model_improved,
    changed_by_level_shape = changed_by_level_shape,
    changed_by_persistent_level = changed_by_persistent_level,
    changed_by_strong_level = changed_by_strong_level,
    changed_by_pattern = changed_by_pattern
  )
}


############################################################
# Clean plot function
############################################################

plot_one_point_base <- function(sample_id_value, ts_daily, ts_monthly, result_df) {
  
  d0 <- ts_daily[ts_daily$sample_id == sample_id_value, ]
  m0 <- ts_monthly[ts_monthly$sample_id == sample_id_value, ]
  r0 <- result_df[result_df$sample_id == sample_id_value, ]
  
  fmt_num <- function(x, digits = 3) {
    if (length(x) == 0 || is.na(x) || !is.finite(x)) return("NA")
    format(round(x, digits), nsmall = digits)
  }
  
  fmt_chr <- function(x) {
    if (length(x) == 0 || is.na(x) || x == "") return("NA")
    as.character(x)
  }
  
  fmt_date <- function(x) {
    if (length(x) == 0 || is.na(x)) return("NA")
    as.character(as.Date(x))
  }
  
  pass_fail <- function(value, threshold, direction = "ge") {
    if (is.na(value) || !is.finite(value)) return("NA")
    
    if (direction == "ge") {
      return(ifelse(value >= threshold, "Y", "N"))
    }
    
    if (direction == "le") {
      return(ifelse(value <= threshold, "Y", "N"))
    }
    
    "NA"
  }
  
  bool_text <- function(x) {
    if (length(x) == 0 || is.na(x)) return("NA")
    ifelse(isTRUE(x), "TRUE", "FALSE")
  }
  
  draw_text_lines <- function(lines, x = 0.02, y_start = 0.96, line_gap = 0.12, cex = 0.78) {
    y <- y_start
    
    for (ln in lines) {
      wrapped <- strwrap(ln, width = 150)
      
      for (w in wrapped) {
        text(
          x = x,
          y = y,
          labels = w,
          adj = c(0, 1),
          cex = cex,
          family = "mono"
        )
        y <- y - line_gap
      }
    }
  }
  
  if (nrow(d0) == 0 || nrow(m0) == 0 || nrow(r0) == 0) {
    plot.new()
    title(main = paste("Sample", sample_id_value, "- no data"))
    return(invisible(NULL))
  }
  
  y_all <- c(d0$evi2, m0$evi2)
  y_lim <- range(y_all, na.rm = TRUE)
  
  if (!all(is.finite(y_lim))) {
    y_lim <- c(0, 1)
  }
  
  y_pad <- diff(y_lim) * 0.1
  if (!is.finite(y_pad) || y_pad == 0) y_pad <- 0.05
  y_lim <- c(y_lim[1] - y_pad, y_lim[2] + y_pad)
  
  old_par <- par(no.readonly = TRUE)
  on.exit(par(old_par))
  
  layout(matrix(c(1, 2), nrow = 2), heights = c(0.32, 0.68))
  
  ############################################################
  # Top information panel
  ############################################################
  
  par(mar = c(0.2, 0.5, 0.5, 0.5))
  plot.new()
  plot.window(xlim = c(0, 1), ylim = c(0, 1))
  
  title_line <- paste0(
    "Sample ", sample_id_value,
    " | ", r0$change_binary[1],
    " | reason: ", fmt_chr(r0$reason[1])
  )
  
  text(
    x = 0.02,
    y = 0.98,
    labels = title_line,
    adj = c(0, 1),
    font = 2,
    cex = 1.15
  )
  
  info_lines <- c(
    paste0(
      "Break=", fmt_date(r0$break_month[1]),
      " | BIC=", fmt_num(r0$bic_gain[1], 2), " >= ", min_bic_gain,
      "? ", pass_fail(r0$bic_gain[1], min_bic_gain),
      " | model_improved=", bool_text(r0$model_improved[1])
    ),
    
    paste0(
      "Level=", fmt_num(r0$level_shift[1]), " >= ", level_shift_thr,
      "? ", pass_fail(r0$level_shift[1], level_shift_thr),
      " | Shape=", fmt_num(r0$shape_rmse[1]), " >= ", shape_rmse_thr,
      "? ", pass_fail(r0$shape_rmse[1], shape_rmse_thr),
      " | level_shape=", bool_text(r0$changed_by_level_shape[1])
    ),
    
    paste0(
      "Persistent level=", fmt_num(r0$level_shift[1]), " >= ", moderate_level_shift_thr,
      "? ", pass_fail(r0$level_shift[1], moderate_level_shift_thr),
      " | persistent=", bool_text(r0$changed_by_persistent_level[1]),
      " | strong level=", fmt_num(r0$level_shift[1]), " >= ", strong_level_shift_thr,
      "? ", pass_fail(r0$level_shift[1], strong_level_shift_thr),
      " | strong=", bool_text(r0$changed_by_strong_level[1])
    ),
    
    paste0(
      "Amp before=", fmt_num(r0$amp_before[1]),
      " | after=", fmt_num(r0$amp_after[1]),
      " | absent<=", pattern_absent_amp_thr,
      ", present>=", pattern_present_amp_thr,
      " | pattern=", fmt_chr(r0$pattern_transition_type[1]),
      " | pattern_change=", bool_text(r0$changed_by_pattern[1])
    ),
    
    paste0(
      "QC: valid=", r0$valid_months[1], ">=", min_valid_months,
      "? ", pass_fail(r0$valid_months[1], min_valid_months),
      " | max_gap=", r0$max_missing_gap_months[1], "<=", max_missing_gap_months,
      "? ", pass_fail(r0$max_missing_gap_months[1], max_missing_gap_months, "le"),
      " | ", fmt_chr(r0$qc_reason[1])
    ),
    
    paste0(
      "Rule: changed if BIC passes AND ",
      "[level+shape OR persistent level OR strong level OR pattern appeared/disappeared]"
    )
  )
  
  draw_text_lines(info_lines, y_start = 0.78, line_gap = 0.115, cex = 0.76)
  
  ############################################################
  # Bottom time-series panel
  ############################################################
  
  par(mar = c(4.5, 4.7, 0.4, 1.2))
  
  plot(
    d0$date,
    d0$evi2,
    pch = 16,
    cex = 0.45,
    col = "gray60",
    xlab = "Date",
    ylab = "EVI2",
    ylim = y_lim
  )
  
  lines(
    m0$month,
    m0$evi2,
    lwd = 2,
    col = "blue"
  )
  
  abline(
    v = as.Date(paste0(2018:2026, "-01-01")),
    lty = 3,
    col = "gray80"
  )
  
  if (!is.na(r0$break_month[1])) {
    abline(
      v = as.Date(r0$break_month[1]),
      lty = 2,
      col = "red"
    )
  }
  
  legend(
    "topright",
    legend = c("Raw EVI2", "Monthly median", "Best harmonic break"),
    pch = c(16, NA, NA),
    lty = c(NA, 1, 2),
    col = c("gray60", "blue", "red"),
    bty = "n",
    cex = 0.85
  )
  
  invisible(NULL)
}

############################################################
# 3. Read site random points
############################################################

pts <- st_read(gpkg_path, quiet = TRUE)

if (!("sample_id" %in% names(pts))) {
  if ("point_id" %in% names(pts)) {
    pts$sample_id <- pts$point_id
  } else {
    pts$sample_id <- seq_len(nrow(pts))
  }
}

base_r <- raster(base_img_path)

############################################################
# 4. Ensure cell, chunk_id, row_in_chunk exist
############################################################

if (!all(c("cell", "chunk_id", "row_in_chunk") %in% names(pts))) {
  
  cat("cell/chunk_id/row_in_chunk columns not found. Calculating from geometry...\n")
  
  pts_base <- st_transform(pts, st_crs(projection(base_r)))
  xy <- st_coordinates(pts_base)
  
  pts$cell <- cellFromXY(base_r, xy)
  
  pts <- pts[!is.na(pts$cell), ]
  
  ck <- cell_to_chunk(
    cell = pts$cell,
    ncell_total = ncell(base_r),
    numChunks = numChunks
  )
  
  pts$chunk_id <- ck$chunk_id
  pts$row_in_chunk <- ck$row_in_chunk
}

pts_df <- pts %>%
  sf::st_drop_geometry() %>%
  dplyr::select(sample_id, cell, chunk_id, row_in_chunk)

cat("Number of input points:", nrow(pts_df), "\n")

############################################################
# 5. Extract daily EVI2 time series
############################################################

ts_daily <- extract_timeseries_from_chunks(pts_df)

ts_daily <- ts_daily %>%
  filter(date >= analysis_start, date <= analysis_end) %>%
  mutate(
    evi2 = ifelse(
      is.finite(evi2) &
        evi2 >= evi2_min_allowed &
        evi2 <= evi2_max_allowed,
      evi2,
      NA_real_
    )
  )

ts_daily_clean <- ts_daily %>%
  group_by(sample_id, date) %>%
  summarise(
    blue = safe_median(blue),
    green = safe_median(green),
    red = safe_median(red),
    nir = safe_median(nir),
    ndvi = safe_median(ndvi),
    gndvi = safe_median(gndvi),
    nirv = safe_median(nirv),
    evi2 = safe_median(evi2),
    .groups = "drop"
  ) %>%
  arrange(sample_id, date)

saveRDS(ts_daily_clean, out_ts_rds)

cat("Saved daily time series:", out_ts_rds, "\n")

############################################################
# 6. Monthly median EVI2
############################################################

ts_monthly <- ts_daily_clean %>%
  mutate(month = floor_date(date, "month")) %>%
  group_by(sample_id, month) %>%
  summarise(
    evi2 = safe_median(evi2),
    ndvi = safe_median(ndvi),
    nirv = safe_median(nirv),
    n_obs_valid = sum(!is.na(evi2)),
    .groups = "drop"
  ) %>%
  arrange(sample_id, month)

saveRDS(ts_monthly, out_monthly_rds)

cat("Saved monthly EVI2:", out_monthly_rds, "\n")

############################################################
# 7. Binary classification
############################################################

split_monthly <- split(ts_monthly, ts_monthly$sample_id)

result_list <- lapply(names(split_monthly), function(sid) {
  classify_one_point_harmonic(
    d = split_monthly[[sid]],
    sid = as.integer(sid)
  )
})

result_df <- bind_rows(result_list)

############################################################
# 8. Join results and save
############################################################

pts_result <- pts %>%
  left_join(result_df, by = "sample_id")

pts_changed <- pts_result %>%
  filter(change_binary == "changed")

pts_nochange <- pts_result %>%
  filter(change_binary == "no_change" | is.na(change_binary))

write.csv(st_drop_geometry(pts_result), out_all_csv, row.names = FALSE)
write.csv(st_drop_geometry(pts_changed), out_changed_csv, row.names = FALSE)
write.csv(st_drop_geometry(pts_nochange), out_nochange_csv, row.names = FALSE)

if (file.exists(out_all_gpkg)) file.remove(out_all_gpkg)

st_write(
  pts_result,
  out_all_gpkg,
  layer = paste0(site, "_harmonic_change_status"),
  delete_layer = TRUE,
  quiet = TRUE
)

if (file.exists(out_changed_gpkg)) file.remove(out_changed_gpkg)

st_write(
  pts_changed,
  out_changed_gpkg,
  layer = paste0(site, "_changed_points"),
  delete_layer = TRUE,
  quiet = TRUE
)

if (file.exists(out_nochange_gpkg)) file.remove(out_nochange_gpkg)

st_write(
  pts_nochange,
  out_nochange_gpkg,
  layer = paste0(site, "_stable_points"),
  delete_layer = TRUE,
  quiet = TRUE
)

############################################################
# 9. Plot all changed and stable time series
############################################################

changed_ids <- pts_changed$sample_id
nochange_ids <- pts_nochange$sample_id

if (length(changed_ids) > 0) {
  
  pdf(out_pdf_changed, width = 12.5, height = 7.2)
  
  for (sid in changed_ids) {
    plot_one_point_base(
      sample_id_value = sid,
      ts_daily = ts_daily_clean,
      ts_monthly = ts_monthly,
      result_df = result_df
    )
  }
  
  dev.off()
  
  cat("Saved changed PDF:", out_pdf_changed, "\n")
} else {
  cat("No changed points. Changed PDF was not created.\n")
}

if (length(nochange_ids) > 0) {
  
  pdf(out_pdf_nochange, width = 12.5, height = 7.2)
  
  for (sid in nochange_ids) {
    plot_one_point_base(
      sample_id_value = sid,
      ts_daily = ts_daily_clean,
      ts_monthly = ts_monthly,
      result_df = result_df
    )
  }
  
  dev.off()
  
  cat("Saved stable PDF:", out_pdf_nochange, "\n")
} else {
  cat("No stable points. Stable PDF was not created.\n")
}

############################################################
# 10. Summary
############################################################

cat("\n========== SUMMARY ==========\n")
cat("Site:", site, "\n")
cat("Input points:", nrow(pts_result), "\n")
cat("Changed points:", nrow(pts_changed), "\n")
cat("Stable points:", nrow(pts_nochange), "\n\n")

cat("Binary change summary:\n")
print(table(pts_result$change_binary, useNA = "ifany"))

cat("\nReason summary:\n")
print(table(pts_result$reason, useNA = "ifany"))

cat("\nQC summary:\n")
print(table(pts_result$qc_reason, useNA = "ifany"))

cat("\nPattern transition summary:\n")
print(table(pts_result$pattern_transition_type, useNA = "ifany"))

cat("\nModel improved summary:\n")
print(table(pts_result$model_improved, useNA = "ifany"))

cat("\nChanged by level+shape summary:\n")
print(table(pts_result$changed_by_level_shape, useNA = "ifany"))

cat("\nChanged by persistent level summary:\n")
print(table(pts_result$changed_by_persistent_level, useNA = "ifany"))

cat("\nChanged by strong level summary:\n")
print(table(pts_result$changed_by_strong_level, useNA = "ifany"))

cat("\nChanged by pattern summary:\n")
print(table(pts_result$changed_by_pattern, useNA = "ifany"))

cat("\nSaved all status GPKG:\n")
cat(out_all_gpkg, "\n")

cat("\nSaved changed points GPKG:\n")
cat(out_changed_gpkg, "\n")

cat("\nSaved stable points GPKG:\n")
cat(out_nochange_gpkg, "\n")

cat("\nSaved changed PDF:\n")
cat(out_pdf_changed, "\n")

cat("\nSaved stable PDF:\n")
cat(out_pdf_nochange, "\n")

cat("\nDone.\n")