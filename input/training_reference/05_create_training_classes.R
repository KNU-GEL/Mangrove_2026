library(sf)
library(dplyr)
library(readxl)
library(stringr)

############################################################
# Paths
############################################################

training_dir <- "/home/hjlee/project/mangrove/training_data"
out_dir <- "/home/hjlee/project/mangrove/shp"

in_gpkg <- file.path(training_dir, "training_points.gpkg")
in_xlsx <- file.path(training_dir, "training_points_labeling_all.xlsx")

if (!dir.exists(out_dir)) {
  dir.create(out_dir, recursive = TRUE)
}

############################################################
# Label classes
############################################################

label_classes <- c(
  "mangrove",
  "terrestrial",
  "crop",
  "built-up",
  "deep water",
  "shallow water"
)

out_names <- c(
  "mangrove"      = "mangrove.gpkg",
  "terrestrial"   = "terrestrial.gpkg",
  "crop"          = "crop.gpkg",
  "built-up"      = "built-up.gpkg",
  "deep water"    = "deep_water.gpkg",
  "shallow water" = "shallow_water.gpkg"
)

layer_names <- c(
  "mangrove"      = "mangrove",
  "terrestrial"   = "terrestrial",
  "crop"          = "crop",
  "built-up"      = "built_up",
  "deep water"    = "deep_water",
  "shallow water" = "shallow_water"
)

############################################################
# Read GPKG
############################################################

layers <- st_layers(in_gpkg)

pts <- st_read(
  in_gpkg,
  layer = layers$name[1],
  quiet = TRUE
)

if (!("final_id" %in% names(pts))) {
  stop("training_points.gpkg에 final_id 컬럼이 없습니다.")
}

############################################################
# Read Excel labels
############################################################

labels <- read_excel(in_xlsx)

# 엑셀 컬럼명 정리
names(labels) <- names(labels) %>%
  str_trim() %>%
  str_to_lower()

if (!all(c("id", "label") %in% names(labels))) {
  stop(
    "엑셀에 id 또는 label 컬럼이 없습니다. 현재 컬럼: ",
    paste(names(labels), collapse = ", ")
  )
}

# Excel의 label을 label_excel로 이름 바꿔서 GPKG의 기존 label과 충돌 방지
labels_clean <- labels %>%
  transmute(
    final_id = as.integer(id),
    label_excel = str_trim(as.character(label)),
    label_excel = str_to_lower(label_excel),
    label_excel = na_if(label_excel, "")
  )

############################################################
# Join Excel labels to GPKG points
############################################################

pts_labeled <- pts %>%
  left_join(labels_clean, by = "final_id")

############################################################
# Check labels
############################################################

cat("\nLabel summary:\n")
print(table(pts_labeled$label_excel, useNA = "ifany"))

unknown_labels <- pts_labeled %>%
  st_drop_geometry() %>%
  filter(!is.na(label_excel), !(label_excel %in% label_classes)) %>%
  distinct(label_excel) %>%
  pull(label_excel)

if (length(unknown_labels) > 0) {
  warning(
    "허용된 6개 라벨 외의 값이 있습니다: ",
    paste(unknown_labels, collapse = ", ")
  )
}

############################################################
# Split and save
############################################################

for (lab in label_classes) {
  
  pts_class <- pts_labeled %>%
    filter(label_excel == lab)
  
  out_gpkg <- file.path(out_dir, out_names[[lab]])
  layer_nm <- layer_names[[lab]]
  
  cat("\nLabel:", lab, "\n")
  cat("Number of points:", nrow(pts_class), "\n")
  
  if (nrow(pts_class) == 0) {
    cat("No points. Skipping:", out_gpkg, "\n")
    next
  }
  
  # 여기서 기존 모든 속성 컬럼 제거
  # 최종 GPKG 속성 테이블에는 id 하나만 남고, 값은 전부 1
  pts_class_out <- pts_class %>%
    transmute(id = 1L)
  
  if (file.exists(out_gpkg)) {
    file.remove(out_gpkg)
  }
  
  st_write(
    pts_class_out,
    out_gpkg,
    layer = layer_nm,
    delete_layer = TRUE,
    quiet = TRUE
  )
  
  cat("Saved:", out_gpkg, "\n")
}

cat("\nDone.\n")
cat("Output directory:\n")
cat(out_dir, "\n")