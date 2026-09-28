library(sf)
library(rjson)

############################################################
# Parameters
############################################################

params <- fromJSON(file = "/home/hjlee/project/mangrove/code/PLCM_Parameters.json")

outDir <- params$setup$outDir

in_gpkg <- file.path(outDir, "training_random_points.gpkg")

site_names <- c(
  # "sanur",
  # "gilimanuk",
  # "jembrana",
  # "lembongan"
  # "hungyen",
  # "taklong",
  "aceh"
)

############################################################
# Read input GPKG
############################################################

layers <- st_layers(in_gpkg)
print(layers)

# 첫 번째 layer 읽기
pts <- st_read(
  in_gpkg,
  layer = layers$name[1],
  quiet = TRUE
)

cat("Input points:", nrow(pts), "\n")
cat("Sites in file:\n")
print(table(pts$site, useNA = "ifany"))

############################################################
# Save each site as site_random.gpkg
############################################################

for (site in site_names) {
  
  pts_site <- pts[pts$site == site, ]
  
  cat("\nSite:", site, "\n")
  cat("Number of points:", nrow(pts_site), "\n")
  
  if (nrow(pts_site) == 0) {
    warning("No points found for site: ", site)
    next
  }
  
  site_dir <- file.path(outDir, site)
  
  if (!dir.exists(site_dir)) {
    dir.create(site_dir, recursive = TRUE)
  }
  
  out_site_gpkg <- file.path(site_dir, paste0(site, "_random.gpkg"))
  
  if (file.exists(out_site_gpkg)) {
    file.remove(out_site_gpkg)
  }
  
  st_write(
    pts_site,
    out_site_gpkg,
    layer = paste0(site, "_random"),
    delete_layer = TRUE,
    quiet = TRUE
  )
  
  cat("Saved:", out_site_gpkg, "\n")
}

cat("\nDone.\n")