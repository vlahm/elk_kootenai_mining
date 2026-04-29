library(tidyverse)
library(sf)
library(terra)
library(elevatr)

## logging setup ####

log_file <- 'logs/watershed_attributes.log'
dir.create(dirname(log_file), recursive = TRUE, showWarnings = FALSE)

log_msg <- function(..., console = TRUE) {

  msg <- paste0(...)
  cat(msg, '\n', file = log_file, append = TRUE)
  if (console) message(msg)
}

cat('=== watershed_attributes.R ===\n',
    format(Sys.time()), '\n\n', file = log_file)

## setup ####

ws_dir   <- 'data/watersheds_nhdhr'
out_file <- file.path(ws_dir, 'watershed_attributes.csv')
dem_dir  <- 'data/dem_cache'
dir.create(dem_dir, recursive = TRUE, showWarnings = FALSE)

# DEM zoom level for elevatr: 9 ≈ ~30m resolution, good balance of
# speed and accuracy.  Increase to 10–11 for finer resolution if needed.
dem_zoom <- 10

## find watershed files ####

gpkg_files <- list.files(ws_dir, pattern = 'comid_.*\\.gpkg$',
                         full.names = TRUE)
log_msg('Found ', length(gpkg_files), ' watershed files in ', ws_dir)

if (length(gpkg_files) == 0) stop('No watershed gpkg files found in ', ws_dir)

## compute attributes ####

results <- vector('list', length(gpkg_files))

for (i in seq_along(gpkg_files)) {

  f <- gpkg_files[i]
  comid <- str_extract(basename(f), 'comid_([0-9e+.]+)', group = 1)

  if (i %% 10 == 0 || i == 1) {
    log_msg('Processing ', i, ' / ', length(gpkg_files),
            ' (COMID ', comid, ')')
  }

  ws <- tryCatch(st_read(f, quiet = TRUE), error = function(e) NULL)
  if (is.null(ws) || nrow(ws) == 0) {
    log_msg('  WARNING: could not read ', basename(f), ' — skipping')
    results[[i]] <- tibble(comid = comid, area_km2 = NA_real_,
                           mean_slope_deg = NA_real_, status = 'read_error')
    next
  }

  # ── 1. Area ──────────────────────────────────────────────────────────────
  area_m2 <- as.numeric(st_area(st_transform(ws, 5070)))
  area_km2 <- area_m2 / 1e6


  # ── 2. Mean slope from DEM ───────────────────────────────────────────────
  mean_slope <- tryCatch({

    # get_elev_raster returns a RasterLayer; convert to SpatRaster
    dem_rl <- get_elev_raster(
      locations = st_transform(ws, 4326),
      z = dem_zoom,
      clip = 'bbox',
      override_size_check = TRUE
    )
    dem <- rast(dem_rl)

    # compute slope in degrees
    slp <- terrain(dem, v = 'slope', unit = 'degrees')

    # mask to watershed boundary
    ws_vect <- vect(st_transform(ws, crs(dem)))
    slp_masked <- mask(crop(slp, ws_vect), ws_vect)

    # mean slope (ignore NAs at edges)
    as.numeric(global(slp_masked, fun = 'mean', na.rm = TRUE))

  }, error = function(e) {
    log_msg('  WARNING: slope computation failed for COMID ', comid,
            ': ', conditionMessage(e))
    NA_real_
  })

  results[[i]] <- tibble(
    comid          = comid,
    area_km2       = round(area_km2, 4),
    mean_slope_deg = round(mean_slope, 2),
    status         = 'success'
  )
}

## write results ####

attrs <- bind_rows(results)

n_ok   <- sum(attrs$status == 'success', na.rm = TRUE)
n_fail <- sum(attrs$status != 'success', na.rm = TRUE)
log_msg('\nDone. ', n_ok, ' succeeded, ', n_fail, ' failed.')
log_msg('Wrote: ', out_file)
log_msg('\nArea summary (km²):')
log_msg('  Min:    ', round(min(attrs$area_km2, na.rm = TRUE), 2))
log_msg('  Median: ', round(median(attrs$area_km2, na.rm = TRUE), 2))
log_msg('  Max:    ', round(max(attrs$area_km2, na.rm = TRUE), 2))
log_msg('Mean slope summary (degrees):')
log_msg('  Min:    ', round(min(attrs$mean_slope_deg, na.rm = TRUE), 2))
log_msg('  Median: ', round(median(attrs$mean_slope_deg, na.rm = TRUE), 2))
log_msg('  Max:    ', round(max(attrs$mean_slope_deg, na.rm = TRUE), 2))

if(all(attrs$status == 'success')) attrs$status <- NULL
write_csv(attrs, out_file)
