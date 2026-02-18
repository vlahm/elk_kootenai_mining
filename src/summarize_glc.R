library(tidyverse)
library(sf)
library(terra)
library(exactextractr)

## logging setup ####

log_file <- 'logs/summarize_glc.log'
dir.create(dirname(log_file), recursive = TRUE, showWarnings = FALSE)

log_msg <- function(..., console = TRUE) {
  msg <- paste0(...)
  cat(msg, '\n', file = log_file, append = TRUE)
  if (console) message(msg)
}

cat('=== summarize_glc.R ===\n',
    format(Sys.time()), '\n\n', file = log_file)

## configuration ####

glc_dir    <- 'data/glc'
ws_dir     <- 'data/watersheds_nhdhr/sheds_good'
out_file   <- 'data/glc_watershed_summary.csv'
dir.create(glc_dir, recursive = TRUE, showWarnings = FALSE)

# years are determined by the two time periods (see band mappings below)
# early: 1985, 1990, 1995 (3 bands, 5-year intervals)
# late:  2000-2022 (23 bands, annual)

# GLC_FCS30D tiles needed for HUC 1701 study area (~N46-49, W114-118)
# Tile naming: the number is the UPPER boundary of the 5-degree cell.
#   W115 covers lon -120 to -115; W120 covers lon -125 to -120
#   N50  covers lat  45 to  50;  N55  covers lat  50 to  55
# For lat ~46-49 we need N50. Include N55 in case watersheds extend
# into southern BC.
tile_lons <- c('W115', 'W120')
tile_lats <- c('N50', 'N55')

# Zenodo record 15063683 (v2, March 2025)
# Each zip covers a 5-degree longitude band and contains ALL latitude tiles.
# Zip naming: GLC_FCS30D_19852022maps_W115-W120.zip
# Inside each zip, tiles are split into two time periods:
#   GLC_FCS30D_19851995_5years_W**N**.tif  — 3 bands: 1985, 1990, 1995
#   GLC_FCS30D_20002022_W**N**.tif         — 23 bands: 2000, 2001, ..., 2022
zenodo_base <- 'https://zenodo.org/records/15063683/files'

# band-to-year mappings
# early files (19852000_*_5years) have 3 bands: 1985, 1990, 1995
# late files (20002022_*_Annual) have 23 bands: 2000, 2001, ..., 2022
early_years <- c(1985L, 1990L, 1995L)            # 3 bands
early_all   <- c(1985L, 1990L, 1995L)            # all 3 bands in early file
late_years  <- 2000L:2022L                        # 23 bands

## broad class reclassification table ####

# GLC_FCS30D class values -> broad categories
# See: https://data.casearth.cn/sdo/detail/5fbc7904819aec1ea2dd7061
#
# Broad classes encoded as integers:
#   1 = forest, 2 = shrubland, 3 = grassland, 4 = cropland,
#   5 = wetland, 6 = water, 7 = urban, 8 = barren_snow_ice, NA = nodata

broad_class_labels <- c(
  '1' = 'forest',
  '2' = 'shrubland',
  '3' = 'grassland',
  '4' = 'cropland',
  '5' = 'wetland',
  '6' = 'water',
  '7' = 'urban',
  '8' = 'barren_snow_ice'
)

rcl <- matrix(c(
  10,  4,   # Rainfed cropland
  20,  4,   # Irrigated cropland
  51,  1,   # Open evergreen broadleaf forest
  52,  1,   # Closed evergreen broadleaf forest
  61,  1,   # Open deciduous broadleaf forest
  62,  1,   # Closed deciduous broadleaf forest
  71,  1,   # Open evergreen needleleaf forest
  72,  1,   # Closed evergreen needleleaf forest
  81,  1,   # Open deciduous needleleaf forest
  82,  1,   # Closed deciduous needleleaf forest
  91,  1,   # Open mixed-leaf forest
  92,  1,   # Closed mixed-leaf forest
  120, 2,   # Shrubland (evergreen)
  121, 2,   # Shrubland (deciduous)
  122, 2,   # Shrubland (mixed)
  130, 3,   # Grassland
  140, 3,   # Lichens and mosses (treat as grassland)
  150, 8,   # Sparse vegetation (barren-like)
  152, 2,   # Sparse shrubland
  153, 3,   # Sparse herbaceous
  180, 5,   # Wetland
  190, 7,   # Impervious surface
  200, 8,   # Bare land
  201, 8,   # Consolidated bare land
  202, 8,   # Unconsolidated bare land
  210, 6,   # Water body
  220, 8,   # Permanent ice and snow
  250, NA   # Filled / no data
), ncol = 2, byrow = TRUE)


## download and extract GLC_FCS30D tiles ####

# Each zip covers a 5-degree longitude band.  We need to figure out
# which zip(s) to download based on tile_lons.
# Zip naming groups adjacent 5-degree bands: W115-W120.zip contains
# tiles for W115** and W120** (i.e. the band from W115 to W120).
# But from the Zenodo listing, W115-W120 is one zip.
# Our tile_lons W115 and W120 are actually in the SAME zip.

# Build the set of zip files needed.  The zip name format is:
#   GLC_FCS30D_19852022maps_{lower}-{upper}.zip
# where lower and upper are the two 5-degree boundaries.
# W115-W120 covers tiles whose SW corner is W115 or W120.
# Actually from the listing: W115-W120.zip is one file.
# Let's determine the needed zips from tile_lons.

get_zip_name <- function(tile_lon) {
  # parse the longitude value from e.g. "W115"
  dir <- substr(tile_lon, 1, 1)
  val <- as.integer(substr(tile_lon, 2, nchar(tile_lon)))

  # zips group in pairs: W5-W10, W15-W20, ..., W115-W120
  # the lower bound of the pair is the smaller number (closer to 0)
  # each zip spans 10 degrees (two 5-degree tiles)
  lower <- (val %/% 10) * 10 + 5
  upper <- lower + 5
  # but if val is already the lower of a pair:
  if (val %% 10 == 5) {
    lower <- val
    upper <- val + 5
  } else if (val %% 10 == 0) {
    lower <- val - 5
    upper <- val
  }

  paste0('GLC_FCS30D_19852022maps_', dir, lower, '-', dir, upper, '.zip')
}

needed_zips <- unique(sapply(tile_lons, get_zip_name))
log_msg('Need zip file(s): ', paste(needed_zips, collapse = ', '))

# download zips
for (zf in needed_zips) {
  zip_dest <- file.path(glc_dir, zf)
  if (file.exists(zip_dest)) {
    log_msg('  Already have ', zf)
    next
  }
  url <- paste0(zenodo_base, '/', zf, '?download=1')
  log_msg('  Downloading ', zf, ' from ', url, '...')
  tryCatch({
    download.file(url, zip_dest, mode = 'wb', quiet = FALSE)
    log_msg('  Downloaded ', zf, ' (',
            round(file.size(zip_dest) / 1e9, 2), ' GB)')
  }, error = function(e) {
    log_msg('  FAILED to download ', zf, ': ', conditionMessage(e))
    if (file.exists(zip_dest)) file.remove(zip_dest)
    stop('Could not download ', zf, '. Try manually from:\n  ', url)
  })
}

# extract zips
for (zf in needed_zips) {
  zip_path <- file.path(glc_dir, zf)
  # check if already extracted by looking for expected TIF files
  # list contents of zip to find what's inside
  zip_contents <- unzip(zip_path, list = TRUE)$Name
  tif_contents <- zip_contents[grepl('\\.tif$', zip_contents, ignore.case = TRUE)]
  log_msg('  ', zf, ' contains ', length(tif_contents), ' TIF files')

  # check if first TIF already exists (extracted previously)
  if (length(tif_contents) > 0) {
    first_tif <- file.path(glc_dir, basename(tif_contents[1]))
    if (file.exists(first_tif)) {
      log_msg('  Already extracted (found ', basename(first_tif), ')')
      next
    }
  }

  log_msg('  Extracting ', zf, '...')
  unzip(zip_path, exdir = glc_dir, junkpaths = TRUE)
  log_msg('  Extracted ', length(tif_contents), ' TIF files')
}

# Now find the specific TIF files for our tiles.
# Actual naming convention (v1.1):
#   GLC_FCS30D_19852000_W115N45_5years_V1.1.tif  (4 bands: 1985, 1990, 1995, 2000)
#   GLC_FCS30D_20002022_W115N45_Annual_V1.1.tif  (23 bands: 2000-2022)

# Build a manifest of tile files: period × tile_lon × tile_lat
tile_manifest <- expand_grid(
  period = c('early', 'late'),
  tile_lon = tile_lons,
  tile_lat = tile_lats
) %>%
  mutate(
    tile_name = paste0(tile_lon, tile_lat),
    fname = case_when(
      period == 'early' ~ paste0('GLC_FCS30D_19852000_', tile_name, '_5years_V1.1.tif'),
      period == 'late'  ~ paste0('GLC_FCS30D_20002022_', tile_name, '_Annual_V1.1.tif')
    ),
    filepath = file.path(glc_dir, fname),
    exists = file.exists(filepath)
  )

log_msg('\nTile manifest:')
for (i in seq_len(nrow(tile_manifest))) {
  log_msg('  ', tile_manifest$fname[i], ' — ',
          ifelse(tile_manifest$exists[i], 'FOUND', 'MISSING'))
}

if (!all(tile_manifest$exists)) {
  # try to find files by glob pattern as fallback
  for (i in which(!tile_manifest$exists)) {
    tn <- tile_manifest$tile_name[i]
    per <- tile_manifest$period[i]
    pattern <- if (per == 'early') {
      paste0('GLC_FCS30D_1985.*', tn, '.*5years.*\\.tif$')
    } else {
      paste0('GLC_FCS30D_2000.*', tn, '.*Annual.*\\.tif$')
    }
    matches <- list.files(glc_dir, pattern = pattern, full.names = TRUE,
                          ignore.case = TRUE)
    if (length(matches) > 0) {
      log_msg('  Found alternate: ', basename(matches[1]))
      tile_manifest$filepath[i] <- matches[1]
      tile_manifest$exists[i] <- TRUE
    }
  }

  # if still missing, list what's actually in glc_dir
  if (!all(tile_manifest$exists)) {
    all_tifs <- list.files(glc_dir, pattern = '\\.tif$', ignore.case = TRUE)
    log_msg('All TIF files in ', glc_dir, ':')
    for (tf in all_tifs) log_msg('  ', tf)
    stop('Missing required GLC tiles. See log for details.')
  }
}

# verify band counts
for (i in seq_len(nrow(tile_manifest))) {
  r <- rast(tile_manifest$filepath[i])
  nb <- nlyr(r)
  expected <- if (tile_manifest$period[i] == 'early') length(early_all) else length(late_years)
  log_msg('  ', basename(tile_manifest$filepath[i]), ': ', nb, ' bands',
          ifelse(nb == expected, ' (OK)', paste0(' (EXPECTED ', expected, '!)')))
  if (nb != expected) {
    warning('Unexpected band count in ', tile_manifest$filepath[i],
            ': got ', nb, ', expected ', expected)
  }
}


## load watersheds ####

ws_files <- list.files(ws_dir, pattern = '\\.gpkg$', full.names = TRUE)
log_msg('Found ', length(ws_files), ' watershed files in ', ws_dir)

if (length(ws_files) == 0) {
  stop('No watershed gpkg files found in ', ws_dir)
}

# read all watersheds into one sf object for batch extraction
ws_list <- lapply(ws_files, function(f) {
  tryCatch(st_read(f, quiet = TRUE), error = function(e) NULL)
})
ws_list <- ws_list[!sapply(ws_list, is.null)]
watersheds <- bind_rows(ws_list)
log_msg('Loaded ', nrow(watersheds), ' watersheds')

# match CRS to GLC raster (EPSG:4326) to avoid on-the-fly transformation
watersheds <- st_transform(watersheds, 4326)
watersheds <- st_make_valid(watersheds)


## extract land cover fractions ####

log_msg('\n--- Extracting land cover fractions ---')

# helper: for a given year, determine which period and band index
# early file has 4 bands (1985,1990,1995,2000) but we only use bands 1-3
# for 1985/1990/1995; year 2000+ comes from the late (annual) file
get_year_info <- function(yr) {
  if (yr %in% early_years) {
    list(period = 'early', band = which(early_all == yr))
  } else if (yr %in% late_years) {
    list(period = 'late', band = which(late_years == yr))
  } else {
    NULL
  }
}

# all unique years to process (early_years + late_years, no duplicates)
all_years <- sort(unique(c(early_years, late_years)))
log_msg('Processing ', length(all_years), ' years x ',
        nrow(watersheds), ' watersheds')

all_results <- list()
t_start <- proc.time()

# for (yr in all_years[-1]) {
for (yr in all_years) {

  yr_info <- get_year_info(yr)
  if (is.null(yr_info)) next

  # get the tile files for this period
  yr_tiles <- tile_manifest %>%
    filter(period == yr_info$period)

  # load the correct band from each tile and mosaic
  rasters <- lapply(yr_tiles$filepath, function(f) {
    rast(f, lyrs = yr_info$band)
  })

  if (length(rasters) > 1) {
    yr_rast <- do.call(merge, rasters)
  } else {
    yr_rast <- rasters[[1]]
  }

  # reclassify to broad classes
  yr_rcl <- classify(yr_rast, rcl, others = NA)

  # extract fractional coverage of each broad class per watershed
  fracs <- exact_extract(yr_rcl, watersheds, function(values, coverage_fractions) {
    vals <- values[!is.na(values)]
    cov  <- coverage_fractions[!is.na(values)]
    total_weight <- sum(cov)
    if (total_weight == 0) {
      row <- as.data.frame(as.list(setNames(rep(NA_real_, length(broad_class_labels)),
                                            broad_class_labels)))
      return(row)
    }

    class_fracs <- vapply(seq_along(broad_class_labels), function(ci) {
      sum(cov[vals == ci]) / total_weight
    }, numeric(1))
    row <- as.data.frame(as.list(setNames(class_fracs, broad_class_labels)))
    return(row)
  })

  # exact_extract with a fun returning data.frames returns a list of data.frames
  frac_df <- bind_rows(fracs)
  frac_df$comid <- watersheds$comid
  frac_df$year  <- yr

  all_results[[as.character(yr)]] <- frac_df

  elapsed <- (proc.time() - t_start)[['elapsed']]
  yr_idx <- which(all_years == yr)
  log_msg('  Year ', yr, ' complete (', yr_idx, '/', length(all_years),
          ', ', round(elapsed / yr_idx, 1), 's avg per year)')

  # free memory
  rm(yr_rast, yr_rcl, rasters, fracs, frac_df)
  gc()
}

## combine and write results ####

results <- bind_rows(all_results) %>%
  select(comid, year, everything()) %>%
  arrange(comid, year)

log_msg('\nFinal results: ', nrow(results), ' rows (', n_distinct(results$comid),
        ' watersheds x ', n_distinct(results$year), ' years)')

# summary check
log_msg('Column means across all watershed-years:')
for (cl in names(broad_class_labels)) {
  nm <- broad_class_labels[cl]
  log_msg('  ', nm, ': ', round(mean(results[[nm]], na.rm = TRUE) * 100, 1), '%')
}

write_csv(results, out_file)
log_msg('Wrote results to ', out_file)
log_msg('Done.')
