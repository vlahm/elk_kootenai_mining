library(tidyverse)
library(sf)
library(mapview)
library(readxl)
library(nhdplusTools)

## logging setup ####

log_file <- 'logs/delineate_watersheds.log'
dir.create(dirname(log_file), recursive = TRUE, showWarnings = FALSE)

# helper: write to log file AND print a short summary to console
log_msg <- function(..., console = TRUE) {
  msg <- paste0(...)
  cat(msg, '\n', file = log_file, append = TRUE)
  if (console) message(msg)
}

# helper: log a data.frame/tibble summary (head + dims) to file only
log_df <- function(df, label = 'data', n = 6) {
  log_msg(sprintf('%s: %d rows x %d cols', label, nrow(df), ncol(df)),
          console = FALSE)
  utils::capture.output(print(utils::head(df, n)), file = log_file,
                        append = TRUE)
}

# start fresh log
cat('=== delineate_watersheds.R ===\n',
    format(Sys.time()), '\n\n', file = log_file)

# collect warnings, print summary at end
options(warn = 1, nwarnings = 500)

## helper: remove interior holes from polygons ####

# Faster than nngeo::st_remove_holes() — works directly on sf geometry.
# For each polygon, keeps only the exterior ring (drops all holes).
remove_holes <- function(x) {
  geom <- st_geometry(x)
  crs <- st_crs(geom)

  no_holes <- lapply(geom, function(g) {
    if (st_is(g, 'POLYGON')) {
      # a POLYGON is a list of rings; first ring is exterior
      st_polygon(list(g[[1]]))
    } else if (st_is(g, 'MULTIPOLYGON')) {
      # each sub-polygon: keep only exterior ring
      # g is a list of polygons; each polygon is a list of rings
      parts <- lapply(g, function(p) list(p[[1]]))
      st_multipolygon(parts)
    } else {
      g
    }
  })

  st_geometry(x) <- st_sfc(no_holes, crs = crs)
  return(x)
}

## setup ####

ws_dir <- 'data/watersheds_nhdhr'
hr_dir <- 'data/nhdplushr'
dir.create(ws_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(hr_dir, recursive = TRUE, showWarnings = FALSE)

# HUC4s covering the Elk-Kootenai study area
# Add more HUC4s here if sites fall outside these
study_huc4s <- c('1701')


## load data ####

# macro sites: unique SiteNumber
macro <- read_csv('data/macro_snapped0.csv', show_col_types = FALSE) %>%
  distinct(site_id = SiteNumber, .keep_all = TRUE) %>%
  rename(lat = Latitude, lon = Longitude) %>% 
  mutate(
    # lat = coalesce(Snap_Y, Latitude),
    # lon = coalesce(Snap_X, Longitude),
    source = 'macro'
  ) %>%
  select(site_id, lat, lon, source)

# chem sites: unique MonitoringLocationIdentifier
conc <- read_xlsx('data/ChemConLoc.xlsx') %>%
  distinct(site_id = MonitoringLocationIdentifier, .keep_all = TRUE) %>%
  rename(lat = LatitudeMeasure, lon = LongitudeMeasure) %>% 
  mutate(source = 'chem') %>%
  select(site_id, lat, lon, source)

# combine all sites
sites <- bind_rows(macro, conc)


## filter sites to Elk-Kootenai study area ####

sites_sf <- st_as_sf(sites, coords = c('lon', 'lat'), crs = 4326)

# use HUC4 boundary to keep only sites within the study basin
# try cached boundary first, then read from GDB
huc4_boundary_file <- file.path(hr_dir, 'huc4_boundary.gpkg')
if (file.exists(huc4_boundary_file)) {
  huc4_boundary <- st_read(huc4_boundary_file, quiet = TRUE)
} else {
  # find the GDB to read WBDHU4 from
  gdb_path <- list.dirs(hr_dir, full.names = TRUE, recursive = TRUE)
  gdb_path <- gdb_path[grepl('\\.gdb$', gdb_path)][1]
  if (!is.na(gdb_path) && 'WBDHU4' %in% st_layers(gdb_path)$name) {
    huc4_boundary <- st_read(gdb_path, layer = 'WBDHU4', quiet = TRUE)
    st_write(huc4_boundary, huc4_boundary_file, quiet = TRUE)
  } else {
    # fallback: use a bounding box for HUC 1701 (approximate)
    log_msg('WARNING: WBDHU4 layer not found — using approximate bbox for 1701')
    huc4_boundary <- st_as_sf(
      st_as_sfc(st_bbox(c(xmin = -117.5, ymin = 46.5,
                           xmax = -114.0, ymax = 49.5),
                         crs = 4326))
    )
  }
}

huc4_boundary <- st_transform(huc4_boundary, st_crs(sites_sf))
in_basin <- st_intersects(sites_sf, st_union(huc4_boundary), sparse = FALSE)[, 1]
n_before <- nrow(sites)
sites <- sites[in_basin, ]
sites_sf <- sites_sf[in_basin, ]
log_msg('Filtered sites to study area: ', nrow(sites), ' of ', n_before,
        ' sites retained')


## link sites to NHDPlus COMIDs ####

# discover_nhdplus_id finds the nearest NHDPlus COMID for each point
# NOTE: this uses NHDPlus Medium Resolution (1:100K). Small headwater
# streams may not exist in this network, causing sites to snap to a
# distant/wrong reach. Check snap distances below and consider NHDPlus
# High-Res (1:24K) for sites with large snap distances.

# --- inspect sites_sf ---
log_msg('sites_sf has ', nrow(sites_sf), ' features')
log_msg('CRS: ', st_crs(sites_sf)$input)
log_msg('First 3 coordinates:', console = FALSE)
log_df(as.data.frame(st_coordinates(sites_sf[1:min(3, nrow(sites_sf)), ])),
       label = 'coords')

# --- download NHDPlus High-Res flowlines & catchments ---
# NHDPlus HR (1:24K) includes small headwater streams and extends
# across the US-Canada border, unlike NHDPlus MR (1:100K).
hr_gpkg <- file.path(hr_dir, 'nhdplushr_merged.gpkg')

if (file.exists(hr_gpkg)) {
  log_msg('Reading cached NHDPlus HR flowlines from ', hr_gpkg)
  flines <- st_read(hr_gpkg, layer = 'NHDFlowline', quiet = TRUE)
  if (nrow(flines) == 0) {
    log_msg('Cached file has 0 flowlines — deleting and re-downloading')
    file.remove(hr_gpkg)
  } else {
    # check if cache is missing critical columns — if so, delete and rebuild
    cached_cols <- names(flines)
    need_cols <- c('REACHCODE', 'ReachCode')
    has_reachcode <- any(need_cols %in% cached_cols)
    if (!has_reachcode) {
      log_msg('WARNING: cached flowlines missing REACHCODE — deleting stale cache')
      rm(flines)
      file.remove(hr_gpkg)
    } else {
      # ensure COMID column exists (needed by get_flowline_index)
      if (!'COMID' %in% cached_cols && 'NHDPlusID' %in% cached_cols) {
        log_msg('Adding COMID column as alias for NHDPlusID (from cache)')
        flines$COMID <- flines$NHDPlusID
      }
      # normalise ReachCode -> REACHCODE if needed
      if (!'REACHCODE' %in% names(flines) && 'ReachCode' %in% names(flines)) {
        log_msg('Renaming ReachCode -> REACHCODE (from cache)')
        flines$REACHCODE <- flines$ReachCode
      }
      # ensure FromMeas/ToMeas exist
      if (!'FromMeas' %in% names(flines)) {
        log_msg('WARNING: cached flowlines missing FromMeas — defaulting to 0')
        flines$FromMeas <- 0
      }
      if (!'ToMeas' %in% names(flines)) {
        log_msg('WARNING: cached flowlines missing ToMeas — defaulting to 100')
        flines$ToMeas <- 100
      }

      # ensure VAA network columns needed by get_UT are present
      needed_vaa_cols <- c('HydroSeq', 'DnHydroSeq', 'LevelPathI')
      missing_net <- setdiff(needed_vaa_cols, names(flines))
      if (length(missing_net) > 0) {
        log_msg('Cached flowlines missing VAA network columns: ',
                paste(missing_net, collapse = ', '),
                ' — joining from source...')
        gdb_path <- list.dirs(hr_dir, full.names = TRUE, recursive = TRUE)
        gdb_path <- gdb_path[grepl('\\.gdb$', gdb_path)][1]
        if (!is.na(gdb_path)) {
          vaa_net <- st_read(gdb_path, layer = 'NHDPlusFlowlineVAA',
                             quiet = TRUE) %>%
            select(any_of(c('NHDPlusID', needed_vaa_cols)))
          # determine join key on flines side
          join_key <- if ('NHDPlusID' %in% names(flines)) 'NHDPlusID' else 'COMID'
          join_by <- setNames('NHDPlusID', join_key)
          flines <- flines %>%
            left_join(st_drop_geometry(vaa_net), by = join_by)
          log_msg('Joined VAA network columns from GDB: ',
                  paste(intersect(needed_vaa_cols, names(flines)),
                        collapse = ', '))
        } else {
          log_msg('WARNING: no GDB found to read VAA network columns')
        }
      }

      log_msg('Cached flowlines: ', nrow(flines), ' rows, ',
              ncol(flines), ' columns')
      log_msg('Columns: ', paste(names(flines), collapse = ', '),
              console = FALSE)
    }
  }
}

if (!file.exists(hr_gpkg)) {
  log_msg('Downloading NHDPlus HR data for HUC4s: ',
          paste(study_huc4s, collapse = ', '))

  # download each HUC4
  hr_files <- character()
  for (huc4 in study_huc4s) {
    # download_nhdplushr puts files under hr_dir
    # look for already-downloaded .gdb directories
    existing_gdb <- list.files(hr_dir,
                               pattern = paste0('NHDPLUS_H_', huc4, '.*\\.gdb$'),
                               full.names = TRUE, recursive = TRUE)
    if (length(existing_gdb) == 0) {
      log_msg('  Downloading HUC4 ', huc4, '...')
      download_nhdplushr(hr_dir, huc4)
    } else {
      log_msg('  HUC4 ', huc4, ' already downloaded: ', existing_gdb[1])
    }

    # search broadly for .gdb directories or .gpkg files related to this HUC4
    gdb_files <- list.files(hr_dir,
                            pattern = paste0('.*', huc4, '.*\\.gdb$'),
                            full.names = TRUE, recursive = TRUE,
                            include.dirs = TRUE)
    # also check for directories ending in .gdb (ESRI GDB is a directory)
    all_dirs <- list.dirs(hr_dir, full.names = TRUE, recursive = TRUE)
    gdb_dirs <- all_dirs[grepl(paste0('.*', huc4, '.*\\.gdb$'), all_dirs)]
    gdb_files <- unique(c(gdb_files, gdb_dirs))

    # also look for .gpkg files
    gpkg_files <- list.files(hr_dir,
                             pattern = paste0('.*', huc4, '.*\\.gpkg$'),
                             full.names = TRUE, recursive = TRUE)

    log_msg('  Found GDB: ', paste(gdb_files, collapse = ', '),
            console = FALSE)
    log_msg('  Found GPKG: ', paste(gpkg_files, collapse = ', '),
            console = FALSE)

    hr_file <- NULL
    if (length(gdb_files) > 0) {
      hr_file <- gdb_files[1]
    } else if (length(gpkg_files) > 0) {
      hr_file <- gpkg_files[1]
    }

    if (is.null(hr_file)) {
      # last resort: look for ANY .gdb or .gpkg anywhere under hr_dir
      any_gdb <- list.dirs(hr_dir, full.names = TRUE, recursive = TRUE)
      any_gdb <- any_gdb[grepl('\\.gdb$', any_gdb)]
      any_gpkg <- list.files(hr_dir, pattern = '\\.gpkg$',
                             full.names = TRUE, recursive = TRUE)
      log_msg('  Any GDB dirs: ', paste(any_gdb, collapse = ', '),
              console = FALSE)
      log_msg('  Any GPKG files: ', paste(any_gpkg, collapse = ', '),
              console = FALSE)
      if (length(any_gdb) > 0) {
        hr_file <- any_gdb[1]
      } else if (length(any_gpkg) > 0) {
        hr_file <- any_gpkg[1]
      } else {
        stop('Could not find downloaded data for HUC4 ', huc4,
             '. Check contents of ', hr_dir)
      }
    }

    log_msg('  Using: ', hr_file)
    hr_files <- c(hr_files, hr_file)
  }

  # list available layers so we pick the right names
  log_msg('Available layers in ', hr_files[1], ':', console = FALSE)
  lyrs <- st_layers(hr_files[1])
  utils::capture.output(print(lyrs), file = log_file, append = TRUE)

  # NHDPlus HR GDBs typically have:
  #   NHDFlowline (flowlines with NHDPlusID)
  #   NHDPlusCatchment (catchment polygons)
  # The flowline layer name varies: NHDFlowline or NHDFlowline_NonNetwork etc.
  fl_layer <- grep('^NHDFlowline$', lyrs$name, value = TRUE)
  if (length(fl_layer) == 0) {
    fl_layer <- grep('NHDFlowline', lyrs$name, value = TRUE)[1]
  }
  catch_layer <- grep('NHDPlusCatchment', lyrs$name, value = TRUE)[1]

  log_msg('Using flowline layer: ', fl_layer)
  log_msg('Using catchment layer: ', catch_layer)

  # read and merge flowlines from all HUC4s
  log_msg('Reading HR flowlines...')
  flines_list <- lapply(hr_files, function(f) {
    st_read(f, layer = fl_layer, quiet = TRUE)
  })
  flines <- bind_rows(flines_list)
  log_msg('Total HR flowlines: ', nrow(flines))

  # read the NHDPlusFlowlineVAA table and join ToMeas/FromMeas
  # (required by get_flowline_index)
  log_msg('Reading NHDPlusFlowlineVAA for measure attributes...')
  vaa_list <- lapply(hr_files, function(f) {
    st_read(f, layer = 'NHDPlusFlowlineVAA', quiet = TRUE)
  })
  vaa <- bind_rows(vaa_list)
  log_msg('VAA rows: ', nrow(vaa))

  # column names to log file only (too long for console)
  log_msg('Flowline columns: ', paste(names(flines), collapse = ', '),
          console = FALSE)
  log_msg('VAA columns: ', paste(names(vaa), collapse = ', '),
          console = FALSE)

  # NHDPlus HR uses NHDPlusID as the join key
  # VAA has ToMeas/FromMeas (or Pathlength, etc.)
  # get_flowline_index needs: COMID, REACHCODE, ToMeas, FromMeas
  wanted_vaa <- c('NHDPlusID', 'ToMeas', 'FromMeas',
                   'HydroSeq', 'DnHydroSeq',
                   'UpHydroSeq', 'LevelPathI', 'StartFlag',
                   'StreamOrde', 'StreamCalc',
                   'TotDASqKm', 'DivDASqKm')
  vaa_cols <- intersect(names(vaa), wanted_vaa)
  missing_vaa <- setdiff(wanted_vaa, names(vaa))
  if (length(missing_vaa) > 0) {
    log_msg('WARNING: VAA missing expected columns: ',
            paste(missing_vaa, collapse = ', '))
    # check for alternate names (NHDPlus HR sometimes uses different cases)
    vaa_lower <- setNames(names(vaa), tolower(names(vaa)))
    for (mc in missing_vaa) {
      alt <- vaa_lower[tolower(mc)]
      if (!is.na(alt)) {
        log_msg('  Found alternate: ', alt, ' for ', mc)
        vaa_cols <- c(vaa_cols, alt)
      }
    }
  }
  vaa_cols <- unique(vaa_cols)
  # drop VAA columns that already exist on flines (except the join key)
  existing_fline_cols <- setdiff(names(flines), 'NHDPlusID')
  vaa_cols_to_join <- setdiff(vaa_cols, existing_fline_cols)
  # always keep the join key
  vaa_cols_to_join <- union('NHDPlusID', vaa_cols_to_join)
  log_msg('Joining VAA columns: ', paste(vaa_cols_to_join, collapse = ', '))
  if (length(setdiff(vaa_cols, vaa_cols_to_join)) > 0) {
    log_msg('  Skipped (already on flowlines): ',
            paste(setdiff(vaa_cols, vaa_cols_to_join), collapse = ', '))
  }
  flines <- flines %>%
    left_join(st_drop_geometry(vaa[, vaa_cols_to_join]), by = 'NHDPlusID')

  # ensure FromMeas/ToMeas exist — fill with 0/100 defaults if missing
  if (!'FromMeas' %in% names(flines)) {
    log_msg('WARNING: FromMeas not found after VAA join — defaulting to 0')
    flines$FromMeas <- 0
  }
  if (!'ToMeas' %in% names(flines)) {
    log_msg('WARNING: ToMeas not found after VAA join — defaulting to 100')
    flines$ToMeas <- 100
  }

  # normalise ReachCode -> REACHCODE if needed
  if (!'REACHCODE' %in% names(flines) && 'ReachCode' %in% names(flines)) {
    log_msg('Renaming ReachCode -> REACHCODE')
    flines$REACHCODE <- flines$ReachCode
    flines$ReachCode <- NULL
  }
  # if still missing, try to get it from VAA
  if (!'REACHCODE' %in% names(flines)) {
    rc_col <- intersect(c('ReachCode', 'REACHCODE', 'Reachcode'), names(vaa))
    if (length(rc_col) > 0) {
      log_msg('Joining REACHCODE from VAA column: ', rc_col[1])
      rc_df <- st_drop_geometry(vaa)[, c('NHDPlusID', rc_col[1])]
      names(rc_df)[2] <- 'REACHCODE'
      flines <- flines %>% left_join(rc_df, by = 'NHDPlusID')
    } else {
      log_msg('WARNING: REACHCODE not found in flowlines or VAA')
    }
  }

  # get_flowline_index expects a COMID column
  if (!'COMID' %in% names(flines) && 'NHDPlusID' %in% names(flines)) {
    log_msg('Adding COMID column as alias for NHDPlusID')
    flines$COMID <- flines$NHDPlusID
  }

  log_msg('Flowlines after VAA join: ', nrow(flines), ' rows, ',
          ncol(flines), ' columns')

  if (nrow(flines) == 0) {
    stop('Read 0 flowlines from HR data. Check the GDB contents.')
  }

  # deduplicate column names before writing (st_write/GPKG is case-insensitive)
  fline_names <- names(flines)
  fline_names_lower <- tolower(fline_names)
  dupes_exact <- fline_names[duplicated(fline_names)]
  dupes_ci <- fline_names[duplicated(fline_names_lower)]
  all_dupes <- unique(c(dupes_exact, dupes_ci))
  if (length(all_dupes) > 0) {
    log_msg('WARNING: removing duplicate columns before caching: ',
            paste(all_dupes, collapse = ', '))
    keep <- !duplicated(fline_names_lower)
    flines <- flines[, keep]
  }

  # cache merged flowlines
  st_write(flines, hr_gpkg, layer = 'NHDFlowline',
           delete_dsn = TRUE, quiet = TRUE)
  log_msg('Cached merged HR flowlines to ', hr_gpkg)

  # read and cache catchments for watershed delineation
  log_msg('Reading HR catchments...')
  catch_list <- lapply(hr_files, function(f) {
    st_read(f, layer = catch_layer, quiet = TRUE)
  })
  catchments <- bind_rows(catch_list)
  log_msg('Total HR catchments: ', nrow(catchments))

  if (nrow(catchments) == 0) {
    stop('Read 0 catchments from HR data. Check the GDB contents.')
  }

  st_write(catchments, hr_gpkg, layer = 'NHDPlusCatchment',
           append = TRUE, quiet = TRUE)
  log_msg('Cached HR catchments to ', hr_gpkg)
}

log_msg('Using ', nrow(flines), ' NHDPlus HR flowlines for COMID matching')

# diagnostic: confirm critical columns exist before get_flowline_index
required_cols <- c('COMID', 'REACHCODE', 'FromMeas', 'ToMeas')
present_cols <- intersect(required_cols, names(flines))
missing_cols <- setdiff(required_cols, names(flines))
log_msg('get_flowline_index required columns present: ',
        paste(present_cols, collapse = ', '))
if (length(missing_cols) > 0) {
  log_msg('ERROR: still missing columns: ',
          paste(missing_cols, collapse = ', '))
  log_msg('All flowline columns: ', paste(names(flines), collapse = ', '),
          console = FALSE)
  stop('Flowlines missing required columns for get_flowline_index: ',
       paste(missing_cols, collapse = ', '),
       '. See ', log_file, ' for details.')
}

# --- snap sites to HR flowlines using get_flowline_index ---

fi_cache <- 'data/flowline_indices.rds'

if (file.exists(fi_cache)) {
  log_msg('Loading cached flowline indices from ', fi_cache)
  fi <- readRDS(fi_cache)
  log_msg('Loaded ', nrow(fi), ' flowline index rows from cache')
} else {
  # pre-process flowlines: match CRS, drop Z, cast to LINESTRING
  flines_idx <- flines %>%
    st_transform(st_crs(sites_sf)) %>%
    st_zm(drop = TRUE, what = 'ZM') %>%
    st_cast('LINESTRING')

  log_msg('Running get_flowline_index (search_radius = 50m)...')
  fi <- tryCatch(
    get_flowline_index(
      flines_idx,
      sites_sf,
      search_radius = units::set_units(500, 'm')
    ),
    error = function(e) {
      log_msg('ERROR in get_flowline_index: ', conditionMessage(e))
      log_msg('Flowline columns at call time: ',
              paste(names(flines), collapse = ', '), console = FALSE)
      log_msg('Flowline nrow: ', nrow(flines), ', CRS: ',
              st_crs(flines)$input, console = FALSE)
      stop('get_flowline_index failed. See ', log_file, ' for details.',
           call. = FALSE)
    }
  )
  rm(flines_idx)  # free memory
  saveRDS(fi, fi_cache)
  log_msg('Saved flowline indices to ', fi_cache)
}
log_msg('get_flowline_index returned ', nrow(fi), ' rows')

# fi has columns: id (row index), COMID, REACHCODE, REACH_meas, offset
# offset is in CRS units; id corresponds to row number in sites_sf
comids <- rep(NA_integer_, nrow(sites))
comids[fi$id] <- fi$COMID
sites$comid <- comids

# some may fail; report
n_found <- sum(!is.na(sites$comid))
log_msg(n_found, ' of ', nrow(sites), ' sites matched to a COMID')

# --- snap-distance diagnostic ---
# use the HR flowlines already in memory (no re-download needed)
matched_comids <- unique(na.omit(sites$comid))
flines_sf <- flines %>%
  filter(NHDPlusID %in% matched_comids | COMID %in% matched_comids) %>%
  st_transform(st_crs(sites_sf))
log_msg('Matched ', nrow(flines_sf), ' HR flowlines for snap-distance check')

sites_with_comid <- sites %>%
  filter(!is.na(comid)) %>%
  st_as_sf(coords = c('lon', 'lat'), crs = 4326)

# HR flowlines may use NHDPlusID or COMID — find the right column
fline_id_col <- if ('NHDPlusID' %in% names(flines_sf)) 'NHDPlusID' else 'COMID'

snap_dist_m <- map_dbl(seq_len(nrow(sites_with_comid)), function(i) {
  fline <- flines_sf[flines_sf[[fline_id_col]] == sites_with_comid$comid[i], ]
  if (nrow(fline) == 0) return(NA_real_)
  suppressMessages(as.numeric(st_distance(sites_with_comid[i, ], fline)[1]))
})

sites_with_comid$snap_dist_m <- snap_dist_m

log_msg('\n--- Snap-distance summary (meters) ---')
log_msg('  Min:    ', round(min(snap_dist_m, na.rm = TRUE), 1))
log_msg('  Median: ', round(median(snap_dist_m, na.rm = TRUE), 1))
log_msg('  Mean:   ', round(mean(snap_dist_m, na.rm = TRUE), 1))
log_msg('  Max:    ', round(max(snap_dist_m, na.rm = TRUE), 1))

snap_threshold_m <- 500
n_suspect <- sum(snap_dist_m > snap_threshold_m, na.rm = TRUE)
log_msg('  Sites > ', snap_threshold_m, 'm from nearest flowline: ', n_suspect)
if (n_suspect > 0) {
  log_msg('  ** Review suspect sites in log file **')
  suspect_sites <- sites_with_comid %>%
    st_drop_geometry() %>%
    filter(snap_dist_m > snap_threshold_m) %>%
    arrange(desc(snap_dist_m))
  log_df(suspect_sites, label = 'suspect_sites', n = 20)
}

write_csv(st_drop_geometry(sites_with_comid),
          file.path(ws_dir, 'site_snap_distances.csv'))
log_msg('Wrote snap-distance diagnostics: ',
        file.path(ws_dir, 'site_snap_distances.csv'))

sites_valid <- filter(sites, !is.na(comid), !site_id %in% suspect_sites$site_id)

# build lookup: site_id -> comid (many sites may share a reach)
site_comid_lookup <- sites_valid %>%
  select(site_id, comid, source)

write_csv(site_comid_lookup, file.path(ws_dir, 'site_comid_lookup.csv'))
log_msg('Wrote site-to-COMID lookup: ', file.path(ws_dir, 'site_comid_lookup.csv'))

# unique reaches to delineate
unique_comids <- sort(unique(sites_valid$comid))
log_msg(length(unique_comids), ' unique COMIDs to delineate')


## load HR catchments for local delineation ####

if (!exists('catchments')) {
  log_msg('Reading HR catchments from ', hr_gpkg)
  catchments <- st_read(hr_gpkg, layer = 'NHDPlusCatchment', quiet = TRUE)
}
log_msg('Loaded ', nrow(catchments), ' HR catchments')

# identify the catchment ID column
catch_id_col <- if ('NHDPlusID' %in% names(catchments)) 'NHDPlusID' else 'FEATUREID'

# build a simple network lookup from the flowlines for upstream traversal
# get_UT (hydroloom) requires: id, levelpath, topo_sort, dn_topo_sort
# These map from NHDPlus HR columns as follows:
#   id          <- COMID (= NHDPlusID)
#   topo_sort   <- HydroSeq
#   dn_topo_sort <- DnHydroSeq
#   levelpath   <- LevelPathI

# first, read LevelPathI from VAA if not already on flines
if (!'LevelPathI' %in% names(flines)) {
  # try to get it from the cached gpkg or GDB
  log_msg('LevelPathI not on flines — reading from source...')
  if (file.exists(hr_gpkg)) {
    # check if VAA was cached as a layer
    gpkg_lyrs <- st_layers(hr_gpkg)$name
    if ('NHDPlusFlowlineVAA' %in% gpkg_lyrs) {
      vaa_lp <- st_read(hr_gpkg, layer = 'NHDPlusFlowlineVAA', quiet = TRUE) %>%
        select(any_of(c('NHDPlusID', 'LevelPathI')))
    } else {
      # read from original GDB
      gdb_path <- list.dirs(hr_dir, full.names = TRUE, recursive = TRUE)
      gdb_path <- gdb_path[grepl('\\.gdb$', gdb_path)][1]
      vaa_lp <- st_read(gdb_path, layer = 'NHDPlusFlowlineVAA', quiet = TRUE) %>%
        select(any_of(c('NHDPlusID', 'LevelPathI')))
    }
    if ('LevelPathI' %in% names(vaa_lp)) {
      flines <- flines %>%
        left_join(st_drop_geometry(vaa_lp), by = c('COMID' = 'NHDPlusID'))
      log_msg('Joined LevelPathI: ', sum(!is.na(flines$LevelPathI)),
              ' of ', nrow(flines), ' flowlines matched')
    } else {
      log_msg('WARNING: LevelPathI not found in VAA table')
    }
  }
}

# Build a fast adjacency list for upstream traversal.
# Each flowline has a HydroSeq (its own topo sort) and DnHydroSeq (its
# downstream neighbour's topo sort).  We invert this to get a mapping
# from each HydroSeq → all COMIDs whose DnHydroSeq equals it, i.e. its
# upstream children.  Then we also need a COMID → HydroSeq lookup so we
# can start the walk from a target COMID.

log_msg('Building upstream adjacency list...')

fl_net <- flines %>%
  st_drop_geometry() %>%
  select(COMID,
         any_of(c('HydroSeq', 'Hydroseq', 'HydroSeq.x',
                   'DnHydroSeq', 'DnHydroseq', 'DnHydroSeq.x'))) %>%
  rename_with(~ 'HydroSeq',   .cols = any_of(c('Hydroseq', 'HydroSeq.x'))) %>%
  rename_with(~ 'DnHydroSeq', .cols = any_of(c('DnHydroseq', 'DnHydroSeq.x'))) %>%
  filter(!is.na(HydroSeq), !is.na(DnHydroSeq))

# COMID → HydroSeq lookup (named vector, character keys for env lookup)
comid_to_hs <- setNames(fl_net$HydroSeq, as.character(fl_net$COMID))

# HydroSeq → COMID lookup
hs_to_comid <- setNames(fl_net$COMID, as.character(fl_net$HydroSeq))

# Build children map: DnHydroSeq → vector of upstream COMIDs
# (i.e. all COMIDs that drain INTO a given HydroSeq)
children_env <- new.env(hash = TRUE, parent = emptyenv(),
                        size = as.integer(nrow(fl_net) * 1.3))
for (i in seq_len(nrow(fl_net))) {
  key <- as.character(fl_net$DnHydroSeq[i])
  cid <- fl_net$COMID[i]
  if (exists(key, envir = children_env, inherits = FALSE)) {
    assign(key, c(get(key, envir = children_env), cid),
           envir = children_env)
  } else {
    assign(key, cid, envir = children_env)
  }
}

log_msg('Adjacency list built: ', length(comid_to_hs), ' flowlines, ',
        length(ls(children_env)), ' unique downstream nodes')

# --- Pre-compute ALL upstream COMIDs in one topological sweep ---
# Process flowlines in ascending HydroSeq order (headwaters first).
# Each node's upstream set = itself + union of its children's upstream sets.
# This turns per-site BFS into a single O(n) pass + O(1) lookup.

log_msg('Pre-computing upstream COMID sets for all ',
        nrow(fl_net), ' flowlines (topological sweep)...')
t_precomp <- proc.time()

# In NHDPlus HR, headwaters have the HIGHEST HydroSeq values.
# Process in descending order so children (upstream) are resolved
# before their downstream parents.
fl_sorted <- fl_net[order(fl_net$HydroSeq, decreasing = TRUE), ]

# pre-allocate list keyed by character COMID
upstream_list <- vector('list', nrow(fl_sorted))
names(upstream_list) <- as.character(fl_sorted$COMID)

# also build a fast index: character COMID → position in upstream_list
ul_idx <- setNames(seq_along(upstream_list), names(upstream_list))

for (i in seq_len(nrow(fl_sorted))) {
  cid <- fl_sorted$COMID[i]
  hs_key <- as.character(fl_sorted$HydroSeq[i])

  # get immediate upstream children (COMIDs whose DnHydroSeq == this HydroSeq)
  kids <- if (exists(hs_key, envir = children_env, inherits = FALSE)) {
    get(hs_key, envir = children_env)
  } else {
    NULL
  }

  if (is.null(kids) || length(kids) == 0) {
    # headwater: upstream set is just itself
    upstream_list[[i]] <- cid
  } else {
    # gather upstream sets of all children + self
    kid_keys <- as.character(kids)
    kid_positions <- ul_idx[kid_keys]
    kid_positions <- kid_positions[!is.na(kid_positions)]

    if (length(kid_positions) == 0) {
      upstream_list[[i]] <- c(cid, kids)
    } else {
      upstream_list[[i]] <- c(cid, unique(unlist(
        upstream_list[kid_positions], use.names = FALSE
      )))
    }
  }
}

elapsed <- (proc.time() - t_precomp)[['elapsed']]
total_entries <- sum(vapply(upstream_list, length, integer(1)))
log_msg('Upstream pre-computation complete in ', round(elapsed, 1), 's')
log_msg('  Total entries across all upstream sets: ',
        format(total_entries, big.mark = ','),
        ' (approx ', round(total_entries * 8 / 1e9, 2), ' GiB)')

# diagnostic: distribution of upstream set sizes
us_lengths <- vapply(upstream_list, length, integer(1))
log_msg('  Upstream set size distribution:')
log_msg('    Min: ', min(us_lengths), '  Median: ', median(us_lengths),
        '  Mean: ', round(mean(us_lengths), 1), '  Max: ', max(us_lengths))
log_msg('    Sets with >10 COMIDs: ', sum(us_lengths > 10),
        '  >100: ', sum(us_lengths > 100),
        '  >1000: ', sum(us_lengths > 1000))

# O(1) lookup wrapper (keeps the same interface for delineate_hr_watershed)
get_upstream_comids <- function(target_comid, comid_to_hs, children_env) {
  key <- as.character(target_comid)
  us <- upstream_list[[key]]
  if (is.null(us)) return(target_comid)
  return(us)
}

# function to delineate a watershed by dissolving ALL upstream catchments
# get_UT walks the full network upstream via hydrosequence, so the result
# is the entire contributing area above the target reach, not just the
# immediate local catchment.
delineate_hr_watershed <- function(target_comid, comid_to_hs, children_env,
                                   catchments, catch_id_col) {
  upstream_ids <- get_upstream_comids(target_comid, comid_to_hs, children_env)

  # upstream_ids come from fl_net$id which holds NHDPlusID values
  upstream_catches <- catchments[catchments[[catch_id_col]] %in% upstream_ids, ]

  if (nrow(upstream_catches) == 0) {
    warning('No catchments found for COMID ', target_comid,
            ' (', length(upstream_ids), ' upstream COMIDs)')
    return(NULL)
  }

  # dissolve all upstream catchments into one polygon, then remove
  # interior holes (DEM artifacts from roads, railroads, etc.)
  # Two-pass approach:
  #   1. remove_holes() strips interior rings (large holes)
  #   2. buffer out then back in to close single-pixel / diagonal-adjacency
  #      gaps that st_union leaves between adjacent catchments
  basin_geom <- upstream_catches %>%
    st_make_valid() %>%
    st_union() %>%
    st_make_valid()

  # close micro-gaps: buffer out 1m, then back in 1m
  # must project to a metric CRS first so buffer distance is in meters
  orig_crs <- st_crs(basin_geom)
  basin_geom <- basin_geom %>%
    st_transform(32611) %>%
    st_buffer(1) %>%
    st_union() %>%
    st_make_valid() %>%
    st_buffer(-1) %>%
    st_make_valid() %>%
    st_transform(orig_crs)

  basin <- st_sf(comid = target_comid,
                 n_upstream = length(upstream_ids),
                 n_catchments = nrow(upstream_catches),
                 geometry = basin_geom) %>%
    remove_holes()

  return(basin)
}


## single-site test ####

# from_col is now 'id' to match fl_net; catch_id_col matches catchments
test_comid <- unique_comids[1]
test_site <- sites_valid %>% filter(comid == test_comid) %>% slice(1)

log_msg('--- TEST: site_id = ', test_site$site_id,
        ', COMID = ', test_comid,
        ' (', test_site$lat, ', ', test_site$lon, ') ---')

test_basin <- delineate_hr_watershed(test_comid, comid_to_hs, children_env,
                                     catchments, catch_id_col)

if (!is.null(test_basin) && nrow(test_basin) > 0) {
  log_msg('  Upstream COMIDs: ', test_basin$n_upstream,
          ', catchments dissolved: ', test_basin$n_catchments)
  test_file <- file.path(ws_dir, paste0('comid_', test_comid, '.gpkg'))
  st_write(test_basin, test_file, delete_dsn = TRUE, quiet = TRUE)
  log_msg('TEST SUCCESS: wrote ', test_file)

  test_pt <- st_as_sf(test_site, coords = c('lon', 'lat'), crs = 4326)
  print(mapview(test_basin, layer.name = 'watershed') +
          mapview(test_pt, col.regions = 'red', layer.name = 'site'))
} else {
  log_msg('TEST FAILED: no basin returned for COMID ', test_comid)
}

log_msg('Test complete. Review the map, then run the full loop below.')
log_msg('Full log written to: ', log_file)


## delineate all unique reaches ####
## (uncomment the block below to run the full loop)

# read previous outcomes to identify failures that should be retried
prev_outcomes_file <- file.path(ws_dir, 'delineation_outcomes.csv')
prev_failures <- character()
if (file.exists(prev_outcomes_file)) {
  prev_outcomes <- read_csv(prev_outcomes_file, show_col_types = FALSE)
  prev_failures <- prev_outcomes %>%
    filter(status == 'error') %>%
    pull(comid) %>%
    as.character()
  log_msg('Found ', length(prev_failures), ' previous failures to retry')
}

already_done <- list.files(ws_dir, pattern = 'comid_.*\\.gpkg$') %>%
  str_extract('comid_([0-9e+]+)', group = 1)

# remove previous failures from already_done so they get retried
already_done <- setdiff(already_done, prev_failures)

todo_comids <- setdiff(unique_comids, already_done)
log_msg(length(todo_comids), ' COMIDs remaining to delineate')

outcomes <- tibble(comid = integer(), status = character(), msg = character())

for(i in seq_along(todo_comids)){

  cid <- todo_comids[i]
  outfile <- file.path(ws_dir, paste0('comid_', cid, '.gpkg'))

  if(i %% 50 == 0) message('  progress: ', i, ' / ', length(todo_comids))

  z <- try({
    basin <- delineate_hr_watershed(cid, comid_to_hs, children_env,
                                    catchments, catch_id_col)
    if(is.null(basin) || nrow(basin) == 0) stop('empty basin')
    st_write(basin, outfile, delete_dsn = TRUE, quiet = TRUE)
  }, silent = TRUE)

  if(inherits(z, 'try-error')){
    outcomes <- add_row(outcomes, comid = cid, status = 'error',
                        msg = as.character(z))
  } else {
    outcomes <- add_row(outcomes, comid = cid, status = 'success', msg = '')
  }
}

write_csv(outcomes, file.path(ws_dir, 'delineation_outcomes.csv'))
message('Done. ', sum(outcomes$status == 'success'), ' succeeded, ',
        sum(outcomes$status == 'error'), ' failed.')


## retroactively fix holes in existing gpkg files ####

# log_msg('\n--- Removing holes from existing watershed gpkg files ---')
# existing_gpkgs <- list.files(ws_dir, pattern = 'comid_.*\\.gpkg$',
#                              full.names = TRUE)
# log_msg('Found ', length(existing_gpkgs), ' existing watershed files to fix')
# 
# n_fixed <- 0
# n_had_holes <- 0
# for (gpkg_f in existing_gpkgs) {
#   ws <- tryCatch(st_read(gpkg_f, quiet = TRUE), error = function(e) NULL)
#   if (is.null(ws) || nrow(ws) == 0) next
# 
#   # count holes before
#   geom <- st_geometry(ws)
#   n_holes_before <- sum(vapply(geom, function(g) {
#     if (st_is(g, 'POLYGON')) {
#       max(0L, length(g) - 1L)
#     } else if (st_is(g, 'MULTIPOLYGON')) {
#       sum(vapply(g, function(p) max(0L, length(p) - 1L), integer(1)))
#     } else {
#       0L
#     }
#   }, integer(1)))
# 
#   # always reprocess: previous fix attempts may have been incomplete
#   # (broken remove_holes or lon/lat buffer)
#   orig_crs <- st_crs(ws)
#   ws_fixed <- ws %>%
#     st_transform(32611) %>%
#     st_buffer(1) %>%
#     st_union() %>%
#     st_make_valid() %>%
#     st_buffer(-1) %>%
#     st_make_valid() %>%
#     st_transform(orig_crs) %>%
#     st_sf(geometry = .) %>%
#     mutate(comid = ws$comid,
#            n_upstream = ws$n_upstream,
#            n_catchments = ws$n_catchments) %>%
#     remove_holes()
# 
#   # count holes after to report
#   geom_after <- st_geometry(ws_fixed)
#   n_holes_after <- sum(vapply(geom_after, function(g) {
#     if (st_is(g, 'POLYGON')) {
#       max(0L, length(g) - 1L)
#     } else if (st_is(g, 'MULTIPOLYGON')) {
#       sum(vapply(g, function(p) max(0L, length(p) - 1L), integer(1)))
#     } else {
#       0L
#     }
#   }, integer(1)))
# 
#   if (n_holes_before > 0) n_had_holes <- n_had_holes + 1
#   if (n_holes_after > 0) {
#     log_msg('  WARNING: ', basename(gpkg_f), ' still has ',
#             n_holes_after, ' holes after fix')
#   }
# 
#   st_write(ws_fixed, gpkg_f, delete_dsn = TRUE, quiet = TRUE)
#   n_fixed <- n_fixed + 1
#   if (n_fixed %% 50 == 0) {
#     log_msg('  Fixed ', n_fixed, ' / ', length(existing_gpkgs), ' files...')
#   }
# }
# log_msg('Hole removal complete: ', n_fixed, ' files rewritten (',
#         n_had_holes, ' originally had holes)')

## verify (optional) ####

dir.create('ws_png', showWarnings = FALSE)
sheds <- list.files(ws_dir, pattern = '*.gpkg', full.names = TRUE)
# lame <- c()
for(f in sheds[1:10]){
  ws <- st_read(f, quiet = TRUE)
  if(!st_is_longlat(ws)) ws <- st_transform(ws, 4326)
  pngfile <- file.path('ws_png', paste0(basename(f), '.png'))
  png(pngfile, width = 900, height = 900)
  tryres <- try(plot(st_geometry(ws), main = basename(f)), silent = TRUE)
  dev.off()
  # if(inherits(tryres, 'try-error')) lame <- c(lame, f)
}

# for(f in lame){
#   file.remove(f)
# }
# for(f in lame){
#   file.copy(from = f, to = paste0('data/nhdhr_bak/', basename(f)))
# }

## investigate questionable snaps ####

questionables <- sites_with_comid %>% 
  filter(snap_dist_m > 400)

redo <- c()
for(i in seq_len(nrow(questionables))){
  
  comid <- questionables$comid[i]
  sitename <- questionables$site_id[i]
  print(comid)
  zshed <- try(st_read(paste0('data/watersheds_nhdhr/comid_', comid, '.gpkg')))
  if(inherits(zshed, 'try-error')) next
  zsite <- filter(sites_sf, site_id == sitename)
  
  # test_pt <- st_as_sf(test_site, coords = c('lon', 'lat'), crs = 4326)
  print(mapview(zshed, layer.name = 'sitename') +
        mapview(zsite, col.regions = 'red', layer.name = 'comid'))
  
  # ws <- st_read(f, quiet = TRUE)
  # mv <- mapview(ws)
  # print(mv)
  message('enter "x" to mark as "redo", or any other key to see next')
  input <- readLines(n = 1)
  if(input == 'x') redo <- c(redo, f)
}


suspect_sites <- bind_rows(
  suspect_sites,
  st_drop_geometry(filter(sites_with_comid, comid %in% redo))
)

write_csv(suspect_sites, 'data/watersheds_nhdhr/questionable_sheds.csv')

for(cc in suspect_sites$comid){
  file.rename(paste0('data/watersheds_nhdhr/comid_', cc, '.gpkg'),
              paste0('data/watersheds_nhdhr/sheds_questionable/comid_', cc, '.gpkg'))
}
