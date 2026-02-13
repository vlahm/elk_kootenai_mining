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


## link sites to NHDPlus COMIDs ####

sites_sf <- st_as_sf(sites, coords = c('lon', 'lat'), crs = 4326)

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
                   'Hydroseq', 'DnHydroseq',
                   'UpHydroseq', 'StartFlag',
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
log_msg('Running get_flowline_index (search_radius = 500m)...')
fi <- tryCatch(
  get_flowline_index(
    flines,
    sites_sf,
    search_radius = units::set_units(50, 'm')
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
# saveRDS(fi, 'data/flowline_indices.rds')
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
  as.numeric(st_distance(sites_with_comid[i, ], fline)[1])
})

sites_with_comid$snap_dist_m <- snap_dist_m

log_msg('\n--- Snap-distance summary (meters) ---')
log_msg('  Min:    ', round(min(snap_dist_m, na.rm = TRUE), 1))
log_msg('  Median: ', round(median(snap_dist_m, na.rm = TRUE), 1))
log_msg('  Mean:   ', round(mean(snap_dist_m, na.rm = TRUE), 1))
log_msg('  Max:    ', round(max(snap_dist_m, na.rm = TRUE), 1))

snap_threshold_m <- 20
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

sites_valid <- filter(sites, !is.na(comid))

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
# HR flowlines have a ToNHDPID or similar downstream pointer
fl_net <- flines %>%
  st_drop_geometry()

# determine column names for network traversal
from_col <- if ('NHDPlusID' %in% names(fl_net)) 'NHDPlusID' else 'COMID'
to_col <- if ('DnHydroSeq' %in% names(fl_net)) {
  # use hydrosequence-based navigation
  NULL
} else if ('ToNHDPID' %in% names(fl_net)) {
  'ToNHDPID'
} else {
  NULL
}

log_msg('Network columns: from=', from_col,
        ', to=', if (is.null(to_col)) 'hydroseq-based' else to_col)

# function to get all upstream COMIDs (including the target)
get_upstream_comids <- function(target_comid, fl_net, from_col) {
  # use nhdplusTools built-in network navigation
  ut <- tryCatch(
    get_UT(fl_net, target_comid),
    error = function(e) target_comid
  )
  return(ut)
}

# function to delineate a watershed by dissolving upstream catchments
delineate_hr_watershed <- function(target_comid, fl_net, catchments,
                                   from_col, catch_id_col) {
  upstream_ids <- get_upstream_comids(target_comid, fl_net, from_col)
  upstream_catches <- catchments[catchments[[catch_id_col]] %in% upstream_ids, ]

  if (nrow(upstream_catches) == 0) {
    warning('No catchments found for COMID ', target_comid)
    return(NULL)
  }

  # dissolve all upstream catchments into one polygon
  basin <- upstream_catches %>%
    st_union() %>%
    st_sf(comid = target_comid, geometry = .)

  return(basin)
}


## single-site test ####

test_comid <- unique_comids[1]
test_site <- sites_valid %>% filter(comid == test_comid) %>% slice(1)

log_msg('--- TEST: site_id = ', test_site$site_id,
        ', COMID = ', test_comid,
        ' (', test_site$lat, ', ', test_site$lon, ') ---')

test_basin <- delineate_hr_watershed(test_comid, fl_net, catchments,
                                     from_col, catch_id_col)

if (!is.null(test_basin) && nrow(test_basin) > 0) {
  test_file <- file.path(ws_dir, paste0('comid_', test_comid, '.shp'))
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

# already_done <- list.files(ws_dir, pattern = 'comid_.*\\.shp$') %>%
#   str_extract('comid_([0-9]+)', group = 1) %>%
#   as.numeric()
# 
# todo_comids <- setdiff(unique_comids, already_done)
# message(length(todo_comids), ' COMIDs remaining to delineate')
# 
# outcomes <- tibble(comid = integer(), status = character(), msg = character())
# 
# for(i in seq_along(todo_comids)){
# 
#   cid <- todo_comids[i]
#   outfile <- file.path(ws_dir, paste0('comid_', cid, '.shp'))
# 
#   if(i %% 50 == 0) message('  progress: ', i, ' / ', length(todo_comids))
# 
#   z <- try({
#     basin <- delineate_hr_watershed(cid, fl_net, catchments,
#                                     from_col, catch_id_col)
#     if(is.null(basin) || nrow(basin) == 0) stop('empty basin')
#     st_write(basin, outfile, delete_dsn = TRUE, quiet = TRUE)
#   }, silent = TRUE)
# 
#   if(inherits(z, 'try-error')){
#     outcomes <- add_row(outcomes, comid = cid, status = 'error',
#                         msg = as.character(z))
#   } else {
#     outcomes <- add_row(outcomes, comid = cid, status = 'success', msg = '')
#   }
# }
# 
# write_csv(outcomes, file.path(ws_dir, 'delineation_outcomes.csv'))
# message('Done. ', sum(outcomes$status == 'success'), ' succeeded, ',
#         sum(outcomes$status == 'error'), ' failed.')


## verify (optional) ####

# dir.create('ws_png', showWarnings = FALSE)
# sheds <- list.files(ws_dir, pattern = '*.shp', full.names = TRUE)
# for(f in sheds){
#   ws <- st_read(f, quiet = TRUE)
#   if(!st_is_longlat(ws)) ws <- st_transform(ws, 4326)
#   pngfile <- file.path('ws_png', paste0(basename(f), '.png'))
#   png(pngfile, width = 900, height = 900)
#   plot(st_geometry(ws), main = basename(f))
#   dev.off()
# }
