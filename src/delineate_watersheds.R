library(tidyverse)
library(sf)
library(mapview)
library(readxl)
library(nhdplusTools)

## setup ####

ws_dir <- 'data/watersheds_nhd'
dir.create(ws_dir, recursive = TRUE, showWarnings = FALSE)


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
message('sites_sf has ', nrow(sites_sf), ' features')
message('CRS: ', st_crs(sites_sf)$input)
message('First 3 coordinates:')
print(st_coordinates(sites_sf[1:3, ]))
message('')

# --- download NHDPlus flowlines for the study area ---
# discover_nhdplus_id and the NLDI position endpoint are unreliable,
# so we download flowlines locally and use get_flowline_index instead.
nhd_gpkg <- file.path(ws_dir, 'nhdplus_subset.gpkg')

if (file.exists(nhd_gpkg)) {
  message('Reading cached NHDPlus flowlines from ', nhd_gpkg)
  flines <- sf::st_read(nhd_gpkg, layer = 'NHDFlowline_Network', quiet = TRUE)
} else {
  message('Downloading NHDPlus flowlines for study area bounding box...')
  # subset_nhdplus needs COMIDs, so we use get_nhdplus to query by bbox
  # First, get the start COMID from any point using navigate_nldi or
  # download by bounding box via nhdplusTools
  bbox <- st_bbox(sites_sf)

  # Build a small polygon from the bbox for the download
  bbox_poly <- st_as_sfc(bbox)

  # Use get_nhdplus to download flowlines intersecting the bbox
  flines <- tryCatch(
    get_nhdplus(AOI = bbox_poly, realization = 'flowline'),
    error = function(e) {
      message('get_nhdplus failed: ', conditionMessage(e))
      # try alternate: download_nhdplusv2 or direct query
      NULL
    }
  )

  if (is.null(flines) || nrow(flines) == 0) {
    stop('Could not download NHDPlus flowlines. Check network connectivity.')
  }

  message('Downloaded ', nrow(flines), ' flowlines')
  st_write(flines, nhd_gpkg, layer = 'NHDFlowline_Network',
           delete_dsn = TRUE, quiet = TRUE)
  message('Cached flowlines to ', nhd_gpkg)
}

message('Using ', nrow(flines), ' NHDPlus flowlines for COMID matching')

# --- snap sites to flowlines using get_flowline_index ---
message('Running get_flowline_index (search_radius = 500m)...')
fi <- get_flowline_index(flines, sites_sf,
                         search_radius = units::set_units(20, 'm'))
message('get_flowline_index returned ', nrow(fi), ' rows')

# fi has columns: id (row index), COMID, REACHCODE, REACH_meas, offset
# offset is in CRS units; id corresponds to row number in sites_sf
comids <- rep(NA_integer_, nrow(sites))
comids[fi$id] <- as.integer(fi$COMID)
sites$comid <- comids

# some may fail; report
n_found <- sum(!is.na(sites$comid))
message(n_found, ' of ', nrow(sites), ' sites matched to a COMID')

# --- snap-distance diagnostic ---
# download the matched flowlines and compute distance from each site
# to its assigned NHDPlus reach
matched_comids <- unique(na.omit(sites$comid))
matched_flines <- subset_nhdplus(comids = matched_comids,
                                  output_file = tempfile(fileext = ".gpkg"),
                                  nhdplus_data = "download",
                                  flowline_only = TRUE,
                                  return_data = TRUE,
                                  overwrite = TRUE)
flines_sf <- st_transform(matched_flines$NHDFlowline_Network, st_crs(sites_sf))

sites_with_comid <- sites %>%
  filter(!is.na(comid)) %>%
  st_as_sf(coords = c('lon', 'lat'), crs = 4326)

snap_dist_m <- map_dbl(seq_len(nrow(sites_with_comid)), function(i) {
  fline <- filter(flines_sf, comid == sites_with_comid$comid[i])
  if (nrow(fline) == 0) return(NA_real_)
  as.numeric(st_distance(sites_with_comid[i, ], fline)[1])
})

sites_with_comid$snap_dist_m <- snap_dist_m

message('\n--- Snap-distance summary (meters) ---')
message('  Min:    ', round(min(snap_dist_m, na.rm = TRUE), 1))
message('  Median: ', round(median(snap_dist_m, na.rm = TRUE), 1))
message('  Mean:   ', round(mean(snap_dist_m, na.rm = TRUE), 1))
message('  Max:    ', round(max(snap_dist_m, na.rm = TRUE), 1))

snap_threshold_m <- 20
n_suspect <- sum(snap_dist_m > snap_threshold_m, na.rm = TRUE)
message('  Sites > ', snap_threshold_m, 'm from nearest flowline: ', n_suspect)
if (n_suspect > 0) {
  message('  ** These sites likely need NHDPlus High-Res delineation **')
  suspect_sites <- sites_with_comid %>%
    st_drop_geometry() %>%
    filter(snap_dist_m > snap_threshold_m) %>%
    arrange(desc(snap_dist_m))
  print(suspect_sites)
}

write_csv(st_drop_geometry(sites_with_comid),
          file.path(ws_dir, 'site_snap_distances.csv'))
message('Wrote snap-distance diagnostics: ',
        file.path(ws_dir, 'site_snap_distances.csv'))

sites_valid <- filter(sites, !is.na(comid))

# build lookup: site_id -> comid (many sites may share a reach)
site_comid_lookup <- sites_valid %>%
  select(site_id, comid, source)

write_csv(site_comid_lookup, file.path(ws_dir, 'site_comid_lookup.csv'))
message('Wrote site-to-COMID lookup: ', file.path(ws_dir, 'site_comid_lookup.csv'))

# unique reaches to delineate
unique_comids <- sort(unique(sites_valid$comid))
message(length(unique_comids), ' unique COMIDs to delineate')


## single-site test ####

test_comid <- unique_comids[1]
test_site <- sites_valid %>% filter(comid == test_comid) %>% slice(1)

message('--- TEST: site_id = ', test_site$site_id,
        ', COMID = ', test_comid,
        ' (', test_site$lat, ', ', test_site$lon, ') ---')

test_basin <- get_nldi_basin(list(featureSource = 'comid', featureID = test_comid))

if(nrow(test_basin) > 0){
  test_file <- file.path(ws_dir, paste0('comid_', test_comid, '.shp'))
  st_write(test_basin, test_file, delete_dsn = TRUE, quiet = TRUE)
  message('TEST SUCCESS: wrote ', test_file)
  
  test_pt <- st_as_sf(test_site, coords = c('lon', 'lat'), crs = 4326)
  print(mapview(test_basin, layer.name = 'watershed') +
          mapview(test_pt, col.regions = 'red', layer.name = 'site'))
} else {
  message('TEST FAILED: no basin returned for COMID ', test_comid)
}

message('Test complete. Review the map, then run the full loop below.')


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
#     basin <- get_nldi_basin(list(featureSource = 'comid', featureID = cid))
#     if(nrow(basin) == 0) stop('empty basin')
#     st_write(basin, outfile, delete_dsn = TRUE, quiet = TRUE)
#   }, silent = TRUE)
# 
#   if(inherits(z, 'try-error')){
#     outcomes <- add_row(outcomes, comid = cid, status = 'error',
#                         msg = as.character(z))
#   } else {
#     outcomes <- add_row(outcomes, comid = cid, status = 'success', msg = '')
#   }
# 
#   Sys.sleep(0.3)  # be polite to the NLDI API
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
