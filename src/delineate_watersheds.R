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
  rename(lat = Snap_Y, lon = Snap_X) %>% 
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
sites$comid <- discover_nhdplus_id(sites_sf)

# some may fail; report
n_found <- sum(!is.na(sites$comid))
message(n_found, ' of ', nrow(sites), ' sites matched to a COMID')

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
