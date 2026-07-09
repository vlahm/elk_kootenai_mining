# snap_harmonized_sites.R
# Rebuild the site -> watershed-comid crosswalk for the HARMONIZED site list
# (ElkChem chemistry + ElkMacro macro), snapping to the CORRECT
# Elk-Kootenai NHDPlus-HR flowlines (data/nhdplushr/elk_kootenai/, VPUID 1701;
# the top-level nhdplushr_merged.gpkg has been overwritten with NH/ME data).
#
# This is the same snapping concept as delineate_watersheds.R (nearest HR
# flowline -> its COMID), done with sf::st_nearest_feature so it needs no
# nhdplusTools. Output goes to a NEW path so existing crosswalks are untouched.

library(tidyverse)
library(sf)
library(readxl)
sf_use_s2(FALSE)

out_dir <- 'data/watersheds_nhdhr/harmonized'
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# ── harmonized monitoring sites (distinct site_id + coordinates) ──────────────
chem_sites <- suppressWarnings(read_xlsx('data/harmonized/ElkChem.xlsx')) %>%
  transmute(site_id = STATION_NUMBER, lon = longitude, lat = latitude, source = 'chem') %>%
  filter(!is.na(lon), !is.na(lat)) %>% distinct(site_id, .keep_all = TRUE)

# ElkMacro only (matches kitchen_sink_plots.R; ElkMacro carries the full Elk macro
# network, and MacroC was found to add only redundant/Kootenai-side sites).
macro_sites <- suppressWarnings(read_xlsx('data/harmonized/ElkMacro.xlsx')) %>%
  transmute(site_id = STATION_NUMBER, lon = LongitudeMeasure, lat = LatitudeMeasure) %>%
  filter(!is.na(lon), !is.na(lat)) %>%
  distinct(site_id, .keep_all = TRUE) %>% mutate(source = 'macro')

# union of all sites (a site_id in both chem & macro keeps one row / one comid)
sites <- bind_rows(chem_sites, macro_sites) %>%
  group_by(site_id) %>%
  summarize(lon = first(lon), lat = first(lat),
            source = paste(sort(unique(source)), collapse = '+'), .groups = 'drop')
message('harmonized distinct sites to snap: ', nrow(sites),
        ' (chem ', nrow(chem_sites), ', macro ', nrow(macro_sites), ')')

sites_sf <- st_as_sf(sites, coords = c('lon', 'lat'), crs = 4326) %>% st_transform(3005)

# ── correct Elk-Kootenai HR flowlines ─────────────────────────────────────────
flines <- st_read('data/nhdplushr/elk_kootenai/nhdplushr_merged.gpkg',
                  layer = 'NHDFlowline', quiet = TRUE) %>%
  st_zm(drop = TRUE) %>% st_transform(3005)
if (!'COMID' %in% names(flines) && 'NHDPlusID' %in% names(flines)) flines$COMID <- flines$NHDPlusID
message('flowlines: ', nrow(flines))

# ── snap: nearest flowline -> its COMID + snap distance (m) ───────────────────
idx <- st_nearest_feature(sites_sf, flines)
sites$comid <- as.numeric(flines$COMID[idx])
sites$snap_dist_m <- as.numeric(st_distance(sites_sf, flines[idx, ], by_element = TRUE))

write_csv(select(sites, site_id, comid, snap_dist_m, source),
          file.path(out_dir, 'site_comid_lookup.csv'))

message('\n--- snap-distance summary (m) ---')
message('  min ', round(min(sites$snap_dist_m), 1),
        ' | median ', round(median(sites$snap_dist_m), 1),
        ' | mean ', round(mean(sites$snap_dist_m), 1),
        ' | max ', round(max(sites$snap_dist_m), 1))
message('  sites > 500 m from a flowline (suspect): ', sum(sites$snap_dist_m > 500))
message('  distinct comids hit: ', n_distinct(sites$comid))

# ── verification against the two figure watersheds ────────────────────────────
for (cc in c(55001100327683, 55001100073395)) {
  n <- sum(sites$comid == cc)
  message('  sites snapping to comid ', format(cc, scientific = FALSE), ': ', n,
          if (n > 0) paste0(' -> ', paste(sites$site_id[sites$comid == cc], collapse = ', ')) else '')
}
message('\nWrote ', file.path(out_dir, 'site_comid_lookup.csv'))
