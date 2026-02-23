library(tidyverse)
library(sf)
library(readxl)
library(lubridate)
library(trend)
library(ggplot2)
library(terra)
library(gifski)

# gages <- read_xlsx('data/allgagesS.xlsx')
# chem <- read_xlsx('data/ChemConLoc.xlsx')
# chem <- read_csv('data/chemcon_snapped.csv')
chem_se <- read_xlsx('data/SeF.xlsx')
chem_so4 <- read_xlsx('data/SulfateF.xlsx')
chem_no3 <- read_xlsx('data/NitrateF.xlsx') %>% 
  filter(! Unit == 'mg N/l******')
chem_cond <- read_xlsx('data/conductivity.xlsx') %>% 
  rename(MonitoringLocationIdentifier = STATION_NUMBER)
glc_smry <- read_csv('data/glc_watershed_summary.csv')
macro <- read_csv('data/macro_snapped0.csv')
crosswalk <- read_csv('data/watersheds_nhdhr/site_snap_distances.csv') %>% 
  mutate(site_id = if_else(grepl('^[0-9]+$', site_id),
                           str_pad(site_id, pad = '0', side = 'left', width = 7),
                           site_id))

## sort out unit basis ambiguity for nitrate

#SAMPLING_AGENCY "Contractor" is ambiguous
no3_smry <- chem_no3 %>%
  group_by(SAMPLING_AGENCY, Unit, CharacteristicName) %>%
  summarize(med = median(Value, na.rm = T),
            max = max(Value, na.rm = T),
            count = n(),
            .groups = 'drop')

chem_no3 <- chem_no3 %>% 
  left_join(no3_smry, by = c('SAMPLING_AGENCY', 'Unit', 'CharacteristicName')) %>% 
  mutate(Value = if_else(Unit == 'mg/l as N' | med < 0.2 | max < 30,
                         Value * 4.426,
                         Value)) %>% 
  select(-Unit, -CharacteristicName)


## a little more cleanup ####

cleanup <- function(d){
  d %>% 
    group_by(MonitoringLocationIdentifier, Year) %>%
    summarize(val = median(Value, na.rm = TRUE),
              .groups = 'drop') %>% 
    mutate(Year = as.numeric(Year)) %>% 
    group_by(MonitoringLocationIdentifier) %>% 
    filter(n() >= 10) %>% 
    arrange(Year) %>% 
    ungroup()
}

chem_no3 <- cleanup(chem_no3)
chem_so4 <- cleanup(chem_so4)
chem_se <- cleanup(chem_se)
chem_cond <- cleanup(chem_cond)

## detect a few watersheds with strong trends ####

get_trends <- function(d){
  trends <- d %>%
    group_by(MonitoringLocationIdentifier) %>%
    summarize(
      sen_slope = sens.slope(val)$estimates,
      p_value   = sens.slope(val)$p.value,
      .groups = 'drop'
    )
  
  top_10 <- trends %>% 
    arrange(desc(sen_slope)) %>% 
    slice(1:10) %>% 
    pull(MonitoringLocationIdentifier)
  
  p <- d %>% 
    filter(MonitoringLocationIdentifier %in% top_10) %>% 
    ggplot(aes(x = Year, y = val)) +
    geom_line() +
    facet_wrap(~MonitoringLocationIdentifier, scales = 'free')
  
  print(p)
  return(arrange(trends, desc(sen_slope)))
}

trends_se <- get_trends(chem_se)
trends_no3 <- get_trends(chem_no3)
trends_so4 <- get_trends(chem_so4)
trends_cond <- get_trends(chem_cond)

selections <- tibble(
  Se = c('E295210', '0200201', '0200311'),
  NO3 = c('0200102', 'E288270', 'E295210'),
  SO4 = c('0200097', 'E295210', '0200209')
  # cond = c()
) %>% 
  pivot_longer(cols = everything(),
               names_to = 'var',
               values_to = 'site_id') %>% 
  arrange(var)


## map selected watersheds ####

comids_to_map <- crosswalk %>%
  filter(site_id %in% selections$site_id) %>% 
  pull(comid) %>% 
  unique()

sheds <- map_dfr(paste0('data/watersheds_nhdhr/sheds_good/comid_', comids_to_map, '.gpkg'), st_read)

# reclassification matrix: GLC_FCS30D class -> broad class integer
# 1=forest, 2=shrubland, 3=grassland, 4=cropland,
# 5=wetland, 6=water, 7=urban, 8=barren_snow_ice
rcl_mat <- matrix(c(
  10,  4,   20,  4,
  51,  1,   52,  1,   61,  1,   62,  1,
  71,  1,   72,  1,   81,  1,   82,  1,
  91,  1,   92,  1,
  120, 2,  121, 2,  122, 2,
  130, 3,  140, 3,
  150, 8,  152, 2,  153, 3,
  180, 5,  190, 7,
  200, 8,  201, 8,  202, 8,
  210, 6,  220, 8,  250, NA
), ncol = 2, byrow = TRUE)

broad_labels <- c(
  'forest', 'shrubland', 'grassland', 'cropland',
  'wetland', 'water', 'impervious', 'barren_snow_ice'
)

broad_colors <- c(
  'forest'          = '#1a9641',
  'shrubland'       = '#a6d96a',
  'grassland'       = '#ffffbf',
  'cropland'        = '#fdae61',
  'wetland'         = '#6baed6',
  'water'           = '#08519c',
  'impervious'      = '#d7191c',
  'barren_snow_ice' = '#bdbdbd'
)

# --- animated GLC land cover map for one watershed ---

poc_comid <- comids_to_map[1]
poc_shed  <- sheds %>% filter(comid == poc_comid) %>% st_transform(4326)
ws_vect   <- vect(poc_shed)
ws_ext    <- ext(st_bbox(poc_shed)) + 0.01

poc_site <- crosswalk %>%
  filter(comid == poc_comid) %>%
  pull(site_id) %>%
  first()

# load Elk-Kootenai basin boundary (for inset)
basin_outline <- st_read('data/elk_kootenai_basin/KootenaiShape.shp', quiet = TRUE) %>%
  st_transform(4326) %>%
  st_make_valid()

# tile configuration (same as summarize_glc.R)
tile_lons <- c('W115', 'W120')
tile_lats <- c('N50', 'N55')
glc_dir   <- 'data/glc'

# early files: 4 bands = 1985, 1990, 1995, 2000
early_years <- c(1985L, 1990L, 1995L, 2000L)
# late files: 23 bands = 2000:2022
late_years  <- 2000L:2022L

# all years to animate (use 2000 from late file to avoid duplication)
all_years <- sort(unique(c(1985L, 1990L, 1995L, late_years)))

# build tile file lookup for both periods
early_files <- expand_grid(lon = tile_lons, lat = tile_lats) %>%
  mutate(path = file.path(glc_dir,
    paste0('GLC_FCS30D_19852000_', lon, lat, '_5years_V1.1.tif'))) %>%
  filter(file.exists(path)) %>%
  pull(path)

late_files <- expand_grid(lon = tile_lons, lat = tile_lats) %>%
  mutate(path = file.path(glc_dir,
    paste0('GLC_FCS30D_20002022_', lon, lat, '_Annual_V1.1.tif'))) %>%
  filter(file.exists(path)) %>%
  pull(path)

# helper: load, crop, mosaic, reclassify, mask for one year
get_lc_year <- function(yr) {
  if (yr %in% late_years) {
    band  <- which(late_years == yr)
    files <- late_files
  } else {
    band  <- which(early_years == yr)
    files <- early_files
  }

  crops <- lapply(files, function(f) {
    r  <- rast(f, lyrs = band)
    ov <- intersect(ext(r), ws_ext)
    if (!is.null(ov)) crop(r, ws_ext) else NULL
  })
  crops <- crops[!sapply(crops, is.null)]
  if (length(crops) == 0) return(NULL)

  yr_rast <- if (length(crops) > 1) do.call(merge, crops) else crops[[1]]
  yr_rcl  <- classify(yr_rast, rcl_mat, others = NA)
  yr_mask <- mask(crop(yr_rcl, ws_vect), ws_vect)
  levels(yr_mask) <- data.frame(id = 1:8, label = broad_labels)
  yr_mask
}

# write individual PNG frames
frame_dir <- file.path(tempdir(), 'glc_frames')
dir.create(frame_dir, showWarnings = FALSE, recursive = TRUE)
frame_paths <- character()

for (yr in all_years) {
  lc <- get_lc_year(yr)
  if (is.null(lc)) next

  fpath <- file.path(frame_dir, sprintf('frame_%04d.png', yr))
  png(fpath, width = 1200, height = 800, res = 150)

  # inset on left (panel 1), land cover on right (panel 2)
  layout(matrix(c(1, 2), nrow = 1), widths = c(1.2, 3))

  # panel 1: basin inset
  par(mar = c(3.1, 2, 3.1, 0.5))
  plot(st_geometry(basin_outline), col = 'grey90', border = 'grey50',
       main = 'Location in\nElk-Kootenai Basin', cex.main = 0.9)
  plot(st_geometry(poc_shed), col = '#d7191c55', border = '#d7191c',
       lwd = 2, add = TRUE)

  # panel 2: land cover
  plot(lc, col = broad_colors[broad_labels],
       main = paste0('GLC Land Cover ', yr,
                     '\nComid ', poc_comid, ' (site ', poc_site, ')'),
       mar = c(3.1, 3.1, 3.1, 8))

  layout(1)
  dev.off()

  frame_paths <- c(frame_paths, fpath)
  message('  ', yr)
}

# assemble GIF
gif_out <- paste0('figs/glc_landcover_comid_', poc_comid, '.gif')
gifski(frame_paths, gif_file = gif_out, width = 1200, height = 800, delay = 0.5)
