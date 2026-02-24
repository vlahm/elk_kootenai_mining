
library(tidyverse)
library(sf)
library(lubridate)
library(ggplot2)
library(terra)
library(exactextractr)

kootenai_basin <- st_read('data/elk_kootenai_basin/h8/ek_h8.shp')

elk_basin <- filter(kootenai_basin, huc8 == 17010106)
kootenay_basin <- filter(kootenai_basin, huc8 %in% c(17010107, 17010108))

## build 3-basin sf for GLC extraction ####

basins <- bind_rows(
  elk_basin      %>% summarise(geometry = st_union(geometry)) %>% mutate(basin = 'Elk'),
  kootenai_basin  %>% summarise(geometry = st_union(geometry)) %>% mutate(basin = 'Kootenai'),
  kootenay_basin %>% summarise(geometry = st_union(geometry)) %>% mutate(basin = 'Kootenay')
) %>%
  select(basin, geometry) %>%
  st_make_valid() %>%
  st_transform(4326)

## GLC reclassification table ####

broad_class_labels <- c(
  '1' = 'forest',       '2' = 'shrubland',  '3' = 'grassland',
  '4' = 'cropland',     '5' = 'wetland',    '6' = 'water',
  '7' = 'impervious',        '8' = 'barren_snow_ice'
)

rcl <- matrix(c(
  10,  4,  20,  4,
  51,  1,  52,  1,  61,  1,  62,  1,
  71,  1,  72,  1,  81,  1,  82,  1,  91,  1,  92,  1,
  120, 2, 121, 2, 122, 2,
  130, 3, 140, 3,
  150, 8, 152, 2, 153, 3,
  180, 5, 190, 7,
  200, 8, 201, 8, 202, 8,
  210, 6, 220, 8, 250, NA
), ncol = 2, byrow = TRUE)

## locate GLC tiles ####

glc_dir <- 'data/glc'

# tile grid covering the study area
tile_combos <- expand_grid(lon = c('W115', 'W120'), lat = c('N50', 'N55'))

early_all   <- c(1985L, 1990L, 1995L, 2000L)
early_years <- c(1985L, 1990L, 1995L)
late_years  <- 2000L:2022L
all_years   <- sort(unique(c(early_years, late_years)))

find_tile <- function(period, lon, lat) {
  tag <- paste0(lon, lat)
  pattern <- if (period == 'early') {
    paste0('GLC_FCS30D_1985.*', tag, '.*5years.*\\.tif$')
  } else {
    paste0('GLC_FCS30D_2000.*', tag, '.*Annual.*\\.tif$')
  }
  list.files(glc_dir, pattern = pattern, full.names = TRUE, ignore.case = TRUE)[1]
}

## extract land cover fractions per basin per year ####

ws_extent <- ext(st_bbox(basins)) + 0.01
all_results <- list()

message('Extracting GLC fractions for ', nrow(basins), ' basins x ',
        length(all_years), ' years')

for (yr in all_years) {

  if (yr %in% early_years) {
    period <- 'early'
    band   <- which(early_all == yr)
  } else {
    period <- 'late'
    band   <- which(late_years == yr)
  }

  # load, crop, and mosaic tiles
  rasters <- list()
  for (i in seq_len(nrow(tile_combos))) {
    f <- find_tile(period, tile_combos$lon[i], tile_combos$lat[i])
    if (is.na(f) || !file.exists(f)) next
    r <- rast(f, lyrs = band)
    ov <- intersect(ext(r), ws_extent)
    if (!is.null(ov)) rasters <- c(rasters, list(crop(r, ws_extent)))
  }
  if (length(rasters) == 0) { message('  No tiles for ', yr); next }

  yr_rast <- if (length(rasters) > 1) do.call(merge, rasters) else rasters[[1]]
  yr_rcl  <- classify(yr_rast, rcl, others = NA)

  # extract fractions
  fracs <- exact_extract(yr_rcl, basins, function(values, coverage_fractions) {
    vals <- values[!is.na(values)]
    cov  <- coverage_fractions[!is.na(values)]
    tw   <- sum(cov)
    if (tw == 0) {
      return(as.data.frame(as.list(
        setNames(rep(NA_real_, length(broad_class_labels)), broad_class_labels)
      )))
    }
    class_fracs <- vapply(seq_along(broad_class_labels), function(ci) {
      sum(cov[vals == ci]) / tw
    }, numeric(1))
    as.data.frame(as.list(setNames(class_fracs, broad_class_labels)))
  })

  frac_df <- bind_rows(fracs)
  frac_df$basin <- basins$basin
  frac_df$year  <- yr
  all_results[[as.character(yr)]] <- frac_df

  message('  ', yr, ' done')
  rm(yr_rast, yr_rcl, rasters, fracs, frac_df); gc(verbose = FALSE)
}

## combine and write ####

basin_lc <- bind_rows(all_results) %>%
  select(basin, year, everything()) %>%
  arrange(basin, year)

out_file <- 'data/major_basin_summary.csv'
write_csv(basin_lc, out_file)


## plot ####

mines <- st_read('data/skytruth/med_0-10_thresh_ridgeMasked_1985-2024_firstMined.geojson')
