
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

# --- prepare Elk land cover time series (Panel 1) ---

elk_lc <- basin_lc %>%
  filter(basin == 'Elk') %>%
  select(-basin) %>%
  pivot_longer(-year, names_to = 'class', values_to = 'fraction') %>%
  mutate(class = factor(class, levels = broad_class_labels))

lc_colors <- c(
  'forest'          = '#1a9641',
  'shrubland'       = '#a6d96a',
  'grassland'       = '#ffffbf',
  'cropland'        = '#fdae61',
  'wetland'         = '#6baed6',
  'water'           = '#08519c',
  'impervious'      = '#d7191c',
  'barren_snow_ice' = '#bdbdbd'
)

# --- prepare mine area time series (Panel 2) ---

mines <- mines %>% st_transform(st_crs(elk_basin))

# DN = year first mined; cumulative area up to each year
mine_years <- sort(unique(mines$DN))

# disable S2 to avoid duplicate-edge errors in st_union
sf_use_s2(FALSE)

mine_area_ts <- tibble(
  year = mine_years,
  area_km2 = sapply(mine_years, function(y) {
    cum <- mines %>% filter(DN <= y) %>% st_union() %>% st_make_valid()
    as.numeric(st_area(cum)) / 1e6
  })
)

sf_use_s2(TRUE)

# --- prepare barren-mine overlap time series (Panel 3) ---

# For each GLC year, extract barren/snow/ice pixels as a polygon,
# intersect with cumulative mine extent, compute fraction.

elk_basin_4326 <- elk_basin %>% st_transform(4326)
elk_ext <- ext(st_bbox(elk_basin_4326)) + 0.01

overlap_results <- list()

for (yr in all_years) {

  if (yr %in% early_years) {
    period <- 'early'
    band   <- which(early_all == yr)
  } else {
    period <- 'late'
    band   <- which(late_years == yr)
  }

  # load and mosaic tiles for this year
  rasters <- list()
  for (i in seq_len(nrow(tile_combos))) {
    f <- find_tile(period, tile_combos$lon[i], tile_combos$lat[i])
    if (is.na(f) || !file.exists(f)) next
    r <- rast(f, lyrs = band)
    ov <- intersect(ext(r), elk_ext)
    if (!is.null(ov)) rasters <- c(rasters, list(crop(r, elk_ext)))
  }
  if (length(rasters) == 0) next

  yr_rast <- if (length(rasters) > 1) do.call(merge, rasters) else rasters[[1]]
  yr_rcl  <- classify(yr_rast, rcl, others = NA)

  # mask to Elk basin
  elk_vect <- vect(elk_basin_4326)
  yr_elk   <- mask(crop(yr_rcl, elk_vect), elk_vect)

  # total barren/snow/ice area (class 8)
  barren_mask <- yr_elk == 8
  barren_area <- global(barren_mask * cellSize(barren_mask, unit = 'km'),
                        'sum', na.rm = TRUE)[[1]]

  # cumulative mine extent up to this year (use nearest prior mine year)
  valid_mine_yrs <- mine_years[mine_years <= yr]
  if (length(valid_mine_yrs) == 0) {
    overlap_results[[as.character(yr)]] <- tibble(
      year = yr, barren_km2 = barren_area, overlap_km2 = 0, pct_overlap = 0
    )
    rm(yr_rast, yr_rcl, yr_elk, barren_mask); gc(verbose = FALSE)
    next
  }

  cum_mines <- mines %>%
    filter(DN <= yr) %>%
    st_union() %>%
    st_make_valid() %>%
    st_transform(4326)

  # rasterize mine extent onto the barren grid
  mine_vect <- vect(cum_mines)
  mine_rast <- rasterize(mine_vect, yr_elk, field = 1, background = 0)

  # overlap = barren AND mine
  overlap_mask <- barren_mask * mine_rast
  overlap_area <- global(overlap_mask * cellSize(overlap_mask, unit = 'km'),
                         'sum', na.rm = TRUE)[[1]]

  overlap_results[[as.character(yr)]] <- tibble(
    year        = yr,
    barren_km2  = barren_area,
    overlap_km2 = overlap_area,
    pct_overlap = ifelse(barren_area > 0, overlap_area / barren_area * 100, NA_real_)
  )

  message('  Overlap ', yr, ' done')
  rm(yr_rast, yr_rcl, yr_elk, barren_mask, mine_rast, overlap_mask)
  gc(verbose = FALSE)
}

overlap_ts <- bind_rows(overlap_results)

# --- 4-panel figure ---

library(patchwork)

out_png <- 'figs/elk_landcover_mining.png'
dir.create(dirname(out_png), showWarnings = FALSE, recursive = TRUE)

# reorder classes so forest (dominant) is on bottom of stack
class_order <- c('forest', 'shrubland', 'grassland', 'cropland',
                 'wetland', 'water', 'impervious', 'barren_snow_ice')

elk_lc <- elk_lc %>%
  mutate(class = factor(class, levels = class_order),
         pct = fraction * 100)

# Panel 1: stacked bar of full composition
p1 <- ggplot(elk_lc, aes(x = year, y = pct, fill = class)) +
  geom_col(width = 1) +
  scale_fill_manual(values = lc_colors, name = 'Land cover') +
  labs(x = NULL, y = '% of basin', title = 'Elk Basin land cover composition') +
  theme_minimal(base_size = 11) +
  theme(legend.position = 'right',
        axis.text.x = element_blank(),
        plot.title = element_text(face = 'bold', size = 12))

# Panel 2: change from 1985 baseline (percentage points)
baseline <- elk_lc %>%
  filter(year == min(year)) %>%
  select(class, baseline_pct = pct)

elk_lc_delta <- elk_lc %>%
  left_join(baseline, by = 'class') %>%
  mutate(delta = pct - baseline_pct)

p2 <- ggplot(elk_lc_delta, aes(x = year, y = delta, color = class)) +
  geom_hline(yintercept = 0, linetype = 'dashed', color = 'grey50') +
  geom_line(linewidth = 0.8) +
  geom_point(size = 0.8) +
  scale_color_manual(values = lc_colors, name = 'Land cover') +
  labs(x = NULL, y = 'Change (%)',
       title = 'Change from 1985 baseline (%)') +
  theme_minimal(base_size = 11) +
  theme(legend.position = 'right',
        axis.text.x = element_blank(),
        plot.title = element_text(face = 'bold', size = 12))

# Panel 3: cumulative mine area
p3 <- ggplot(mine_area_ts, aes(x = year, y = area_km2)) +
  geom_line(color = '#e41a1c', linewidth = 0.8) +
  geom_point(color = '#e41a1c', size = 0.8) +
  labs(x = NULL, y = expression('Cumulative area (km'^2*')'),
       title = 'Cumulative mining extent') +
  theme_minimal(base_size = 11) +
  theme(axis.text.x = element_blank(),
        plot.title = element_text(face = 'bold', size = 12))

# Panel 4: % barren overlapping mines
p4 <- ggplot(overlap_ts, aes(x = year, y = pct_overlap)) +
  geom_line(color = '#984ea3', linewidth = 0.8) +
  geom_point(color = '#984ea3', size = 0.8) +
  labs(x = 'Year', y = 'overlap (%)',
       title = 'Barren/snow/ice area overlapping mining extent') +
  theme_minimal(base_size = 11) +
  theme(plot.title = element_text(face = 'bold', size = 12))

combined <- (p1 / p2 / p3 / p4) +
  plot_annotation(
    title = 'Elk River Basin — Land Cover & Mining',
    theme = theme(plot.title = element_text(face = 'bold', size = 14, hjust = 0.5))
  ) +
  plot_layout(heights = c(3, 3, 2, 2))

ggsave(out_png, combined, width = 10, height = 14, dpi = 300)
