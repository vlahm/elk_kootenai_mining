library(tidyverse)
library(sf)
library(readxl)
library(lubridate)
library(ggplot2)
library(terra)
library(patchwork)
library(trend)

# ── load data ──────────────────────────────────────────────────────────────────

glc_smry <- read_csv('data/glc_watershed_summary.csv') %>%
  rename(impervious = urban)

macro <- read_csv('data/macro_snapped0.csv') %>%
  rename(site_id = SiteNumber) %>%
  mutate(site_id = if_else(grepl('^[0-9]+$', site_id),
                           str_pad(site_id, pad = '0', side = 'left', width = 7),
                           site_id)) %>%
  select(site_id:Count)

crosswalk <- read_csv('data/watersheds_nhdhr/site_snap_distances.csv') %>%
  mutate(site_id = if_else(grepl('^[0-9]+$', site_id),
                           str_pad(site_id, pad = '0', side = 'left', width = 7),
                           site_id))

macro <- macro %>%
  left_join(filter(crosswalk, source == 'macro'), by = 'site_id')

mines_raw <- st_read('data/skytruth/med_0-10_thresh_ridgeMasked_1985-2024_firstMined.geojson')

# ── orders of interest ────────────────────────────────────────────────────────

focal_orders <- c('Ephemeroptera', 'Plecoptera', 'Diptera',
                  'Trichoptera', 'Coleoptera', 'Trombidiformes')

order_colors <- c(

  'Ephemeroptera'  = '#1b9e77',
  'Plecoptera'     = '#d95f02',
  'Diptera'        = '#7570b3',
  'Trichoptera'    = '#e7298a',
  'Coleoptera'     = '#66a61e',
  'Trombidiformes' = '#e6ab02'
)

# ── find 3 watersheds with strongest declining genus richness trends ──────────

macro_year_counts <- macro %>%
  filter(!is.na(comid)) %>%
  group_by(comid) %>%
  summarize(n_macro_years = length(unique(Year)), .groups = 'drop') %>%
  filter(n_macro_years >= 3)

richness_trends <- macro %>%
  filter(!is.na(comid),
         comid %in% macro_year_counts$comid) %>%
  group_by(comid, Year) %>%
  summarize(richness = length(unique(na.omit(Genus))), .groups = 'drop') %>%
  group_by(comid) %>%
  filter(n() >= 3) %>%
  arrange(Year) %>%
  summarize(
    sen_slope = sens.slope(richness)$estimates,
    p_value   = sens.slope(richness)$p.value,
    n_years   = n(),
    .groups   = 'drop'
  ) %>%
  arrange(sen_slope)

# keep only comids that have a watershed file on disk
richness_trends <- richness_trends %>%
  filter(file.exists(paste0('data/watersheds_nhdhr/sheds_good/comid_', comid, '.gpkg')))

# also require GLC summary data
richness_trends <- richness_trends %>%
  filter(comid %in% glc_smry$comid)

# ── pre-compute mining % for all viable comids ───────────────────────────────

viable_comids <- richness_trends$comid

message('Computing mining % for ', length(viable_comids), ' viable COMIDs ...')

all_sheds <- map_dfr(paste0('data/watersheds_nhdhr/sheds_good/comid_',
                            viable_comids, '.gpkg'), st_read) %>%
  st_transform(4326)

sf_use_s2(FALSE)
mining_pct <- map_dfr(viable_comids, function(cid) {
  shed <- all_sheds %>% filter(comid == cid)
  ws_area_km2 <- as.numeric(st_area(st_union(shed))) / 1e6
  mines_ws <- tryCatch(
    st_intersection(mines_raw, shed) %>% st_make_valid(),
    error = function(e) st_sf(geometry = st_sfc(), crs = 4326)
  )
  if (nrow(mines_ws) > 0) {
    mine_area_km2 <- as.numeric(st_area(st_union(mines_ws))) / 1e6
  } else {
    mine_area_km2 <- 0
  }
  tibble(comid = cid,
         ws_area_km2 = ws_area_km2,
         mine_area_km2 = mine_area_km2,
         mine_pct = mine_area_km2 / ws_area_km2 * 100)
})
sf_use_s2(TRUE)

message('Mining % summary:')
print(summary(mining_pct$mine_pct))

# keep only comids with some mining influence, then pick top 3 declining richness
comids_with_mining <- mining_pct %>%
  filter(mine_pct > 0) %>%
  pull(comid)

message(length(comids_with_mining), ' of ', length(viable_comids),
        ' viable COMIDs have mining influence')

comids_to_map <- richness_trends %>%
  filter(comid %in% comids_with_mining) %>%
  slice_head(n = 10) %>%
  pull(comid)

message('Selected COMIDs (declining richness + mining): ',
        paste(comids_to_map, collapse = ', '))

sheds <- all_sheds %>% filter(comid %in% comids_to_map)

# ── GLC / tile setup ──────────────────────────────────────────────────────────

tile_lons <- c('W115', 'W120')
tile_lats <- c('N50', 'N55')
glc_dir   <- 'data/glc'

early_years <- c(1985L, 1990L, 1995L, 2000L)
late_years  <- 2000L:2022L
all_years   <- sort(unique(c(1985L, 1990L, 1995L, late_years)))

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

broad_labels <- c('forest', 'shrubland', 'grassland', 'cropland',
                  'wetland', 'water', 'impervious', 'barren_snow_ice')

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

# ── SkyTruth mines ─────────────────────────────────────────────────────────────

mines_raw <- st_read('data/skytruth/med_0-10_thresh_ridgeMasked_1985-2024_firstMined.geojson')

# ── basin outline (for inset) ─────────────────────────────────────────────────

basin_outline <- st_read('data/elk_kootenai_basin/KootenaiShape.shp', quiet = TRUE) %>%
  st_transform(4326) %>%
  st_make_valid()

# ── helper: get reclassified land cover raster for one year ───────────────────

get_lc_year <- function(yr, ws_ext, ws_vect) {
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

# ── main loop ─────────────────────────────────────────────────────────────────

dir.create('figs', showWarnings = FALSE, recursive = TRUE)

for (i in seq_along(comids_to_map)) {

  poc_comid <- comids_to_map[i]
  poc_shed  <- sheds %>% filter(comid == poc_comid) %>% st_transform(4326)
  ws_vect   <- vect(poc_shed)
  ws_ext    <- ext(st_bbox(poc_shed)) + 0.01

  message('\n=== Static plot for COMID ', poc_comid, ' ===')

  # ── 1. Land cover change from 1985 baseline ────────────────────────────────

  cur_lc <- glc_smry %>%
    filter(comid == poc_comid) %>%
    select(year, all_of(broad_labels)) %>%
    pivot_longer(-year, names_to = 'class', values_to = 'fraction') %>%
    mutate(class = factor(class, levels = broad_labels),
           pct = fraction * 100)

  baseline <- cur_lc %>%
    filter(year == min(year)) %>%
    select(class, baseline_pct = pct)

  cur_lc_delta <- cur_lc %>%
    left_join(baseline, by = 'class') %>%
    mutate(delta = pct - baseline_pct)

  p1 <- ggplot(cur_lc_delta, aes(x = year, y = delta, color = class)) +
    geom_hline(yintercept = 0, linetype = 'dashed', color = 'grey50') +
    geom_line(linewidth = 0.8) +
    geom_point(size = 0.8) +
    scale_color_manual(values = broad_colors, name = 'Land cover') +
    labs(x = NULL, y = 'Change from 1985 (%)',
         title = 'Land cover change') +
    theme_minimal(base_size = 11) +
    theme(legend.position = 'none',
          axis.text.x = element_blank(),
          plot.title = element_text(face = 'bold', size = 11))

  # ── 2. Cumulative mining extent ─────────────────────────────────────────────

  mines_ws <- st_intersection(mines_raw, poc_shed) %>% st_make_valid()
  ws_area_km2 <- as.numeric(st_area(st_union(poc_shed))) / 1e6

  if (nrow(mines_ws) > 0) {
    mine_yrs_ws <- sort(unique(mines_ws$DN))
    sf_use_s2(FALSE)
    ts_mining <- tibble(
      year = mine_yrs_ws,
      area_km2 = sapply(mine_yrs_ws, function(y) {
        cum <- mines_ws %>% filter(DN <= y) %>% st_union() %>% st_make_valid()
        as.numeric(st_area(cum)) / 1e6
      })
    )
    sf_use_s2(TRUE)
    ts_mining <- tibble(year = as.numeric(all_years)) %>%
      left_join(ts_mining, by = 'year') %>%
      arrange(year) %>%
      mutate(area_km2 = ifelse(year < min(mine_yrs_ws), 0, area_km2)) %>%
      tidyr::fill(area_km2, .direction = 'down')
  } else {
    ts_mining <- tibble(year = as.numeric(all_years), area_km2 = 0)
  }

  p2 <- ggplot(ts_mining, aes(x = year, y = area_km2)) +
    geom_line(color = '#e41a1c', linewidth = 0.8) +
    geom_point(color = '#e41a1c', size = 0.8) +
    labs(x = 'Year', y = expression('Cum. area (km'^2*')'),
         title = 'Cumulative mining extent') +
    theme_minimal(base_size = 11) +
    theme(plot.title = element_text(face = 'bold', size = 11))

  # ── 3. Invertebrate density ─────────────────────────────────────────────────

  ts_macro <- macro %>%
    filter(comid == poc_comid) %>%
    group_by(Year) %>%
    summarize(density  = sum(Count, na.rm = TRUE),
              richness = length(unique(na.omit(Genus))),
              .groups = 'drop') %>%
    arrange(Year)

  p3 <- ggplot(ts_macro, aes(x = Year, y = density)) +
    geom_line(color = '#984ea3', linewidth = 0.8) +
    geom_point(color = '#984ea3', size = 0.8) +
    labs(x = NULL, y = 'Count',
         title = 'Macroinvertebrate density (top 6 orders)') +
    theme_minimal(base_size = 11) +
    theme(axis.text.x = element_blank(),
          plot.title = element_text(face = 'bold', size = 11))

  # ── 4. Invertebrate genus richness ──────────────────────────────────────────

  p4 <- ggplot(ts_macro, aes(x = Year, y = richness)) +
    geom_line(color = '#ff7f00', linewidth = 0.8) +
    geom_point(color = '#ff7f00', size = 0.8) +
    labs(x = 'Year', y = 'Count',
         title = 'Genus richness (unique genera)') +
    theme_minimal(base_size = 11) +
    theme(plot.title = element_text(face = 'bold', size = 11))

  # ── 5. Density by taxonomic order ──────────────────────────────────────────

  ts_order_density <- macro %>%
    filter(comid == poc_comid, Order_ %in% focal_orders) %>%
    group_by(Year, Order_) %>%
    summarize(density = sum(Count, na.rm = TRUE), .groups = 'drop') %>%
    arrange(Year)

  p5 <- ggplot(ts_order_density, aes(x = Year, y = density, color = Order_)) +
    geom_line(linewidth = 0.8) +
    geom_point(size = 0.8) +
    scale_color_manual(values = order_colors, name = 'Order') +
    labs(x = NULL, y = 'Count',
         title = 'Density by order') +
    theme_minimal(base_size = 11) +
    theme(legend.position = 'none',
          axis.text.x = element_blank(),
          plot.title = element_text(face = 'bold', size = 11))

  # ── 6. Richness by taxonomic order ─────────────────────────────────────────

  ts_order_richness <- macro %>%
    filter(comid == poc_comid, Order_ %in% focal_orders) %>%
    group_by(Year, Order_) %>%
    summarize(richness = length(unique(na.omit(Genus))), .groups = 'drop') %>%
    arrange(Year)

  p6 <- ggplot(ts_order_richness, aes(x = Year, y = richness, color = Order_)) +
    geom_line(linewidth = 0.8) +
    geom_point(size = 0.8) +
    scale_color_manual(values = order_colors, name = 'Order') +
    labs(x = 'Year', y = 'Count',
         title = 'Richness by order') +
    theme_minimal(base_size = 11) +
    theme(legend.position = c(0.85, 0.85),
          legend.background = element_rect(fill = 'white', color = NA),
          legend.text = element_text(size = 7),
          legend.title = element_text(size = 8),
          legend.key.size = unit(0.4, 'cm'),
          plot.title = element_text(face = 'bold', size = 11))

  # ── right column: inset map, legend, 1985 & 2022 land cover ────────────────

  # inset map
  bb <- st_bbox(basin_outline)
  p_inset <- ggplot() +
    geom_sf(data = basin_outline, fill = 'grey90', color = 'grey50') +
    geom_sf(data = poc_shed, fill = '#d7191c55', color = '#d7191c', linewidth = 0.8) +
    coord_sf(xlim = c(bb['xmin'], bb['xmax']),
             ylim = c(bb['ymin'], bb['ymax'])) +
    labs(title = paste0('COMID ', poc_comid)) +
    theme_void(base_size = 11) +
    theme(plot.title = element_text(face = 'bold', size = 11, hjust = 0.5))

  # land cover legend
  present_classes <- cur_lc %>%
    filter(pct > 0) %>%
    pull(class) %>%
    unique() %>%
    as.character()

  p_legend <- ggplot(tibble(class = factor(present_classes, levels = broad_labels),
                            y = seq_along(present_classes)),
                     aes(x = 1, y = y, fill = class)) +
    geom_tile(width = 0.3, height = 0.8) +
    geom_text(aes(x = 1.25, label = class), hjust = 0, size = 3) +
    scale_fill_manual(values = broad_colors, guide = 'none') +
    xlim(0.8, 2.5) +
    theme_void() +
    labs(title = 'Land cover classes') +
    theme(plot.title = element_text(face = 'bold', size = 10, hjust = 0.5))

  # static land cover maps for 1985 and 2022
  make_lc_ggplot <- function(yr) {
    lc <- get_lc_year(yr, ws_ext, ws_vect)
    if (is.null(lc)) {
      return(ggplot() + theme_void() +
               labs(title = paste0(yr, ' — no data')) +
               theme(plot.title = element_text(hjust = 0.5, size = 10)))
    }
    lc_df <- as.data.frame(lc, xy = TRUE)
    names(lc_df)[3] <- 'class_id'
    lc_df <- lc_df %>%
      filter(!is.na(class_id)) %>%
      mutate(class = factor(broad_labels[class_id], levels = broad_labels))

    ggplot(lc_df, aes(x = x, y = y, fill = class)) +
      geom_raster() +
      geom_sf(data = poc_shed, fill = NA, color = 'black',
              linewidth = 0.4, inherit.aes = FALSE) +
      scale_fill_manual(values = broad_colors, guide = 'none', drop = FALSE) +
      coord_sf(expand = FALSE) +
      labs(title = yr) +
      theme_void(base_size = 11) +
      theme(plot.title = element_text(face = 'bold', size = 11, hjust = 0.5))
  }

  p_lc1985 <- make_lc_ggplot(1985)
  p_lc2022 <- make_lc_ggplot(2022)

  # ── assemble with patchwork ─────────────────────────────────────────────────

  left_col  <- (p1 / p2 / p5)
  right_col_ts <- (p3 / p4 / p6)
  sidebar   <- (p_inset / p_legend / p_lc1985 / p_lc2022) +
    plot_layout(heights = c(2, 1.5, 2, 2))

  combined <- (left_col | right_col_ts | sidebar) +
    plot_layout(widths = c(2, 2, 1.2)) +
    plot_annotation(
      title = paste0('Watershed ', poc_comid, ' — Land Cover, Mining & Macroinvertebrates'),
      theme = theme(plot.title = element_text(face = 'bold', size = 14, hjust = 0.5))
    )

  out_png <- paste0('figs/macro_lc_comid_', poc_comid, '.png')
  ggsave(out_png, combined, width = 16, height = 10, dpi = 300)
  message('Saved ', out_png)
}
