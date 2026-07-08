library(tidyverse)
library(sf)
library(readxl)
library(lubridate)
library(ggplot2)
library(terra)
library(patchwork)

options(scipen = 999)   # keep 14-digit comids out of scientific notation
sf_use_s2(FALSE)        # planar geometry ops (areas computed in projected CRS below)

# Aggregation mode for the chem & macro panels (land cover / mining are always
# watershed-scale). 'containment' = all monitoring sites inside the shed;
# 'hybrid' = chem from the reach-snapped site(s) + macro from the nearest macro
# site to that reach (site-scale, like the original figure). Set via env var:
#   KS_MODE=hybrid Rscript src/kitchen_sink_plots.R
AGG_MODE   <- Sys.getenv('KS_MODE', 'containment')
out_suffix <- if (AGG_MODE == 'hybrid') '_hybrid' else ''
message('Aggregation mode: ', AGG_MODE)

# ── load data (harmonized inputs) ────────────────────────────────────────────
# Rebuilt to source from the official harmonized files instead of the legacy
# per-analyte pulls:
#   chemistry  -> data/harmonized/ElkChem.xlsx  (carries comid directly)
#   macro      -> data/harmonized/ElkMacro.xlsx (RAEMP) + data/for_reals/MacroC.xlsx (non-RAEMP)
#   mining     -> SkyTruth med_0-7 (Pericak-standard threshold used by the report)
#   watersheds -> data/watersheds_nhdhr/chem_macro_set1/sheds_good (Elk-Kootenai only)
# QA carried over from misc_analyses.R: ElkChem de-duplicated, Se >1000 ug/L dropped.

std_cmp <- function(x) case_when(
  x == 'Selenium' ~ 'Selenium', x == 'Sulfate' ~ 'Sulfate',
  x %in% c('Nitrate', 'NitrateNO3') ~ 'Nitrate',
  x %in% c('Specific conductance', 'Specific Conductivity') ~ 'Conductivity',
  TRUE ~ NA_character_)

# Chemistry (ElkChem) and macro (ElkMacro RAEMP + MacroC non-RAEMP) are kept at
# the SITE level with coordinates. Each watershed's panels aggregate the sites
# that fall INSIDE its delineated shed (spatial containment). This is robust to
# the fact that ElkChem's `comid` key does not match the watershed comids and
# the original site->comid crosswalk (site_snap_distances.csv) no longer exists.
elkchem <- suppressWarnings(read_xlsx('data/harmonized/ElkChem.xlsx')) %>%
  distinct() %>%
  mutate(analyte = std_cmp(Compound)) %>%
  filter(!is.na(analyte), !(analyte == 'Selenium' & Value > 1000),
         !is.na(longitude), !is.na(latitude)) %>%
  transmute(site_id = STATION_NUMBER, lon = longitude, lat = latitude,
            analyte, Value, Year = year)

macro_raemp <- suppressWarnings(read_xlsx('data/harmonized/ElkMacro.xlsx')) %>%
  filter(SAMPLING_AGENCY == 'RAEMP') %>%
  transmute(site_id = STATION_NUMBER, lon = LongitudeMeasure, lat = LatitudeMeasure,
            sample_date = as_date(date), Year = as.integer(year),
            Genus = na_if(as.character(Unit), '-'), Count = Value,
            agency = SAMPLING_AGENCY, program = 'RAEMP')
macro_pub <- suppressWarnings(read_xlsx('data/for_reals/MacroC.xlsx')) %>%
  transmute(site_id = STATION_NUMBER, lon = LongitudeMeasure, lat = LatitudeMeasure,
            sample_date = as_date(SampleDate), Year = year(as_date(SampleDate)),
            Genus = na_if(as.character(Unit), '-'), Count = Value,
            agency = SAMPLING_AGENCY, program = 'non-RAEMP')
macro <- bind_rows(macro_raemp, macro_pub) %>% filter(!is.na(lon), !is.na(lat))

# site points for spatial containment (built once)
chem_pts  <- elkchem %>% distinct(site_id, lon, lat) %>%
  st_as_sf(coords = c('lon', 'lat'), crs = 4326)
macro_pts <- macro %>% distinct(site_id, lon, lat) %>%
  st_as_sf(coords = c('lon', 'lat'), crs = 4326)

# (land cover is computed per-watershed from the GLC rasters in the loop below,
#  so the 110-comid glc_watershed_summary.csv intermediary is no longer needed)

# ── watersheds to render ──────────────────────────────────────────────────────
# The two watersheds behind the existing all_trends figures, rebuilt directly on
# the harmonized inputs for a like-for-like comparison. (To re-run the automatic
# "strongest declining richness" selection instead, restore that block and point
# it at chem_macro_set1/sheds_good.)
comids_to_map <- c(55001100327683, 55001100073395)
if (AGG_MODE == 'hybrid') comids_to_map <- 55001100327683   # the important figure only
message('Rebuilding COMIDs: ', paste(comids_to_map, collapse = ', '))

sheds <- map_dfr(paste0('data/watersheds_nhdhr/chem_macro_set1/sheds_good/comid_',
                        comids_to_map, '.gpkg'), st_read)

# ── GLC / tile setup ──────────────────────────────────────────────────────────

tile_lons <- c('W115', 'W120')
tile_lats <- c('N45', 'N50', 'N55')   # tiles overlapping the Elk-Kootenai watersheds
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

mines_raw <- st_read('data/skytruth/mining_extent_firstMined_1984_2025_med0-7.geojson') %>%
  st_transform(3005) %>% st_make_valid()   # BC Albers -> planar area in km^2, no lwgeom needed

# ── basin outline (for inset) ─────────────────────────────────────────────────

basin_outline <- st_read('data/elk_kootenai_basin/KootenaiShape.shp', quiet = TRUE) %>%
  st_transform(4326) %>%
  st_make_valid()

# ── helper: reclassified land cover raster for one year, cropped to watershed ─
# Reads only the one band needed per year (a whole-stack crop turned out slower:
# merging all 23 bands across tiles forces expensive cross-tile resampling).
# A per-comid memo (below) makes sure each year is computed only once, even
# though 1985 & 2022 are used by both the change panel and the maps.
get_lc_year <- function(yr, ws_ext, ws_vect) {
  band  <- if (yr %in% late_years) which(late_years == yr) else which(early_years == yr)
  files <- if (yr %in% late_years) late_files else early_files
  crops <- lapply(files, function(f) {
    r <- rast(f, lyrs = band)
    if (is.null(intersect(ext(r), ws_ext))) NULL else crop(r, ws_ext)
  })
  crops <- crops[!sapply(crops, is.null)]
  if (length(crops) == 0) return(NULL)
  yr_rast <- if (length(crops) > 1) do.call(merge, crops) else crops[[1]]
  yr_rcl  <- classify(yr_rast, rcl_mat, others = NA)
  yr_mask <- mask(crop(yr_rcl, ws_vect), ws_vect)
  levels(yr_mask) <- data.frame(id = 1:8, label = broad_labels)
  yr_mask
}

# broad-class fractions for one (pre-cropped) year raster
lc_fractions_from <- function(lc, yr) {
  v <- as.vector(terra::values(lc)); v <- v[!is.na(v)]
  if (length(v) == 0) return(NULL)
  ids <- if (is.numeric(v)) v else match(as.character(v), broad_labels)
  tab <- table(factor(ids, levels = 1:8))
  tibble(year = yr, class = factor(broad_labels, levels = broad_labels),
         fraction = as.numeric(tab) / sum(tab))
}

# ── main loop ─────────────────────────────────────────────────────────────────

dir.create('figs', showWarnings = FALSE, recursive = TRUE)

for (i in seq_along(comids_to_map)) {

  poc_comid <- comids_to_map[i]
  poc_shed  <- sheds %>% filter(comid == poc_comid) %>% st_transform(4326)
  ws_vect   <- vect(poc_shed)
  ws_ext    <- ext(st_bbox(poc_shed)) + 0.01
  # per-comid memo so each year's land cover is computed once (reused by the
  # change panel and the 1985/2022 maps)
  lc_env <- new.env()
  get_lc <- function(yr) {
    k <- as.character(yr)
    if (is.null(lc_env[[k]]))
      lc_env[[k]] <- { r <- get_lc_year(yr, ws_ext, ws_vect); if (is.null(r)) 'NODATA' else r }
    v <- lc_env[[k]]; if (identical(v, 'NODATA')) NULL else v
  }

  # monitoring sites physically inside this watershed (containment default)
  chem_here  <- chem_pts$site_id[lengths(st_within(chem_pts, poc_shed)) > 0]
  macro_here <- macro_pts$site_id[lengths(st_within(macro_pts, poc_shed)) > 0]

  if (AGG_MODE == 'hybrid') {
    # chem: the site(s) snapped to THIS reach (harmonized crosswalk); macro: the
    # single nearest macro site to that reach -- reach-scale, like the original.
    xw_h <- read_csv('data/watersheds_nhdhr/harmonized/site_comid_lookup.csv',
                     show_col_types = FALSE)
    reach_sites <- xw_h$site_id[xw_h$comid == poc_comid]
    ch <- intersect(reach_sites, chem_pts$site_id)
    if (length(ch) > 0) chem_here <- ch
    reach_loc <- chem_pts %>% filter(site_id %in% chem_here) %>%
      st_union() %>% st_centroid()
    if (length(reach_loc) == 0) reach_loc <- st_centroid(st_union(poc_shed))
    macro_here <- macro_pts$site_id[st_nearest_feature(reach_loc, macro_pts)]
    message('  hybrid: reach chem sites = ', paste(chem_here, collapse = ', '),
            ' | nearest macro site = ', paste(macro_here, collapse = ', '))
  }

  message('\n=== Static plot for COMID ', poc_comid, ' === (',
          length(chem_here), ' chem sites, ', length(macro_here), ' macro sites)')

  # ── 1. Land cover change from 1985 baseline ────────────────────────────────

  cur_lc <- map_dfr(all_years, ~ { lc <- get_lc(.x)
                                    if (is.null(lc)) NULL else lc_fractions_from(lc, .x) }) %>%
    mutate(pct = fraction * 100)

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
    labs(x = NULL, y = paste0('Change from ', min(cur_lc$year), ' (%)'),
         title = 'Land cover change') +
    theme_minimal(base_size = 11) +
    theme(legend.position = 'none',
          axis.text.x = element_blank(),
          plot.title = element_text(face = 'bold', size = 11))

  # ── 2. Cumulative mining extent ─────────────────────────────────────────────

  poc_shed_p  <- st_make_valid(st_transform(poc_shed, 3005))
  mines_ws    <- suppressWarnings(st_intersection(mines_raw, poc_shed_p)) %>% st_make_valid()
  ws_area_km2 <- as.numeric(st_area(st_union(poc_shed_p))) / 1e6

  if (nrow(mines_ws) > 0) {
    mine_yrs_ws <- sort(unique(mines_ws$DN))
    ts_mining <- tibble(
      year = mine_yrs_ws,
      area_km2 = sapply(mine_yrs_ws, function(y) {
        cum <- mines_ws %>% filter(DN <= y) %>% st_union() %>% st_make_valid()
        as.numeric(st_area(cum)) / 1e6
      })
    )
    ts_mining <- tibble(year = as.numeric(all_years)) %>%
      left_join(ts_mining, by = 'year') %>%
      arrange(year) %>%
      mutate(area_km2 = ifelse(year < min(mine_yrs_ws), 0, area_km2)) %>%
      tidyr::fill(area_km2, .direction = 'down')
  } else {
    ts_mining <- tibble(year = as.numeric(all_years), area_km2 = 0)
  }

  p2 <- ggplot(ts_mining, aes(x = year, y = area_km2)) +
    geom_line(color = 'black', linewidth = 0.8) +
    geom_point(color = 'black', size = 0.8) +
    labs(x = NULL, y = expression('Cum. area (km'^2*')'),
         title = 'Cumulative mining extent') +
    theme_minimal(base_size = 11) +
    theme(axis.text.x = element_blank(),
          plot.title = element_text(face = 'bold', size = 11))

  # ── 3. Selenium ─────────────────────────────────────────────────────────────

  se_here <- function(an) elkchem %>%
    filter(site_id %in% chem_here, analyte == an) %>%
    group_by(Year) %>%
    summarize(val = median(Value, na.rm = TRUE), .groups = 'drop') %>%
    arrange(Year)
  ts_se <- se_here('Selenium')

  p3 <- ggplot(ts_se, aes(x = Year, y = val)) +
    geom_line(color = 'black', linewidth = 0.8) +
    geom_point(color = 'black', size = 0.8) +
    labs(x = 'Year', y = 'Se (µg/L)', title = 'Selenium') +
    theme_minimal(base_size = 11) +
    theme(plot.title = element_text(face = 'bold', size = 11))

  # ── 4. NO3 + SO4 ───────────────────────────────────────────────────────────

  ts_so4 <- se_here('Sulfate') %>% mutate(param = 'SO4')
  ts_no3 <- se_here('Nitrate') %>% mutate(param = 'NO3')

  ts_nutrients <- bind_rows(ts_so4, ts_no3)

  p4 <- ggplot(ts_nutrients, aes(x = Year, y = val, color = param)) +
    geom_line(linewidth = 0.8) +
    geom_point(size = 0.8) +
    scale_color_manual(values = c('SO4' = 'black', 'NO3' = 'grey60'),
                       name = NULL) +
    labs(x = NULL, y = 'Concentration (mg/L)',
         title = 'NO3 & SO4') +
    theme_minimal(base_size = 11) +
    theme(legend.position = c(0.15, 0.85),
          legend.background = element_rect(fill = 'white', color = NA),
          axis.text.x = element_blank(),
          plot.title = element_text(face = 'bold', size = 11))

  # ── 5. Invertebrate density ─────────────────────────────────────────────────

  # per-sample (site x visit) totals, then averaged per year -- intensive metrics
  # that don't blow up when many sites fall inside a large watershed
  ts_macro <- macro %>%
    filter(site_id %in% macro_here) %>%
    group_by(site_id, sample_date, Year) %>%
    summarize(samp_density  = sum(Count, na.rm = TRUE),
              samp_richness = length(unique(na.omit(Genus))), .groups = 'drop') %>%
    group_by(Year) %>%
    summarize(density  = mean(samp_density),
              richness = mean(samp_richness), .groups = 'drop') %>%
    arrange(Year)

  p5 <- ggplot(ts_macro, aes(x = Year, y = density)) +
    geom_line(color = 'black', linewidth = 0.8) +
    geom_point(color = 'black', size = 0.8) +
    labs(x = NULL, y = expression('Density (m'^-2*')'),
         title = 'Invertebrate density') +
    theme_minimal(base_size = 11) +
    theme(axis.text.x = element_blank(),
          plot.title = element_text(face = 'bold', size = 11))

  # ── 6. Invertebrate genus richness ──────────────────────────────────────────

  p6 <- ggplot(ts_macro, aes(x = Year, y = richness)) +
    geom_line(color = 'black', linewidth = 0.8) +
    geom_point(color = 'black', size = 0.8) +
    labs(x = 'Year', y = 'Richness',
         title = 'Genus richness') +
    theme_minimal(base_size = 11) +
    theme(plot.title = element_text(face = 'bold', size = 11))

  # ── right column: inset map, legend, 1985 & 2022 land cover ────────────────

  # inset map, with the chem & macro sampling locations used in this figure
  bb <- st_bbox(basin_outline)
  site_pts <- bind_rows(
    chem_pts  %>% filter(site_id %in% chem_here)  %>% mutate(kind = 'chemistry'),
    macro_pts %>% filter(site_id %in% macro_here) %>% mutate(kind = 'macroinvertebrate'))
  p_inset <- ggplot() +
    geom_sf(data = basin_outline, fill = 'grey90', color = 'grey50') +
    geom_sf(data = poc_shed, fill = '#d7191c55', color = '#d7191c', linewidth = 0.8) +
    geom_sf(data = site_pts, aes(color = kind, shape = kind), size = 2.4, stroke = 1.1) +
    scale_color_manual(values = c(chemistry = '#1f4e79', macroinvertebrate = '#e6550d'), name = NULL) +
    scale_shape_manual(values = c(chemistry = 17, macroinvertebrate = 16), name = NULL) +
    coord_sf(xlim = c(bb['xmin'], bb['xmax']),
             ylim = c(bb['ymin'], bb['ymax'])) +
    labs(title = paste0('COMID ', poc_comid)) +
    theme_void(base_size = 11) +
    theme(plot.title = element_text(face = 'bold', size = 11, hjust = 0.5),
          legend.position = 'bottom', legend.text = element_text(size = 8),
          legend.key.height = unit(0.8, 'lines'))

  # land cover legend
  # find classes present in any year for this watershed
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
    labs(title = 'GLC land cover classes') +
    theme(plot.title = element_text(face = 'bold', size = 10, hjust = 0.5))

  # static land cover maps for 1985 and 2022
  make_lc_ggplot <- function(yr) {
    lc <- get_lc(yr)
    if (is.null(lc)) {
      return(ggplot() + theme_void() +
               labs(title = paste0(yr, ' — no data')) +
               theme(plot.title = element_text(hjust = 0.5, size = 10)))
    }
    # convert to data frame for ggplot
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

  # ── caption: where the invertebrate data come from ─────────────────────────
  macro_src  <- macro %>% filter(site_id %in% macro_here)
  n_events   <- nrow(distinct(macro_src, site_id, sample_date))
  yr_rng     <- paste(range(macro_src$Year, na.rm = TRUE), collapse = '–')
  agency     <- paste(sort(unique(na.omit(macro_src$agency))), collapse = '; ')
  raemp_note <- if (all(macro_src$program == 'RAEMP'))
                  'the coal-industry Regional Aquatic Effects Monitoring Program (RAEMP)'
                else 'independent monitoring, not the coal-industry RAEMP program'
  if (AGG_MODE == 'hybrid') {
    macro_desc <- paste0('macroinvertebrate site ', paste(macro_here, collapse = ', '),
                         ' — the benthic-monitoring station nearest this watershed’s reach')
  } else {
    macro_desc <- paste0(length(macro_here), ' macroinvertebrate stations within the watershed')
  }
  invert_caption <- paste0(
    'Invertebrate density and genus richness (black lines) are from ', macro_desc,
    ', operated by ', agency, ' (', raemp_note, '): ', n_events, ' benthic samples, ',
    yr_rng, '. Water chemistry (Se, SO4, NO3) is from the reach-snapped monitoring station(s); ',
    'land cover and cumulative mining are integrated over the full upstream watershed shown at right.')

  # ── assemble with patchwork ─────────────────────────────────────────────────

  left_col  <- (p1 / p2 / p3) + plot_layout(axes = 'collect_x')
  mid_col   <- (p4 / p5 / p6) + plot_layout(axes = 'collect_x')
  right_col <- (p_inset / p_legend / p_lc1985 / p_lc2022) +
    plot_layout(heights = c(2, 1.5, 2, 2))

  combined <- (left_col | mid_col | right_col) +
    plot_layout(widths = c(2, 2, 1.2)) +
    plot_annotation(
      title = paste0('Watershed ', poc_comid, ': Land Cover, Mining & Water Quality'),
      theme = theme(plot.title = element_text(face = 'bold', size = 14, hjust = 0.5)))

  out_png <- paste0('figs/all_trends_comid_', poc_comid, out_suffix, '.png')
  ggsave(out_png, combined, width = 16, height = 10, dpi = 300)
  message('Saved ', out_png)

  # caption text emitted to stdout (kept off the figure) for copy-paste
  cat('\n===== FIGURE CAPTION (comid ', poc_comid, ') =====\n', sep = '')
  cat(invert_caption, '\n')
  cat('=================================================\n')
}
