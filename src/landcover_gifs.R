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
    ungroup() %>% 
    left_join(filter(crosswalk, source == 'chem'),
              by = c(MonitoringLocationIdentifier = 'site_id'))
}

chem_no3 <- cleanup(chem_no3)
chem_so4 <- cleanup(chem_so4)
chem_se <- cleanup(chem_se)
chem_cond <- cleanup(chem_cond)

macro <- macro %>% 
  left_join(filter(crosswalk, source == 'macro'),
            by = 'site_id')


## detect a few watersheds (COMIDs) with strong trends ####

get_trends <- function(d){
  trends <- d %>%
    group_by(comid, Year) %>%
    summarize(val = median(val, na.rm = TRUE), .groups = 'drop') %>%
    group_by(comid) %>%
    arrange(Year) %>% 
    summarize(
      sen_slope = sens.slope(val)$estimates,
      p_value   = sens.slope(val)$p.value,
      .groups = 'drop'
    ) %>% 
    arrange(desc(sen_slope))
  
  top_10 <- trends %>% 
    slice(1:10) %>% 
    pull(comid)
  
  p <- d %>% 
    filter(comid %in% top_10) %>% 
    ggplot(aes(x = Year, y = val)) +
    geom_line() +
    facet_wrap(~comid, scales = 'free_y')
  
  print(p)
  return(trends)
}

trends_se <- get_trends(chem_se)
trends_no3 <- get_trends(chem_no3)
trends_so4 <- get_trends(chem_so4)
trends_cond <- get_trends(chem_cond)

selections <- tibble(
  Se = c('55001100243742', '55001100286402', '55001100074634'),
  NO3 = c('55001100244332', '55001100286402', '55001100327683'),
  SO4 = c('55001100074634', '55001100244332', '55001100117908')
  # Se = c('E295210', '0200201', '0200311'),
  # NO3 = c('0200102', 'E288270', 'E295210'),
  # SO4 = c('0200097', 'E295210', '0200209')
  # cond = c()
) %>% 
  pivot_longer(cols = everything(),
               names_to = 'var',
               values_to = 'comid') %>% 
  arrange(var)

selections <- selections %>% 
  filter(comid %in% macro$comid)

## map selected watersheds ####

comids_to_map <- unique(selections$comid)

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

# --- animated GLC land cover map + time series for one watershed ---

for(i in seq_along(comids_to_map)){
  
  poc_comid <- comids_to_map[i]
  poc_shed  <- sheds %>% filter(comid == poc_comid) %>% st_transform(4326)
  ws_vect   <- vect(poc_shed)
  ws_ext    <- ext(st_bbox(poc_shed)) + 0.01
  
  # poc_sites <- crosswalk %>%
  #   filter(comid == poc_comid) %>%
  #   pull(site_id)
  
  # --- prepare time-series data for this watershed's site(s) ---
  
  ts_se <- chem_se %>%
    filter(comid == poc_comid) %>%
    group_by(Year) %>%
    summarize(val = median(val, na.rm = TRUE), .groups = 'drop') %>%
    arrange(Year)
  if(poc_comid == '55001100244332') ts_se <- filter(ts_se, Year != 2012)
  
  ts_so4 <- chem_so4 %>%
    filter(comid == poc_comid) %>%
    group_by(Year) %>%
    summarize(val = median(val, na.rm = TRUE), .groups = 'drop') %>%
    arrange(Year)
  
  ts_no3 <- chem_no3 %>%
    filter(comid == poc_comid) %>%
    group_by(Year) %>%
    summarize(val = median(val, na.rm = TRUE), .groups = 'drop') %>%
    arrange(Year)
  
  # ts_macro <- macro %>%
  #   filter(comid == poc_comid) %>%
  #   group_by(Year) %>%
  #   summarize(density  = sum(Count, na.rm = TRUE),
  #             richness = length(unique(na.omit(Genus))),
  #             .groups = 'drop') %>%
  #   arrange(Year)
  
  # global y-axis ranges (so axes stay fixed across frames)
  ts_list <- list(
    Se       = list(data = ts_se,    col = 'val',      ylab = 'Se (µg/L)',       color = '#e41a1c'),
    SO4      = list(data = ts_so4,   col = 'val',      ylab = 'SO4 (mg/L)',      color = '#377eb8'),
    NO3      = list(data = ts_no3,   col = 'val',      ylab = 'NO3 (mg/L)',      color = '#4daf4a')
    # Density  = list(data = ts_macro, col = 'density',  ylab = 'Density (m⁻²)',   color = '#984ea3'),
    # Richness = list(data = ts_macro, col = 'richness', ylab = 'Richness',        color = '#ff7f00')
  )
  
  # compute fixed axis limits
  ts_xlim <- range(c(1985, 2022))
  for (nm in names(ts_list)) {
    vals <- ts_list[[nm]]$data[[ts_list[[nm]]$col]]
    ts_list[[nm]]$ylim <- if (length(vals) > 0 && any(!is.na(vals))) {
      range(vals, na.rm = TRUE) * c(0.9, 1.1)
    } else {
      c(0, 1)
    }
  }
  
  # helper to draw a progressive time-series panel
  draw_ts_panel <- function(ts_info, current_year) {
    d <- ts_info$data %>% filter(Year <= current_year)
    plot(NA, xlim = ts_xlim, ylim = ts_info$ylim,
         xlab = '', ylab = ts_info$ylab,
         main = '', cex.lab = 0.9, cex.axis = 0.8)
    if (nrow(d) > 0) {
      lines(d$Year, d[[ts_info$col]], col = ts_info$color, lwd = 2)
      points(d$Year, d[[ts_info$col]], col = ts_info$color, pch = 16, cex = 0.7)
    }
    # vertical line at current year
    abline(v = current_year, col = 'grey40', lty = 2)
  }
  
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
  
  # layout: 2 columns
  #   left column:  inset (1) on top, legend (2) below
  #   right column: land cover map (3) on top, Se (4), SO4 (5), NO3 (6) stacked below
  # Use lcm() for left column width to keep it narrow and fixed
  n_ts <- length(ts_list)
  lay_mat <- cbind(
    c(1, 2, 4:(3 + n_ts)),
    c(rep(3, 2), 4:(3 + n_ts))
  )
  lay_heights <- c(3, 2, rep(1.2, n_ts - 1), 1.7)
  lay_widths  <- c(1, 1.5)
  
  for (yr in all_years) {
    lc <- get_lc_year(yr)
    if (is.null(lc)) next
  
    # look up land cover percentages from pre-computed summary
    glc_row <- glc_smry %>% filter(comid == poc_comid, year == yr)
    if (nrow(glc_row) > 0) {
      lc_pcts <- as.numeric(glc_row[1, broad_labels]) * 100
    } else {
      lc_pcts <- rep(NA_real_, length(broad_labels))
    }
    # build legend labels with percentages
    lc_legend_labels <- ifelse(is.na(lc_pcts),
                               broad_labels,
                               paste0(broad_labels, ' (', sprintf('%.1f', lc_pcts), '%)'))
    # only show classes that are present (> 0%)
    lc_present <- !is.na(lc_pcts) & lc_pcts > 0
    
    fpath <- file.path(frame_dir, sprintf('frame_%04d.png', yr))
    png(fpath, width = 1400, height = 1200, res = 150)
  
    layout(lay_mat, widths = lay_widths, heights = lay_heights)
  
    # panel 1: basin inset
    par(mar = c(0.5, 10, 2, 0.5), oma = c(0, 0, 0, 0))
    bb <- st_bbox(basin_outline)
    plot(st_geometry(basin_outline), col = 'grey90', border = 'grey50',
         main = 'Elk-Kootenai\nBasin', cex.main = 0.8,
         xlim = c(bb['xmin'], bb['xmax']), ylim = c(bb['ymin'], bb['ymax']),
         setParUsrBB = TRUE)
    plot(st_geometry(poc_shed), col = '#d7191c55', border = '#d7191c',
         lwd = 2, add = TRUE)
  
    # panel 2: legend (below inset)
    par(mar = c(0, 0, 0, 0), xpd = NA)
    plot.new()
    legend('top',
           legend = lc_legend_labels[lc_present],
           fill = broad_colors[broad_labels][lc_present],
           cex = 1, bty = 'n', title = 'Land Cover',
           title.cex = 1, y.intersp = 1.1)
  
    # panel 3: land cover map — use image() for full margin control
    par(mar = c(1, 1, 3, 1), xpd = FALSE)
    lc_ext <- ext(lc)
    lc_m <- as.matrix(lc, wide = TRUE)
    # image() expects rows bottom-to-top, terra gives top-to-bottom
    lc_m <- lc_m[nrow(lc_m):1, ]
    image(seq(lc_ext[1], lc_ext[2], length.out = ncol(lc_m)),
          seq(lc_ext[3], lc_ext[4], length.out = nrow(lc_m)),
          t(lc_m),
          col = broad_colors[broad_labels],
          breaks = seq(0.5, 8.5, by = 1),
          asp = 1, axes = FALSE, xlab = '', ylab = '',
          main = paste0('GLC Land Cover ', yr, '  —  Comid ', poc_comid))
    # overlay watershed boundary
    plot(st_geometry(poc_shed), border = 'black', lwd = 0.5, add = TRUE)
  
    # panels 4+: stacked time-series with shared x-axis
    ts_names <- names(ts_list)
    for (i in seq_along(ts_names)) {
      nm <- ts_names[i]
      is_last <- (i == length(ts_names))
      # bottom margin only on last panel for x-axis labels
      par(mar = c(ifelse(is_last, 3, 0.5), 4, 0.5, 1))
      d <- ts_list[[nm]]$data %>% filter(Year <= yr)
      plot(NA, xlim = ts_xlim, ylim = ts_list[[nm]]$ylim,
           xlab = ifelse(is_last, 'Year', ''), ylab = ts_list[[nm]]$ylab,
           xaxt = ifelse(is_last, 's', 'n'),
           main = '', cex.lab = 1, cex.axis = 1)
      if (nrow(d) > 0) {
        lines(d$Year, d[[ts_list[[nm]]$col]], col = ts_list[[nm]]$color, lwd = 2)
        points(d$Year, d[[ts_list[[nm]]$col]], col = ts_list[[nm]]$color, pch = 16, cex = 0.6)
      }
      abline(v = yr, col = 'grey40', lty = 2)
      # label on the right margin
      # mtext(nm, side = 4, line = 0.3, cex = 0.7, las = 0)
    }
  
    layout(1)
    dev.off()
  
    frame_paths <- c(frame_paths, fpath)
    message('  ', yr)
  }
  
  # assemble GIF
  dir.create('figs', showWarnings = FALSE, recursive = TRUE)
  gif_out <- paste0('figs/glc_landcover_comid_', poc_comid, '.gif')
  gifski(frame_paths, gif_file = gif_out, width = 1600, height = 1200, delay = 0.5)
}