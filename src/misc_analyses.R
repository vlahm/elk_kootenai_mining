# misc_analyses.R
# Investigation of mining impacts on water chemistry and macroinvertebrates in
# the Elk basin (BC), corroborating consulting-report claims. Inputs are the
# more-complete files in data/for_reals/. "Reference" is defined SPATIALLY:
# sites in the Elk River basin (where ~99% of mining sits) are Test; sites
# elsewhere in the Kootenay/Kootenai drainage are Reference (no recent mining).
# This replaces the unreliable Status column.
#
# Basins:
#   Elk (test)  = watersheds_nhdhr/.../comid_55001100169801.gpkg  (4416 km2)
#   Full basin  = elk_kootenai_basin/kootenai_comid_55001100100314.gpkg
# Figures -> figs/misc_analyses/ ; cleaned data + summaries -> data/derived/

library(tidyverse)
library(lubridate)
library(readxl)
library(sf)

set.seed(1)
sf_use_s2(FALSE)
fig_dir <- 'figs/misc_analyses'
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)
dir.create('data/derived', recursive = TRUE, showWarnings = FALSE)

theme_set(theme_minimal(base_size = 13) +
    theme(panel.grid.minor = element_blank(),
          plot.title = element_text(face = 'bold'),
          plot.subtitle = element_text(color = 'grey35'),
          strip.text = element_text(face = 'bold'),
          legend.position = 'bottom'))

class_pal <- c('Test' = '#c0392b', 'Reference' = '#27ae60')
analyte_pal <- c('Selenium' = '#c0392b', 'Sulfate' = '#e67e22',
                 'Nitrate' = '#2980b9', 'Conductivity' = '#16a085')

# ===========================================================================
# 1. Spatial basins + a site classifier
# ===========================================================================
elk_shed <- st_read('data/watersheds_nhdhr/chem_macro_set2/sheds_good/comid_55001100169801.gpkg',
                    quiet = TRUE) %>% st_transform(4326) %>% st_make_valid() %>% st_union()
koot <- st_read('data/elk_kootenai_basin/kootenai_comid_55001100100314.gpkg',
                quiet = TRUE) %>% st_transform(4326) %>% st_make_valid() %>% st_union()

# classify distinct stations -> basin / site_class; returns lookup by STATION_NUMBER
classify_stations <- function(df) {
    st <- df %>% group_by(STATION_NUMBER) %>%             # one row per station
        summarize(lon = first(LongitudeMeasure[!is.na(LongitudeMeasure)]),
                  lat = first(LatitudeMeasure[!is.na(LatitudeMeasure)]), .groups = 'drop') %>%
        filter(!is.na(lon), !is.na(lat))
    pts <- st_as_sf(st, coords = c('lon', 'lat'), crs = 4326)
    in_elk  <- lengths(st_within(pts, elk_shed)) > 0
    in_koot <- lengths(st_within(pts, koot)) > 0
    st %>% mutate(
        basin = case_when(in_elk ~ 'Elk', in_koot ~ 'Kootenai', TRUE ~ 'Outside'),
        site_class = case_when(in_elk ~ 'Test', in_koot ~ 'Reference', TRUE ~ NA_character_)) %>%
        select(STATION_NUMBER, basin, site_class)
}

# ===========================================================================
# 2. Chemistry: load + QA + spatial classification
# ===========================================================================
chem_raw <- suppressWarnings(read_excel('data/for_reals/ChemC.xlsx', sheet = 'Sheet 1'))
n_raw <- nrow(chem_raw)
chem_raw <- distinct(chem_raw)                       # drop exact duplicate rows
n_dup <- n_raw - nrow(chem_raw)

# Supplement the sparse conductivity record with data/conductivity.xlsx (holds
# many sites absent from ChemC). Append only records not already present.
cond_x <- suppressWarnings(read_excel('data/conductivity.xlsx')) %>%
    transmute(STATION_NUMBER, LatitudeMeasure, LongitudeMeasure, Org_Identifier,
              Compound = 'Specific conductance', Value, Unit,
              Year, Month, JulianDay, unique_yr = NA_real_, Status = NA_character_)
existing_cond <- chem_raw %>%
    filter(str_detect(Compound, regex('conduct', ignore_case = TRUE))) %>%
    distinct(STATION_NUMBER, Year, JulianDay, Value)
cond_new <- anti_join(cond_x, existing_cond,
                      by = c('STATION_NUMBER', 'Year', 'JulianDay', 'Value'))
n_cond_add <- nrow(cond_new)
chem_raw <- bind_rows(chem_raw, cond_new)

chem <- chem_raw %>%
    mutate(
        Unit = str_squish(Unit),
        analyte = case_when(
            str_detect(Compound, 'Selenium') ~ 'Selenium',
            str_detect(Compound, 'Sulfate')  ~ 'Sulfate',
            str_detect(Compound, regex('Nitrate', ignore_case = TRUE)) ~ 'Nitrate',
            str_detect(Compound, regex('Conduct', ignore_case = TRUE)) ~ 'Conductivity',
            TRUE ~ 'Other'),
        # nitrate: convert 'as N' records to a common as-NO3 basis (x 62/14)
        no3_factor = if_else(analyte == 'Nitrate' &
                             str_detect(Unit, regex('as N\\b|mg N', ignore_case = TRUE)) &
                             !str_detect(Unit, 'NO3'), 4.426, 1),
        value_std = Value * no3_factor,
        unit_std = case_when(analyte == 'Conductivity' ~ 'uS/cm',
                             analyte == 'Nitrate'      ~ 'mg/L as NO3',
                             analyte == 'Selenium'     ~ 'ug/L',
                             analyte == 'Sulfate'      ~ 'mg/L',
                             TRUE ~ Unit),
        year = Year,
        flag = case_when(value_std <= 0 ~ 'nonpositive',
                         analyte == 'Selenium' & value_std > 1000 ~ 'se_implausible',
                         analyte == 'Sulfate' & value_std > 3000 ~ 'so4_implausible',  # > gypsum saturation
                         TRUE ~ 'ok'))

chem <- left_join(chem, classify_stations(chem), by = 'STATION_NUMBER')
n_no3conv <- sum(chem$no3_factor != 1)

chem_clean <- chem %>%
    filter(flag == 'ok', analyte != 'Other', !is.na(site_class))
write_csv(chem_clean, 'data/derived/ChemC_clean.csv')

message('--- ChemC QA ---')
message('  raw rows: ', n_raw, ' | exact duplicates removed: ', n_dup)
message('  conductivity records added from conductivity.xlsx: ', n_cond_add)
message('  nitrate "as N" records converted to as-NO3 (x4.426): ', n_no3conv)
message('  flagged implausible Se >1000 ug/L: ', sum(chem$flag == 'se_implausible'),
        ' | sulfate >3000 mg/L: ', sum(chem$flag == 'so4_implausible'),
        ' | nonpositive: ', sum(chem$flag == 'nonpositive'))
message('  sites by basin: ',
        paste(names(table(chem$basin[!duplicated(chem$STATION_NUMBER)])),
              table(chem$basin[!duplicated(chem$STATION_NUMBER)]), collapse = '  '))
message('  clean in-basin records: ', nrow(chem_clean),
        ' (Test=', sum(chem_clean$site_class == 'Test'),
        ', Reference=', sum(chem_clean$site_class == 'Reference'), ')')

# ===========================================================================
# 3. Macroinvertebrates: load + QA + spatial classification
# ===========================================================================
macro_raw <- suppressWarnings(read_excel('data/for_reals/MacroC.xlsx', sheet = 1))
nm_raw <- nrow(macro_raw)
macro_raw <- distinct(macro_raw)
nm_dup <- nm_raw - nrow(macro_raw)

macro <- macro_raw %>%
    rename(Order_ = Order, genus = Unit) %>%
    mutate(genus = na_if(genus, '-'),
           date = as_date(SampleDate),
           year = year(SampleDate))
macro <- left_join(macro, classify_stations(macro), by = 'STATION_NUMBER')
macro_clean <- filter(macro, !is.na(site_class))
write_csv(macro_clean, 'data/derived/MacroC_clean.csv')

message('--- MacroC QA ---')
message('  raw rows: ', nm_raw, ' | exact duplicates removed: ', nm_dup)
message('  clean in-basin records: ', nrow(macro_clean),
        ' (Test=', sum(macro_clean$site_class == 'Test'),
        ', Reference=', sum(macro_clean$site_class == 'Reference'), ') across ',
        n_distinct(macro_clean$STATION_NUMBER), ' stations')

# ===========================================================================
# 4. Basin classification map (the new spatial reference definition)
# ===========================================================================
rivers <- st_read('data/streams.shp', quiet = TRUE) %>% st_set_crs(4326)
bb <- st_bbox(koot)
chem_sites <- chem_clean %>% distinct(STATION_NUMBER, LongitudeMeasure, LatitudeMeasure, site_class) %>%
    st_as_sf(coords = c('LongitudeMeasure', 'LatitudeMeasure'), crs = 4326)

p <- ggplot() +
    geom_sf(data = koot, fill = 'grey96', color = 'grey55', linewidth = 0.4) +
    geom_sf(data = st_crop(rivers, bb), color = '#9ecae1', linewidth = 0.25) +
    geom_sf(data = elk_shed, fill = '#c0392b', alpha = 0.12, color = '#c0392b', linewidth = 0.5) +
    geom_sf(data = chem_sites, aes(color = site_class), size = 1, alpha = 0.7) +
    scale_color_manual(values = class_pal, name = NULL) +
    guides(color = guide_legend(override.aes = list(size = 3))) +
    labs(title = 'Reference = Kootenai, Test = Elk',
         subtitle = paste0(sum(!duplicated(chem_clean$STATION_NUMBER) & chem_clean$site_class == 'Test'),
                           ' test / ',
                           sum(!duplicated(chem_clean$STATION_NUMBER) & chem_clean$site_class == 'Reference'),
                           ' reference chem stations'),
         x = NULL, y = NULL) +
    theme(plot.title = element_text(hjust = 0.5), plot.subtitle = element_text(hjust = 0.5))
ggsave(file.path(fig_dir, 'basin_classification_map.png'), p,
       width = 8, height = 9, dpi = 150, bg = 'white')

# ===========================================================================
# 5. Reference-vs-test chemistry summary + guideline exceedance
# ===========================================================================
guides <- tibble(analyte = c('Selenium', 'Sulfate', 'Nitrate'),
                 guideline = c(2, 128, 13.3),
                 guide_txt = c('BC chronic 2 ug/L', 'BC chronic 128 mg/L (low hardness)',
                               'CCME 13.3 mg/L as NO3'))

ref_test <- chem_clean %>%
    filter(analyte %in% c('Selenium', 'Sulfate', 'Nitrate', 'Conductivity')) %>%
    left_join(guides, by = 'analyte') %>%
    group_by(analyte, unit = unit_std, site_class) %>%
    summarize(n = n(), n_sites = n_distinct(STATION_NUMBER),
              median = round(median(value_std), 2), p90 = round(quantile(value_std, .9), 1),
              max = round(max(value_std), 1), guideline = first(guideline),
              pct_exceed = ifelse(is.na(first(guideline)), NA_real_,
                                  round(mean(value_std > first(guideline)) * 100, 1)),
              .groups = 'drop') %>%
    arrange(analyte, site_class)
write_csv(ref_test, 'data/derived/chem_ref_vs_test_summary.csv')
message('\n--- reference vs. test summary (spatial classes) ---')
print(as.data.frame(ref_test))

lev_an <- c('Selenium', 'Sulfate', 'Nitrate', 'Conductivity')
plot_df <- chem_clean %>%
    filter(analyte %in% lev_an, value_std > 0) %>%
    mutate(facet = factor(paste0(analyte, '\n(', unit_std, ')'),
                          levels = paste0(lev_an, '\n(',
                                          c('ug/L', 'mg/L', 'mg/L as NO3', 'uS/cm'), ')')),
           site_class = factor(site_class, levels = c('Reference', 'Test')))
gl_df <- guides %>%
    mutate(facet = factor(paste0(analyte, '\n(', c('ug/L', 'mg/L', 'mg/L as NO3'), ')'),
                          levels = levels(plot_df$facet)))
p <- ggplot(plot_df, aes(site_class, value_std, fill = site_class)) +
    geom_boxplot(outlier.size = 0.3, outlier.alpha = 0.15) +
    geom_hline(data = gl_df, aes(yintercept = guideline),
               color = '#8e44ad', linetype = 'dashed', linewidth = 0.6) +
    scale_y_log10(labels = scales::label_number()) +
    scale_fill_manual(values = class_pal, guide = 'none') +
    facet_wrap(~facet, scales = 'free_y', nrow = 1) +
    labs(title = 'Cleaned chemistry: reference (Kootenai) vs. test (Elk), with guidelines',
         subtitle = 'dashed line = aquatic-life guideline; reference now defined spatially, not by Status',
         x = NULL, y = 'concentration (log scale)')
ggsave(file.path(fig_dir, 'chem_ref_vs_test.png'), p,
       width = 11, height = 5, dpi = 150, bg = 'white')

# ===========================================================================
# 6. Chemistry trends over time (test vs reference)
# ===========================================================================
trend_an <- c('Selenium', 'Sulfate', 'Nitrate')
sy <- chem_clean %>%
    filter(analyte %in% trend_an, !is.na(year)) %>%
    group_by(analyte, site_class, STATION_NUMBER, year) %>%
    summarize(val = median(value_std), .groups = 'drop') %>%
    mutate(analyte = factor(analyte, levels = trend_an))

sen_slope <- function(x, y) {
    if (length(x) < 3) return(NA_real_)
    cb <- combn(length(x), 2)
    s <- (y[cb[2, ]] - y[cb[1, ]]) / (x[cb[2, ]] - x[cb[1, ]])
    median(s[is.finite(s)])
}
mk <- function(x, y) suppressWarnings(cor.test(x, y, method = 'kendall'))

reg <- sy %>% group_by(analyte, site_class, year) %>%
    summarize(med = median(val), q25 = quantile(val, .25), q75 = quantile(val, .75),
              n_sta = n(), .groups = 'drop') %>% filter(n_sta >= 3)
reg_stats <- reg %>% group_by(analyte, site_class) %>%
    summarize(pct_yr = (10^sen_slope(year, log10(med)) - 1) * 100,
              tau = mk(year, med)$estimate, p = mk(year, med)$p.value, .groups = 'drop') %>%
    mutate(lab = sprintf('%s: %+.1f%%/yr (tau=%.2f, p=%.3f)', site_class, pct_yr, tau, p))

p <- ggplot(reg, aes(year, med, color = site_class, fill = site_class)) +
    geom_ribbon(aes(ymin = q25, ymax = q75), alpha = 0.12, color = NA) +
    geom_line(linewidth = 0.8) + geom_point(aes(size = n_sta), alpha = 0.8) +
    geom_text(data = filter(reg_stats, site_class == 'Test'),
              aes(x = -Inf, y = Inf, label = lab), inherit.aes = FALSE,
              hjust = -0.05, vjust = 1.6, size = 3.1, color = class_pal['Test']) +
    geom_text(data = filter(reg_stats, site_class == 'Reference'),
              aes(x = -Inf, y = Inf, label = lab), inherit.aes = FALSE,
              hjust = -0.05, vjust = 3.0, size = 3.1, color = class_pal['Reference']) +
    scale_color_manual(values = class_pal, name = NULL) +
    scale_fill_manual(values = class_pal, guide = 'none') +
    scale_size_continuous(range = c(0.6, 3), name = 'n stations') +
    scale_y_log10(labels = scales::label_number()) +
    facet_wrap(~analyte, scales = 'free_y', nrow = 1) +
    labs(title = 'Chemistry trends: regional annual-median concentration, test vs. reference',
         subtitle = 'per-station annual medians pooled by year (>=3 stations/yr); Mann-Kendall on the annual series',
         x = NULL, y = 'concentration (log scale)')
ggsave(file.path(fig_dir, 'chem_trends_regional.png'), p,
       width = 12, height = 5, dpi = 150, bg = 'white')

st_tr <- sy %>% group_by(analyte, site_class, STATION_NUMBER) %>%
    filter(n_distinct(year) >= 8) %>%
    summarize(tau = mk(year, val)$estimate, p = mk(year, val)$p.value, .groups = 'drop') %>%
    mutate(direction = factor(case_when(p < 0.05 & tau > 0 ~ 'increasing',
                                        p < 0.05 & tau < 0 ~ 'decreasing',
                                        TRUE ~ 'no trend'),
                              levels = c('increasing', 'no trend', 'decreasing')))
dir_summary <- st_tr %>% count(analyte, site_class, direction) %>%
    group_by(analyte, site_class) %>% mutate(pct = round(n / sum(n) * 100)) %>% ungroup()
write_csv(st_tr, 'data/derived/chem_station_trends.csv')

p <- ggplot(dir_summary, aes(site_class, pct, fill = direction)) +
    geom_col(width = 0.7) +
    geom_text(aes(label = ifelse(pct >= 7, paste0(pct, '%'), '')),
              position = position_stack(vjust = 0.5), size = 3, color = 'white') +
    scale_fill_manual(values = c(increasing = '#c0392b', 'no trend' = 'grey70',
                                 decreasing = '#2980b9'), name = NULL) +
    facet_wrap(~analyte, nrow = 1) +
    labs(title = 'Per-station chemistry trends (stations with >=8 years)',
         subtitle = 'Mann-Kendall direction at p<0.05; share of stations by direction',
         x = NULL, y = '% of stations')
ggsave(file.path(fig_dir, 'chem_trends_per_station.png'), p,
       width = 10, height = 5, dpi = 150, bg = 'white')

message('\n--- chemistry trends (regional annual medians) ---')
reg_stats %>% arrange(analyte, site_class) %>%
    transmute(t = paste0('  ', analyte, ' ', lab)) %>% pull(t) %>% walk(message)

# ===========================================================================
# 7. Overview: records per year, macro composition
# ===========================================================================
p <- chem_clean %>% filter(!is.na(year)) %>% count(year, analyte) %>%
    ggplot(aes(year, n, fill = analyte)) + geom_col() +
    scale_fill_manual(values = analyte_pal, name = NULL) +
    labs(title = 'ChemC: cleaned water-chemistry records per year',
         subtitle = paste0(nrow(chem_clean), ' records | ',
                           n_distinct(chem_clean$STATION_NUMBER), ' in-basin stations'),
         x = NULL, y = 'records')
ggsave(file.path(fig_dir, 'chem_records_per_year.png'), p, width = 9, height = 5, dpi = 150, bg = 'white')

p <- macro_clean %>% filter(!is.na(year)) %>% distinct(STATION_NUMBER, date, year) %>%
    count(year) %>% ggplot(aes(year, n)) + geom_col(fill = '#2c3e50') +
    labs(title = 'MacroC: macroinvertebrate sampling events per year',
         subtitle = paste0(nrow(distinct(macro_clean, STATION_NUMBER, date)), ' events | ',
                           n_distinct(macro_clean$STATION_NUMBER), ' in-basin stations'),
         x = NULL, y = 'sampling events')
ggsave(file.path(fig_dir, 'macro_samples_per_year.png'), p, width = 9, height = 5, dpi = 150, bg = 'white')

p <- macro_clean %>% filter(!is.na(Order_)) %>%
    group_by(Order_) %>% summarize(abundance = sum(Value, na.rm = TRUE), .groups = 'drop') %>%
    slice_max(abundance, n = 12) %>% mutate(Order_ = fct_reorder(Order_, abundance)) %>%
    ggplot(aes(abundance, Order_)) + geom_col(fill = '#34495e') +
    scale_x_continuous(labels = scales::label_number(scale_cut = scales::cut_short_scale())) +
    labs(title = 'MacroC: total abundance by order (top 12)', x = 'total count', y = NULL)
ggsave(file.path(fig_dir, 'macro_composition_by_order.png'), p, width = 9, height = 5, dpi = 150, bg = 'white')

# ===========================================================================
# 8. Biology vs. selenium (watershed-matched), test vs reference
# ===========================================================================
# Match macro to selenium by comid + year via the crosswalk. With the spatial
# reference definition, reference sites now anchor a genuinely low-Se baseline.
xw <- read_csv('data/watersheds_nhdhr/chem_macro_set1/site_comid_lookup.csv',
               show_col_types = FALSE) %>%
    group_by(site_id) %>% summarize(comid = first(comid), .groups = 'drop')  # one comid per site

se_cy <- chem_clean %>% filter(analyte == 'Selenium') %>%
    left_join(xw, by = c('STATION_NUMBER' = 'site_id')) %>% filter(!is.na(comid), !is.na(year)) %>%
    group_by(comid, year) %>% summarize(se_ugL = median(value_std), .groups = 'drop')

build_bio <- function(macro_df) {
    macro_df %>% filter(!is.na(year)) %>%
        left_join(xw, by = c('STATION_NUMBER' = 'site_id')) %>% filter(!is.na(comid)) %>%
        group_by(STATION_NUMBER, comid, year, date, site_class) %>%
        summarize(total = sum(Value, na.rm = TRUE),
                  rich = n_distinct(genus[!is.na(genus)]), .groups = 'drop') %>%
        group_by(STATION_NUMBER, comid, year, site_class) %>%
        summarize(density = mean(total), richness = mean(rich), .groups = 'drop') %>%
        inner_join(se_cy, by = c('comid', 'year')) %>%
        mutate(site_class = factor(site_class, levels = c('Reference', 'Test')))
}

rho_val <- function(x, y) {
    ct <- suppressWarnings(cor.test(x, y, method = 'spearman'))
    sprintf('rho=%.2f, p=%.3f, n=%d', unname(ct$estimate), ct$p.value, length(x))
}
bio_plot <- function(df, yvar, taxon, metric, logy, outfile) {
    sub <- df %>% group_by(site_class) %>%
        summarize(l = rho_val(se_ugL, .data[[yvar]]), .groups = 'drop') %>%
        mutate(t = paste0(site_class, ': ', l)) %>% pull(t) %>% paste(collapse = '      ')
    message('  ', taxon, ' ', metric, ' -> ', sub)
    p <- ggplot(df, aes(se_ugL, .data[[yvar]], color = site_class, fill = site_class)) +
        geom_smooth(method = 'lm', linewidth = 0.8, alpha = 0.15) +
        geom_point(size = 2, alpha = 0.7) +
        scale_x_log10() +
        scale_color_manual(values = class_pal, name = NULL) +
        scale_fill_manual(values = class_pal, guide = 'none') +
        labs(title = paste0('Selenium vs. ', taxon, ' ', metric, ': test vs. reference'),
             subtitle = paste0('watershed (comid+year) matching\n', sub),
             x = 'selenium (ug/L, annual median in watershed)', y = paste0(taxon, ' ', metric))
    if (logy) p <- p + scale_y_log10(labels = scales::label_number(scale_cut = scales::cut_short_scale()))
    ggsave(file.path(fig_dir, outfile), p, width = 8.5, height = 6, dpi = 150, bg = 'white')
}

comm <- build_bio(macro_clean)
ephe <- build_bio(filter(macro_clean, Order_ == 'Ephemeroptera'))
message('\n--- biology vs selenium (watershed-matched) ---')
message('  matched station-years: ', nrow(comm),
        ' (Test=', sum(comm$site_class == 'Test'), ', Reference=', sum(comm$site_class == 'Reference'), ')')
bio_plot(filter(comm, density > 0), 'density', 'macroinvertebrate', 'density', TRUE,
         'selenium_vs_density.png')
bio_plot(comm, 'richness', 'macroinvertebrate', 'genus richness', FALSE,
         'selenium_vs_richness.png')
bio_plot(ephe, 'richness', 'Ephemeroptera', 'genus richness', FALSE,
         'selenium_vs_richness_ephemeroptera.png')

# ===========================================================================
# 9. Figure 1: chemistry vs. cumulative upstream mining % (Elk, BC)
# ===========================================================================
# Canonical Figure 1, sourced from the official harmonized data/harmonized/
# ElkChem.xlsx (Elk subbasin, all orgs incl. RAEMP), using its precomputed
# per-year cumulative mining % looked up at each sample's OWN year. Our
# independent SkyTruth-derived mining agrees with these values at Spearman 0.99
# (diagnostic in src/fig1_mining_compare.R).
#
# Conventions (mine): pure-log x floored to 0.1 so 0/1/10/100 are evenly spaced
# with no symlog hook; every power of 10 on y; station-year medians; loess
# smoother. Corrections vs the shared pipeline: ElkChem de-duplicated, and
# selenium unit errors (>1000 ug/L, ng/L mislabels) dropped. Note: 'Nitrate'
# and 'NitrateNO3' are pooled as the harmonized file provides them (their
# measurement basis cannot be reconciled from this file alone).
X_FLOOR <- 0.1
fig1_an <- c('Selenium', 'Sulfate', 'Nitrate', 'Conductivity')
fig1_units <- c(Selenium = 'ug/L', Sulfate = 'mg/L', Nitrate = 'mg/L', Conductivity = 'uS/cm')

fig1_std <- function(x) case_when(
    x == 'Selenium' ~ 'Selenium', x == 'Sulfate' ~ 'Sulfate',
    x %in% c('Nitrate', 'NitrateNO3') ~ 'Nitrate',
    x %in% c('Specific conductance', 'Specific Conductivity') ~ 'Conductivity',
    TRUE ~ NA_character_)

elkchem <- suppressWarnings(read_excel('data/harmonized/ElkChem.xlsx')) %>%
    distinct() %>%                                           # correct: drop 1487 dup rows
    mutate(analyte = fig1_std(Compound)) %>%
    filter(!is.na(analyte), year >= 1984, year <= 2025,
           !(analyte == 'Selenium' & Value > 1000))          # correct: drop Se unit errors

cu_cols <- paste0('cu_mining_pct_', 1984:2025)               # mining % at each sample's own year
elkchem$cum_mining_pct <- as.matrix(elkchem[cu_cols])[
    cbind(seq_len(nrow(elkchem)), elkchem$year - 1984L + 1L)]

fig1 <- elkchem %>%
    group_by(analyte, STATION_NUMBER, year) %>%
    summarize(conc = median(Value), cum_mining_pct = first(cum_mining_pct),
              dev = first(urban), .groups = 'drop') %>%
    filter(conc > 0, !is.na(cum_mining_pct), is.na(dev) | dev <= 0.05) %>%  # <=5% development
    mutate(analyte = factor(analyte, levels = fig1_an), x = pmax(cum_mining_pct, X_FLOOR))

pfmt <- function(p) sub('e([+-])0', 'e\\1', formatC(p, format = 'e', digits = 1))
fig1_stats <- fig1 %>% group_by(analyte) %>%
    summarize(rho = cor(cum_mining_pct, conc, method = 'spearman'),
              p = suppressWarnings(cor.test(cum_mining_pct, conc, method = 'spearman')$p.value),
              n = n(), ypos = max(conc), .groups = 'drop') %>%
    mutate(lab = sprintf('Elk (n = %d): ρ = %.2f, p = %s', n, rho, pfmt(p)),
           analyte = factor(analyte, levels = fig1_an), xpos = X_FLOOR)

an_lab <- setNames(paste0(fig1_an, ' (', fig1_units[fig1_an], ')'), fig1_an)
p <- ggplot(fig1, aes(x, conc)) +
    geom_point(color = '#3b6ea5', alpha = 0.4, size = 1.4) +
    geom_smooth(method = 'loess', color = '#1f4e79', fill = '#9ecae1', linewidth = 0.9) +
    geom_label(data = fig1_stats, aes(x = xpos, y = ypos, label = lab),
               inherit.aes = FALSE, hjust = 0, vjust = 1, size = 3.1, color = 'grey20',
               fill = 'white', alpha = 0.85, label.padding = unit(0.18, 'lines')) +
    scale_x_log10(breaks = c(X_FLOOR, 1, 10, 100), labels = c('0', '1', '10', '100')) +
    scale_y_log10(breaks = 10^(-6:6), minor_breaks = NULL,
                  labels = scales::trans_format('log10', scales::math_format(10^.x))) +
    facet_wrap(~analyte, scales = 'free_y', labeller = as_labeller(an_lab)) +
    labs(title = 'Figure 1: water chemistry vs. cumulative upstream-watershed mining (Elk, BC)',
         subtitle = 'ElkChem station-year medians; Se >1000 ug/L dropped; sites with >5% upstream development excluded',
         x = 'Cumulative mining in upstream watershed (%)', y = 'concentration (log scale)')
ggsave(file.path(fig_dir, 'fig1_chem_vs_mining.png'), p, width = 11, height = 8, dpi = 150, bg = 'white')

message('\n--- Figure 1: chemistry vs cumulative upstream mining % (ElkChem) ---')
fig1_stats %>% transmute(t = paste0('  ', analyte, ': ', lab)) %>% pull(t) %>% walk(message)

message('\nWrote figures to ', fig_dir, '; cleaned data + summaries to data/derived/')
