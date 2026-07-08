# build_development_cache.R
# Per-comid upstream-watershed developed (urban) fraction, used by misc_analyses.R
# to exclude urban-confounded sites (<=5% development) in Figure 1. Sourced from
# the GLC watershed summary (preferred) with ElkChem's landcover 'urban' column
# as a fallback. Writes data/derived/comid_development.csv (dev_pct = urban * 100).
suppressMessages({library(dplyr); library(readr); library(readxl)})

glc <- read_csv('data/glc_watershed_summary.csv', show_col_types = FALSE) %>%
    group_by(comid) %>% slice_max(year, n = 1, with_ties = FALSE) %>% ungroup() %>%
    transmute(comid, dev_pct = urban * 100)

old <- suppressWarnings(read_excel('data/harmonized/ElkChem.xlsx')) %>%
    filter(!is.na(comid), !is.na(urban)) %>%
    group_by(comid) %>% summarize(dev_old = first(urban) * 100, .groups = 'drop')

dev <- full_join(glc, old, by = 'comid') %>%
    transmute(comid, dev_pct = coalesce(dev_pct, dev_old),
              source = ifelse(comid %in% glc$comid, 'glc', 'elkchem')) %>%
    distinct(comid, .keep_all = TRUE)

write_csv(dev, 'data/derived/comid_development.csv')
cat('wrote', nrow(dev), 'comids | dev% >5%:', sum(dev$dev_pct > 5, na.rm = TRUE), '\n')
