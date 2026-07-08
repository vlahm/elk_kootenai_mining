# compute_mining.R
# Independently re-derive cumulative upstream-watershed mining % per comid from
# the SkyTruth footprint (all firstMined features = cumulative through the record)
# intersected with each site's delineated upstream watershed. Writes the cache
# data/derived/comid_cumulative_mining.csv that fig1_mining_compare.R reads.
# Prereq: data/derived/ChemC_clean.csv (written by misc_analyses.R).
suppressMessages({library(sf); library(dplyr); library(readr)})
sf_use_s2(FALSE)

xw <- read_csv('data/watersheds_nhdhr/chem_macro_set1/site_comid_lookup.csv', show_col_types = FALSE) %>%
    group_by(site_id) %>% summarize(comid = first(comid), .groups = 'drop')
ch <- read_csv('data/derived/ChemC_clean.csv', show_col_types = FALSE)
needed <- ch %>% left_join(xw, by = c('STATION_NUMBER' = 'site_id')) %>%
    filter(!is.na(comid)) %>% pull(comid) %>% unique()

# map comid -> shed file (dedup, prefer first)
shed_files <- Sys.glob('data/watersheds_nhdhr/*/sheds_good/*.gpkg')
shed_comid <- as.numeric(sub('comid_(\\d+)\\.gpkg', '\\1', basename(shed_files)))
shed_map <- tibble(comid = shed_comid, file = shed_files) %>%
    filter(comid %in% needed) %>% distinct(comid, .keep_all = TRUE)
cat('needed comids:', length(needed), ' with shed:', nrow(shed_map), '\n')

# NOTE: explicit threshold file (a bare *firstMined* glob now matches two).
mining <- st_read('data/skytruth/med_0-10_thresh_ridgeMasked_1985-2024_firstMined.geojson',
                  quiet = TRUE) %>%
    st_transform(3005) %>% st_make_valid()
mining_u <- st_union(mining)
cat('mining unioned\n')

res <- lapply(seq_len(nrow(shed_map)), function(i) {
    f <- shed_map$file[i]
    shed <- tryCatch(st_make_valid(st_union(st_transform(st_read(f, quiet = TRUE), 3005))),
                     error = function(e) NULL)
    if (is.null(shed)) return(NULL)
    a_shed <- as.numeric(st_area(shed))
    inter <- suppressWarnings(tryCatch(st_intersection(shed, mining_u), error = function(e) NULL))
    a_mine <- if (is.null(inter) || length(inter) == 0) 0 else as.numeric(st_area(inter))
    if (i %% 50 == 0) cat('  ', i, '/', nrow(shed_map), '\n')
    tibble(comid = shed_map$comid[i], shed_area_km2 = a_shed / 1e6,
           cum_mining_pct = 100 * a_mine / a_shed)
})
out <- bind_rows(res)
write_csv(out, 'data/derived/comid_cumulative_mining.csv')
cat('wrote', nrow(out), 'comids | mining% range',
    round(min(out$cum_mining_pct), 3), '-', round(max(out$cum_mining_pct), 2),
    '| n>0:', sum(out$cum_mining_pct > 0), '\n')
cat('DONE\n')
