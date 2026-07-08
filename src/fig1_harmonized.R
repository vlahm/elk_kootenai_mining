# fig1_harmonized.R
# Regenerate Figure 1 (chemistry vs. cumulative upstream mining) under the
# collaborator's constraints, to check whether our independent implementation
# reproduces their figure now that we share input files and boundaries:
#   - Elk chemistry + PRECOMPUTED mining from data/harmonized/ElkChem.xlsx,
#     using the mining value for each sample's own year (cu_mining_pct_<year>)
#   - Conductivity = "Specific conductance" + "Specific Conductivity" (the
#     mislabelled-mg/L-but-really-uS/cm mine-site records)
#   - Nitrate = "Nitrate" + "NitrateNO3" pooled (as they do)
#   - station-year median; Spearman(concentration, cumulative mining %)
#   - reference band = Upper Kootenai River HUC8 subbasin, from ChemC.xlsx
#   - symlog x (linthresh = 1): 0/1/10/100 land at equal spacing
# This is our own R code, not a port of their Python.

library(tidyverse)
library(readxl)
library(sf)
sf_use_s2(FALSE)

fig_dir <- 'figs/misc_analyses'
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

CMP <- c('Selenium', 'Sulfate', 'Nitrate', 'Conductivity')
UNITS <- c(Selenium = 'ug/L', Sulfate = 'mg/L', Nitrate = 'mg/L', Conductivity = 'uS/cm')

standardize_compound <- function(x) case_when(
    x == 'Selenium' ~ 'Selenium',
    x == 'Sulfate' ~ 'Sulfate',
    x %in% c('Nitrate', 'NitrateNO3') ~ 'Nitrate',
    x %in% c('Specific conductance', 'Specific Conductivity') ~ 'Conductivity',
    TRUE ~ NA_character_)

# ---- Elk points: ElkChem.xlsx + precomputed per-year mining --------------
ec <- suppressWarnings(read_excel('data/harmonized/ElkChem.xlsx')) %>%
    mutate(cmp = standardize_compound(Compound)) %>%
    filter(!is.na(cmp), year >= 1984, year <= 2025)

cu_cols <- paste0('cu_mining_pct_', 1984:2025)
cu_mat <- as.matrix(ec[cu_cols])
ec$cum_mining_pct <- cu_mat[cbind(seq_len(nrow(ec)), ec$year - 1984L + 1L)]

# station-year median concentration; one mining value per station-year
agg <- ec %>% group_by(STATION_NUMBER, cmp, year) %>%
    summarize(val = median(Value), .groups = 'drop')
mining <- ec %>% group_by(STATION_NUMBER, year) %>%
    summarize(cum_mining_pct = first(cum_mining_pct), .groups = 'drop')
elk <- agg %>% left_join(mining, by = c('STATION_NUMBER', 'year')) %>%
    filter(val > 0, !is.na(cum_mining_pct)) %>%
    mutate(cmp = factor(cmp, levels = CMP))

# ---- Reference band: Upper Kootenai River subbasin from ChemC ------------
huc8 <- st_read('data/gis/Kootenai_Subbasins_HUC8.shp', quiet = TRUE) %>%
    st_transform(4326) %>% select(subbasin = name)
cc <- suppressWarnings(read_excel('data/for_reals/ChemC.xlsx', sheet = 'Sheet 1')) %>%
    mutate(cmp = standardize_compound(Compound)) %>% filter(!is.na(cmp), Value > 0)
cc_sites <- cc %>% distinct(STATION_NUMBER, LongitudeMeasure, LatitudeMeasure) %>%
    filter(!is.na(LongitudeMeasure)) %>%
    st_as_sf(coords = c('LongitudeMeasure', 'LatitudeMeasure'), crs = 4326) %>%
    st_join(huc8, join = st_within) %>% st_drop_geometry() %>%
    distinct(STATION_NUMBER, .keep_all = TRUE)   # one subbasin per station
ref <- cc %>% left_join(cc_sites, by = 'STATION_NUMBER') %>%
    filter(subbasin == 'Upper Kootenai River') %>%
    group_by(cmp) %>%
    summarize(n = n(), lm = mean(log10(Value)), lsd = sd(log10(Value)), .groups = 'drop') %>%
    filter(n >= 20) %>%          # too few -> CI uninformative (their rule); drops conductivity
    mutate(geomean = 10^lm,
           ci_lo = 10^(lm - qt(.975, n - 1) * lsd / sqrt(n)),
           ci_hi = 10^(lm + qt(.975, n - 1) * lsd / sqrt(n)),
           cmp = factor(cmp, levels = CMP))

# ---- per-panel Spearman + log-log fit line -------------------------------
pf <- function(p) sub('e([+-])0', 'e\\1', formatC(p, format = 'e', digits = 1))
stats <- elk %>% group_by(cmp) %>%
    summarize(rho = cor(cum_mining_pct, val, method = 'spearman'),
              p = suppressWarnings(cor.test(cum_mining_pct, val, method = 'spearman')$p.value),
              n = n(), ypos = max(val), .groups = 'drop') %>%
    mutate(lab = sprintf('Elk (n = %d): ρ = %.2f, p = %s', n, rho, pf(p)))

fitline <- elk %>% group_by(cmp) %>% group_modify(~{
    f <- lm(log10(val) ~ log10(pmax(cum_mining_pct, 0.01)), data = .x)
    xr <- range(.x$cum_mining_pct[.x$cum_mining_pct > 0])
    xs <- 10^seq(log10(xr[1]), log10(xr[2]), length.out = 60)
    tibble(cum_mining_pct = xs,
           val = 10^predict(f, newdata = tibble(cum_mining_pct = xs)))
}) %>% ungroup()

# symlog transform (matplotlib symlog, linthresh = 1): 0/1/10/100 -> 0/1/2/3
symlog <- scales::trans_new('symlog',
    transform = function(x) ifelse(abs(x) <= 1, x, sign(x) * (1 + log10(abs(x)))),
    inverse   = function(y) ifelse(abs(y) <= 1, y, sign(y) * 10^(abs(y) - 1)))

an_lab <- setNames(paste0(CMP, ' (', UNITS[CMP], ')'), CMP)
p <- ggplot(elk, aes(cum_mining_pct, val)) +
    geom_rect(data = ref, inherit.aes = FALSE,
              aes(xmin = -Inf, xmax = Inf, ymin = ci_lo, ymax = ci_hi),
              fill = '#2b6cb0', alpha = 0.13) +
    geom_hline(data = ref, aes(yintercept = geomean),
               color = '#2b6cb0', linetype = 'dotted', linewidth = 0.5) +
    geom_point(color = '#2b6cb0', alpha = 0.4, size = 1.4) +
    geom_line(data = fitline, color = '#1f4e79', linewidth = 1) +
    geom_label(data = stats, aes(x = 0, y = ypos, label = lab), inherit.aes = FALSE,
               hjust = 0, vjust = 1, size = 3.1, color = 'grey20', fill = 'white',
               alpha = 0.85, label.size = 0.2, label.padding = unit(0.18, 'lines')) +
    scale_x_continuous(trans = symlog, breaks = c(0, 1, 10, 100),
                       labels = c('0', '1', '10', '100')) +
    scale_y_log10(breaks = 10^(-6:6), minor_breaks = NULL,
                  labels = scales::trans_format('log10', scales::math_format(10^.x))) +
    facet_wrap(~cmp, scales = 'free_y', labeller = as_labeller(an_lab)) +
    labs(title = 'Figure 1 (harmonized inputs): chemistry vs. cumulative upstream mining, Elk (BC)',
         subtitle = 'ElkChem precomputed per-year mining; station-year medians; band = Upper Kootenai River reference',
         x = 'Cumulative mining in upstream watershed (%)', y = 'concentration (log scale)') +
    theme_minimal(base_size = 13) +
    theme(panel.grid.minor = element_blank(), plot.title = element_text(face = 'bold'),
          plot.subtitle = element_text(color = 'grey35'), strip.text = element_text(face = 'bold'))
ggsave(file.path(fig_dir, 'fig1_harmonized.png'), p, width = 11, height = 8, dpi = 150, bg = 'white')

message('--- Figure 1 (harmonized) Spearman rho vs collaborator ---')
theirs <- c(Selenium = 0.651, Sulfate = 0.652, Nitrate = 0.619, Conductivity = 0.608)
theirs_n <- c(Selenium = 1015, Sulfate = 1069, Nitrate = 1011, Conductivity = 263)
stats %>% mutate(their_rho = theirs[as.character(cmp)], their_n = theirs_n[as.character(cmp)]) %>%
    transmute(t = sprintf('  %-13s ours: rho=%.3f n=%d  |  theirs: rho=%.3f n=%d',
                          cmp, rho, n, their_rho, their_n)) %>% pull(t) %>% walk(message)
