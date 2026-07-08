# fig1_mining_compare.R
# Figure 1 built two ways -- with PRECOMPUTED mining (ElkChem cu_mining_pct_2025)
# vs. RE-DERIVED mining (our SkyTruth-intersect-delineated-watershed cache) --
# holding chemistry, aggregation, and every plotting choice fixed so the only
# difference is the mining variable. Uses my conventions and corrects the two
# clear data issues in the shared pipeline:
#   - drop selenium > 1000 ug/L (ng/L-mislabelled unit errors)
#   - de-duplicate ElkChem (1487 exact dup rows)
# My plotting choices: pure log10 x floored to 0.1 (equal 0/1/10/100 ticks,
# no symlog hook), every power of 10 on y, single mark per source, rho boxes.

library(tidyverse)
library(readxl)

fig_dir <- 'figs/misc_analyses'
CMP <- c('Selenium', 'Sulfate', 'Nitrate', 'Conductivity')
UNITS <- c(Selenium = 'ug/L', Sulfate = 'mg/L', Nitrate = 'mg/L', Conductivity = 'uS/cm')
X_FLOOR <- 0.1
src_pal <- c('precomputed (ElkChem)' = '#e07b00', 're-derived (SkyTruth)' = '#2b6cb0')

standardize_compound <- function(x) case_when(
    x == 'Selenium' ~ 'Selenium', x == 'Sulfate' ~ 'Sulfate',
    x %in% c('Nitrate', 'NitrateNO3') ~ 'Nitrate',
    x %in% c('Specific conductance', 'Specific Conductivity') ~ 'Conductivity',
    TRUE ~ NA_character_)

# ---- load ElkChem: dedup, Se QA, compounds, precomputed + re-derived mining
ec <- suppressWarnings(read_excel('data/harmonized/ElkChem.xlsx')) %>%
    distinct() %>%                                             # correct: drop 1487 dup rows
    mutate(cmp = standardize_compound(Compound)) %>%
    filter(!is.na(cmp)) %>%
    filter(!(cmp == 'Selenium' & Value > 1000))               # correct: drop Se unit errors

rederived <- read_csv('data/derived/comid_cumulative_mining.csv', show_col_types = FALSE) %>%
    distinct(comid, cum_mining_pct)
# station-level mining: precomputed (2025) and re-derived (by comid); one row/station
site_mining <- ec %>% group_by(STATION_NUMBER) %>%
    summarize(comid = first(comid), precomputed = first(cu_mining_pct_2025),
              urban = first(urban), .groups = 'drop') %>%
    left_join(rederived, by = 'comid') %>%
    rename(`re-derived` = cum_mining_pct)

# station-year median concentration, joined to both mining values + development
sy <- ec %>% group_by(STATION_NUMBER, cmp, year) %>%
    summarize(val = median(Value), .groups = 'drop') %>%
    left_join(site_mining, by = 'STATION_NUMBER') %>%
    filter(val > 0, urban <= 0.05)                            # my choice: <=5% upstream development

# long over the two mining sources; keep rows where that source has a value
long <- sy %>%
    pivot_longer(c(precomputed, `re-derived`), names_to = 'source', values_to = 'mining') %>%
    filter(!is.na(mining)) %>%
    mutate(source = recode(source, precomputed = 'precomputed (ElkChem)',
                           `re-derived` = 're-derived (SkyTruth)'),
           cmp = factor(cmp, levels = CMP), x = pmax(mining, X_FLOOR))

# ---- per-panel, per-source Spearman + boxed labels -----------------------
pf <- function(p) sub('e([+-])0', 'e\\1', formatC(p, format = 'e', digits = 1))
stats <- long %>% group_by(cmp, source) %>%
    summarize(rho = cor(mining, val, method = 'spearman'),
              p = suppressWarnings(cor.test(mining, val, method = 'spearman')$p.value),
              n = n(), .groups = 'drop') %>%
    mutate(lab = sprintf('%s: ρ=%.2f, n=%d', sub(' .*', '', source), rho, n))
stat_lab <- stats %>% group_by(cmp) %>%
    summarize(lab = paste(lab, collapse = '\n'),
              ypos = max(long$val[long$cmp == cmp[1]]), .groups = 'drop')

an_lab <- setNames(paste0(CMP, ' (', UNITS[CMP], ')'), CMP)
p <- ggplot(long, aes(x, val, color = source, fill = source)) +
    geom_point(alpha = 0.3, size = 1.1) +
    geom_smooth(method = 'loess', linewidth = 1, se = FALSE) +
    geom_label(data = stat_lab, aes(x = X_FLOOR, y = ypos, label = lab), inherit.aes = FALSE,
               hjust = 0, vjust = 1, size = 3, color = 'grey20', fill = 'white', alpha = 0.85,
               label.padding = unit(0.18, 'lines')) +
    scale_color_manual(values = src_pal, name = NULL) +
    scale_fill_manual(values = src_pal, guide = 'none') +
    scale_x_log10(breaks = c(X_FLOOR, 1, 10, 100), labels = c('0', '1', '10', '100')) +
    scale_y_log10(breaks = 10^(-6:6), minor_breaks = NULL,
                  labels = scales::trans_format('log10', scales::math_format(10^.x))) +
    facet_wrap(~cmp, scales = 'free_y', labeller = as_labeller(an_lab)) +
    labs(title = 'Figure 1: precomputed vs. re-derived cumulative mining (Elk, BC)',
         subtitle = 'same chemistry/aggregation; Se >1000 ug/L dropped; ElkChem de-duplicated; <=5% development',
         x = 'Cumulative mining in upstream watershed (%)', y = 'concentration (log scale)') +
    theme_minimal(base_size = 13) +
    theme(panel.grid.minor = element_blank(), plot.title = element_text(face = 'bold'),
          plot.subtitle = element_text(color = 'grey35'), strip.text = element_text(face = 'bold'),
          legend.position = 'bottom')
ggsave(file.path(fig_dir, 'fig1_mining_compare.png'), p, width = 11, height = 8.5, dpi = 150, bg = 'white')

# ---- how well do the two mining derivations agree, per comid? ------------
agree <- site_mining %>% distinct(comid, precomputed, `re-derived`) %>%
    filter(!is.na(precomputed), !is.na(`re-derived`))
message('--- mining-source agreement (per comid) ---')
message('  comids with both values: ', nrow(agree),
        ' | Spearman(precomputed, re-derived) = ',
        round(cor(agree$precomputed, agree$`re-derived`, method = 'spearman'), 3),
        ' | Pearson(log10) = ',
        round(cor(log10(agree$precomputed + 0.01), log10(agree$`re-derived` + 0.01)), 3))
message('--- per-analyte Spearman(mining, concentration) by source ---')
stats %>% arrange(cmp, source) %>%
    transmute(t = sprintf('  %-13s %-24s rho=%.3f  n=%d', cmp, source, rho, n)) %>%
    pull(t) %>% walk(message)
