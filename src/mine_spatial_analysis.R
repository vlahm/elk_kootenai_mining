library(sf)
library(mapview)
library(dplyr)
library(ggplot2)
# library(ggspatial)

d <- st_read('~/Downloads/med_0-10_thresh_ridgeMasked_activeMining_1984-2025.geojson')
mapviewOptions(fgb = TRUE)
for(yr in unique(d$year)){
  mapview(filter(d, year == !!yr))
}

d %>% 
  filter(year %in% seq(1984, 2025, 3)) %>% 
  ggplot() +
  geom_sf(fill = 'steelblue', linewidth = 0.2, alpha = 0.9) +
  scale_fill_viridis_c(na.value = "grey80") +
  coord_sf(expand = FALSE) +
  theme_minimal() +
  theme(panel.grid.major = element_line(linewidth = 0.2, linetype = "dotted")) +
  facet_wrap(~year) +
  theme(panel.spacing = unit(0.8, "lines")) +
  guides(fill = "none") +
  theme(aspect.ratio = NULL)
