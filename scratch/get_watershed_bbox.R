library(sf)

# your bounds
xmin <- -118
ymin <-   47.9
xmax <- -114.6
ymax <-   51.4

# center + square half-size
cx <- (xmin + xmax) / 2
cy <- (ymin + ymax) / 2
half <- max(xmax - xmin, ymax - ymin) / 2

# square bbox
bb_sq <- st_bbox(
  c(
    xmin = cx - half,
    ymin = cy - half,
    xmax = cx + half,
    ymax = cy + half
  ),
  crs = 4326
)

# as an sf polygon (set crs = 4326 or your CRS as appropriate)
sq_sfc <- st_as_sfc(bb_sq)
sq_sf  <- st_sf(geometry = sq_sfc)

mapview::mapview(m_sf, col.regions = 'gray70') +
  mapview::mapview(q_sf) +
  mapview(sq_sf)

macrosheds:::get_osm_roads(sq_sf, 'data/roads.shp')
macrosheds:::get_osm_streams(sq_sf, 'data/streams.shp')
