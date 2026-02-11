library(tidyverse)
library(sf)
library(mapview)
library(readxl)
library(macrosheds)

## will need this if running first time
# whitebox::install_whitebox()


## load data ####

m <- read_csv('data/macro_snapped.csv') %>% 
  distinct(lat = Snap_Y, lon = Snap_X) %>% 
  filter(lat > 10, lon < 10)

q <- read_xlsx('data/allgagesS.xlsx') %>% 
  select(STATION_NUMBER, lat = Latitude, lon = Longitude) %>% 
  distinct(STATION_NUMBER, .keep_all = TRUE)
  # st_as_sf(coords = c('lon', 'lat'),
  #          crs = 4326)


## verify ####

m_sf <- st_as_sf(m, coords = c('lon', 'lat'), crs = 4326)
q_sf <- st_as_sf(q, coords = c('lon', 'lat'), crs = 4326)

mapview::mapview(m_sf, col.regions = 'gray70') +
  mapview::mapview(q_sf)


## delineate ####

deets <- list()
outcomes <- c()
for(i in seq_len(nrow(m))){
  
  site <- m[i, ]
  
  z <- try({
    macrosheds::ms_delineate_watershed(
      site$lat, site$lon,
      write_dir = 'data/watersheds',
      write_name = paste0('w', i),
      spec_buffer_radius_m = 3000,
      spec_snap_distance_m = 10,
      spec_snap_method = 'standard',
      spec_dem_resolution = 11,
      spec_breach_method = 'basic',
      spec_burn_streams = TRUE,
      streams_shapefile = 'data/streams.shp',
      roads_shapefile = 'data/roads.shp',
      verbose = FALSE,
      confirm = FALSE
    )
  })
  
  deets[[i]] <- z
  
  if(inherits(z, 'try-error')){
    outcomes[i] <- 'err'
  } else {
    outcomes[i] <- 'success'
  }
}
    