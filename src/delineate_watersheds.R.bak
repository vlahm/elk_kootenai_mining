library(tidyverse)
library(sf)
library(mapview)
library(readxl)
library(macrosheds)

## setup ####

## will need this if running first time
# whitebox::install_whitebox()

dir.create('data/watersheds/', recursive = TRUE)


## load data ####

m <- read_csv('data/macro_snapped.csv') %>% 
  distinct(lat = Snap_Y, lon = Snap_X) %>% 
  filter(lat > 10, lon < 10)

q <- read_xlsx('data/allgagesS.xlsx') %>% 
  select(STATION_NUMBER, lat = Latitude, lon = Longitude) %>% 
  distinct(STATION_NUMBER, .keep_all = TRUE)


## verify ####

m_sf <- st_as_sf(m, coords = c('lon', 'lat'), crs = 4326)
q_sf <- st_as_sf(q, coords = c('lon', 'lat'), crs = 4326)

mapview::mapview(m_sf, col.regions = 'gray70') +
  mapview::mapview(q_sf)


## delineate ####

#already tried:
# dem 11? or was it 10

z11 <- list.files('data/watersheds/', pattern = '*.shp')
z12 <- list.files('data/watersheds2/', pattern = '*.shp')
z9 <- list.files('data/watersheds3/', pattern = '*.shp')
z8 <- list.files('data/watersheds4/', pattern = '*.shp')
z10 <- list.files('data/watersheds5/', pattern = '*.shp')
completed <- c(z11, z12, z9, z8, z10) %>% 
  str_extract('w([0-9]+)', group = 1) %>% 
  as.numeric() %>% 
  sort()
todo <- setdiff(1:nrow(m), completed)

deets <- list()
outcomes <- c()

# for(i in setdiff(todo, c(151, 171, 198, 234, 246, 394, 440, 446))){
for(i in todo){
  
  site <- m[i, ]
  
  z <- try({
    macrosheds::ms_delineate_watershed(
      site$lat, site$lon,
      write_dir = 'data/watersheds6',
      write_name = paste0('w', i),
      spec_buffer_radius_m = 1000,
      spec_snap_distance_m = 10,
      spec_snap_method = 'standard',
      spec_dem_resolution = 13,
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


## verify ####

dir.create('ws_png', showWarnings = FALSE)

sheds <- list.files('data/watersheds6/', pattern = '*.shp', full.names = TRUE)
for (f in sheds) {
  ws <- sf::st_read(f, quiet = TRUE)
  if (!sf::st_is_longlat(ws)) ws <- sf::st_transform(ws, 4326)
  
  pngfile <- file.path('ws_png', paste0(basename(f), '.png'))
  png(pngfile, width = 900, height = 900)
  plot(sf::st_geometry(ws), main = basename(f))
  dev.off()
}
  
# options(viewer = NULL)
# mapviewOptions(fgb = FALSE)
# 
# redo <- c()
# for(f in sheds){
#   print(f)
#   ws <- st_read(f, quiet = TRUE)
#   mv <- mapview(ws)
#   print(mv)
#   message('enter "x" to mark as "redo", or any other key to see next')
#   input <- readLines(n = 1)
#   if(input == 'x') redo <- c(redo, f)
# }

# redo <- c(102, 103, 104, 106, 107, 108, 109, 110, 111, 112, 114, 116, 117:119,
#           121, 124, 126, 129, 132, 134:147, 149:161, 163:167, 169:171,
#           173, 174, 176, 179, 181:183, 185, 187, 190, 193, 195:199, 201:208,
#           210, 211, 214, 216, 217, 219:225, 227, 228, 230, 232, 234, 235, 237:239,
#           242, 243, 246:249, 252:256, 258:260, 262:266, 269, 270, 273:274, 277, 278,
#           280, 281, 286, 288, 291, 292, 294:296, 299:301, 303, 304, 306, 308:313,
#           315, 316, 318, 320:325, 327:330, 332, 333, 336, 337, 341, 343, 344, 346,
#           351:355, 358, 361:367, 371, 372, 373, 375, 377:381, 383:385, 387:389,
#           391:397, 400:406, 408, 410:412, 414:421, 424, 425, 427, 431, 434, 436:442,
#           445, 446, 4, 6, 9, 10, 11, 12, 13, 14, 18, 23, 26, 28, 33, 35, 37, 38,
#           40:42, 45, 47, 49:51, 53:57, 59, 61:63, 65, 67:68, 71, 74, 77, 78, 83,
#           85, 86:90, 94, 97:99)
# redo <- c(4:38, 41, 45:49, 51:56, 59:98, 102, 106, 107, 110:129, 134:137, 140:146,
#           149, 153, 156:164, 166:176,
#           181:214, 217:223, 225:242, 247:281, 288:304, 308:315, 318:328, 332:346,
#           352:355, 361, 364:365, 371, 373:377, 379:383, 387:392, 395:401, 403,
#           405:410, 412:445)
# redo <- c(9, 10, 12:26, 33:38, 45, 49, 53:56, 61, 63, 67, 77:89, 94:98, 106, 107,
#           111:141, 143, 145:149, 153:157, 160, 163:167, 170, 173:176, 182:185, 190,
#           195:197, 199, 202:204, 207:211, 217:220, 222, 223, 227:232, 237:242,
#           247:260, 263, 265:269, 273:278, 281:288, 292:294, 296, 300, 301, 304:320,
#           322:324, 327:337, 344, 352:354, 361, 365:373, 377:391, 394:400, 406:418,
#           420:427, 434:445)
# redo <- c(12, 14:26, 35:45, 54, 56:67, 86:88, 94, 98, 107, 111, 116, 117, 118,
#           119, 121, 126, 134, 137, 143, 153, 156, 160:164, 170, 174:183, 196,
#           199:203, 207:211, 220:222, 228:232, 237, 239:249, 255:258, 260, 265, 266,
#           277:296, 301:309, 311:322, 324, 328:332, 337, 344, 373:379, 381:389,
#           394, 395, 400:410, 414:417, 420, 421, 434, 436, 440, 441)
# redo <- c(12:18, 35, 38:143, 156, 163, 170:208, 211, 222:277, 281, 292:406, 410:417, 421:441)
wsh <- 414
xx = st_read(paste0('data/watersheds5/w', wsh, '.shp'))
mapview::mapview(xx) + mapview::mapview(m_sf[wsh,])

# dir.create('data/backup_ws')
# to_move <- paste0('data/watersheds/w', redo, '.*')
# for(mov in to_move){
#   system(paste('mv', mov, 'data/backup_ws/'))
# }

to_rm <- paste0('data/watersheds5/w', redo, '.*')
for(mov in to_rm){
  system(paste('rm', mov))
}
