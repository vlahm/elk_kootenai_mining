
library(tidyverse)
library(sf)
library(lubridate)
library(ggplot2)

kootenai_basin <- st_read('data/elk_kootenai_basin/h8/ek_h8.shp')

elk_basin <- filter(kootenai_basin, huc8 == 17010106)
kootenay_basin <- filter(kootenai_basin, huc8 %in% c(17010107, 17010108))
