
library(tidyverse)
library(sf)
library(lubridate)
library(ggplot2)

glc_smry <- read_csv('data/glc_watershed_summary.csv') %>% 
  rename(impervious = urban)


