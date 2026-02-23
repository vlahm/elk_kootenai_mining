library(tidyverse)
library(sf)
library(readxl)
library(lubridate)
library(trend)
library(ggplot2)

# gages <- read_xlsx('data/allgagesS.xlsx')
# chem <- read_xlsx('data/ChemConLoc.xlsx')
# chem <- read_csv('data/chemcon_snapped.csv')
chem_se <- read_xlsx('data/SeF.xlsx')
chem_so4 <- read_xlsx('data/SulfateF.xlsx')
chem_no3 <- read_xlsx('data/NitrateF.xlsx') %>% 
  filter(! Unit == 'mg N/l******')
chem_cond <- read_xlsx('data/conductivity.xlsx') %>% 
  rename(MonitoringLocationIdentifier = STATION_NUMBER)
glc <- read_csv('data/glc_watershed_summary.csv')
macro <- read_csv('data/macro_snapped0.csv')
crosswalk <- read_csv('data/watersheds_nhdhr/site_snap_distances.csv') %>% 
  mutate(site_id = if_else(grepl('^[0-9]+$', site_id),
                           str_pad(site_id, pad = '0', side = 'left', width = 7),
                           site_id))

## sort out unit basis ambiguity for nitrate

#SAMPLING_AGENCY "Contractor" is ambiguous
no3_smry <- chem_no3 %>%
  group_by(SAMPLING_AGENCY, Unit, CharacteristicName) %>%
  summarize(med = median(Value, na.rm = T),
            max = max(Value, na.rm = T),
            count = n(),
            .groups = 'drop')

chem_no3 <- chem_no3 %>% 
  left_join(no3_smry, by = c('SAMPLING_AGENCY', 'Unit', 'CharacteristicName')) %>% 
  mutate(Value = if_else(Unit == 'mg/l as N' | med < 0.2 | max < 30,
                         Value * 4.426,
                         Value)) %>% 
  select(-Unit, -CharacteristicName)


## a little more cleanup ####

cleanup <- function(d){
  d %>% 
    group_by(MonitoringLocationIdentifier, Year) %>%
    summarize(val = median(Value, na.rm = TRUE),
              .groups = 'drop') %>% 
    mutate(Year = as.numeric(Year)) %>% 
    group_by(MonitoringLocationIdentifier) %>% 
    filter(n() >= 10) %>% 
    arrange(Year) %>% 
    ungroup()
}

chem_no3 <- cleanup(chem_no3)
chem_so4 <- cleanup(chem_so4)
chem_se <- cleanup(chem_se)
chem_cond <- cleanup(chem_cond)

## detect a few watersheds with strong trends ####

get_trends <- function(d){
  trends <- d %>%
    group_by(MonitoringLocationIdentifier) %>%
    summarize(
      sen_slope = sens.slope(val)$estimates,
      p_value   = sens.slope(val)$p.value,
      .groups = 'drop'
    )
  
  top_10 <- trends %>% 
    arrange(desc(sen_slope)) %>% 
    slice(1:10) %>% 
    pull(MonitoringLocationIdentifier)
  
  p <- d %>% 
    filter(MonitoringLocationIdentifier %in% top_10) %>% 
    ggplot(aes(x = Year, y = val)) +
    geom_line() +
    facet_wrap(~MonitoringLocationIdentifier, scales = 'free')
  
  print(p)
  return(arrange(trends, desc(sen_slope)))
}

trends_se <- get_trends(chem_se)
trends_no3 <- get_trends(chem_no3)
trends_so4 <- get_trends(chem_so4)
trends_cond <- get_trends(chem_cond)

selections <- tibble(
  Se = c('E295210', '0200201', '0200311'),
  NO3 = c('0200102', 'E288270', 'E295210'),
  SO4 = c('0200097', 'E295210', '0200209')
  # cond = c()
) %>% 
  pivot_longer(cols = everything(),
               names_to = 'var',
               values_to = 'site_id') %>% 
  arrange(var)


## map selected watersheds ####

comids_to_map <- crosswalk %>%
  filter(site_id %in% selections$site_id) %>% 
  pull(comid) %>% 
  unique()

sheds <- map_dfr(paste0('data/watersheds_nhdhr/sheds_good/comid_', comids_to_map, '.gpkg'), st_read)
