install.packages('TITAN2')
library(TITAN2)
library(tidyverse)
library(readxl)

# this is example code from the package, to familiarize ####

data(glades.env)
glimpse(glades.env)
glimpse(glades.taxa)

glades.titan <- titan(glades.env, glades.taxa)

#parameters you might want to mess with:
glades.titan <- titan(glades.env, glades.taxa,
                      minSplt = 5, numPerm = 250, boot = TRUE, nBoot = 500, imax = FALSE,
                      ivTot = FALSE, pur.cut = 0.95, rel.cut = 0.95, ncpus = 8, memory = FALSE
)

data(glades.titan)
glimpse(glades.titan[[1]])

plot_sumz_density(glades.titan)


# Elk-Kootenai ####

d <- read_xls('~/Downloads/Montana_Idaho_Macros_Clipped.xls')

d_ready <- d %>%
    select(TaxonomyID, SplitCount, NAMCSiteID) %>%
    #need to think hard about aggregation/filtering
    group_by(NAMCSiteID, TaxonomyID) %>%
    summarize(count = sum(SplitCount, na.rm = TRUE),
              .groups = 'drop') %>%
    #the pivot. no need to transpose
    pivot_wider(id_cols = NAMCSiteID, names_from = 'TaxonomyID',
                values_from = 'count') %>%
    #site IDs can't be present as a column, so set them to rownames.
    #this will force the data set to class data.frame instead of tibble,
    #which is fine, just less legible when printed
    column_to_rownames(var = 'NAMCSiteID')
