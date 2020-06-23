library(fdrtool)
library(tidyverse)
library(jsonlite)

ifl <- snakemake@input[[1]]
ofl <- snakemake@output[[1]]

ica_raw <- read_csv(ifl)

qval_df <- ica_raw %>%
  gather(module,weight,-X1) %>%
  group_by(module) %>%
  mutate(qval = fdrtool(weight,plot=F,verbose = F)[['qval']]) %>%
  ungroup()

write_csv(qval_df, ofl)
