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

qval_df <- qval_df %>%
  group_by(module) %>%
  mutate(side=ifelse(weight > mean(weight),"a","b")) %>%
  ungroup() %>%
  mutate(module=paste0(module,side)) %>%
  dplyr::select(-side)

write_csv(qval_df, ofl)
