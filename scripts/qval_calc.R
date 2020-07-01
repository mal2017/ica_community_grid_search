library(fdrtool)
library(tidyverse)
library(jsonlite)

ifl <- snakemake@input[[1]]
ofl <- snakemake@output[[1]]
ica_ver <- snakemake@params[["ICAver"]]

ica_raw <- read_csv(ifl) %>% dplyr::rename(X1=index)

qval_df <- ica_raw %>%
  gather(module,weight,-X1) %>%
  group_by(module) %>%
  mutate(qval = fdrtool(weight,plot=F,verbose = F)[['qval']]) %>%
  ungroup()

if (ica_ver == 2) {
  qval_df <- qval_df %>%
    group_by(module) %>%
    mutate(side=ifelse(weight > mean(weight),"a","b")) %>%
    ungroup() %>%
    mutate(module=paste0(module,side)) %>%
    dplyr::select(-side)
}


write_csv(qval_df, ofl)
