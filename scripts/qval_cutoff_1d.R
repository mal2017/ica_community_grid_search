library(fdrtool)
library(tidyverse)
library(jsonlite)

ifl <- snakemake@input[[1]]
ofl <- snakemake@output[[1]]
q <- snakemake@params[['q']]

set.seed(snakemake@params[['seed']])

ica_raw <- read_csv(ifl)

qval_df <- ica_raw %>%
  gather(module,weight,-X1) %>%
  group_by(module) %>%
  mutate(qval = fdrtool(weight,plot=F,verbose = F)[['qval']]) %>%
  filter(qval <= q) %>%
  ungroup()

modules <- qval_df %>% split(.,.$module) %>%
  map(pull,X1)

jsonlite::write_json(modules, ofl)
