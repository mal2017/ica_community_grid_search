library(fdrtool)
library(tidyverse)
library(jsonlite)

ifl <- snakemake@input[[1]]
ofl <- snakemake@output[[1]]
q <- snakemake@params[['q']]

qval_df <- read_csv(ifl)

modules <- qval_df %>% 
  filter(qval <= q) %>%
  split(.,.$module) %>%
  map(pull,X1)

jsonlite::write_json(modules, ofl)
