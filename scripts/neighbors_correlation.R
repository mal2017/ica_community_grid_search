library(tidyverse)
library(jsonlite)
library(feather)

cfl <- snakemake@input[[1]]
ofl <- snakemake@output[[1]]

X.corr <- read_feather(cfl)


# --------------------------------
# process corr mat for easy lookup
# ---------------------------------
X.corr <- X.corr %>%
  gather(g2,s,-index) %>%
	dplyr::rename(g1=index) %>%
  as_tibble() %>%
	group_by(g1) %>%
	filter(g1 != g2) %>%
	top_n(1, s) %>%
	ungroup()

nns <- X.corr %>% dplyr::select(g1,g2) %>% deframe() %>% as.list()

# write nns in JSON
jsonlite::write_json(nns, ofl)
