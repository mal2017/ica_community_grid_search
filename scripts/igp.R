library(tidyverse)
library(jsonlite)
library(feather)

# ------------------------------------
# get variables passed from snakemake
# ------------------------------------
nnfl <- snakemake@input[['nn']]
ifl <- snakemake@input[['communities']]
ofl <- snakemake@output[[1]]
qval <- snakemake@wildcards[['fdr']]
comps <- snakemake@wildcards[['components']]

# --------------
# get input data
# --------------
nn <- jsonlite::read_json(nnfl)
mods <- jsonlite::read_json(ifl)


# ------------------
# define functions
# ------------------
get_clusters <- function(x, mods_df) {
  filter(mods_df, value == x) %>%
    pull(name)
}

nn_is_same_clust <- function(x, cl, mods_df, nn_list) {
  nn <- nn_list[[x]]
  nn_cl <- get_clusters(nn, mods_df)
  cl %in% nn_cl
}

# per community IGP
cIGP <- function(c, mods, nn_list) {
  mods_df <- mods %>% enframe() %>% unnest() %>% mutate(value=unlist(value))
  mod_genes <- mods[[c]] %>% unlist()

  in_grp <- sum(map_lgl(mod_genes,nn_is_same_clust,c, mods_df, nn_list))

  in_grp/length(mod_genes)

}

# ---------------------------------
# run IGP on all
# ---------------------------------
igps <- map_dbl(names(mods), cIGP, mods, nn)

df <- tibble(community = names(mods), igp =igps,qval = qval,components = comps)

write_csv(df, ofl)
