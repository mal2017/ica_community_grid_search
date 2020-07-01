library(tidyverse)
library(feather)

# ------------------------------------
# get variables passed from snakemake
# ------------------------------------
ifl <- snakemake@input[['communities']] #ifl <- "test/ica_fdr2_20comps_rep1_0.01qval.json"
ofl <- snakemake@output[[1]]
qval <- snakemake@wildcards[['fdr']] #qval <- 0.01
comps <- snakemake@wildcards[['components']] #comps <- 20
rep <- snakemake@wildcards[['rep']] #rep <- 1
cfl <- snakemake@input[['corrmat']] #cfl <- "test/bicor-tom-abs.feather"

# --------------
# get input data
# --------------
mods <- jsonlite::read_json(ifl)
corr <- read_feather(cfl)


# -------------
# process mat
# -------------
mat <- corr %>% 
  column_to_rownames("index") %>% 
  as.matrix()

len <- ncol(mat)

# make top half + diagonals NA
for (i in 1:(len-1)) {
  mat[i,(i):len] <- NA
}

# ---------------------------------
# run mgc on all
# ---------------------------------

mods %>%
  map(unlist) %>%
  map(.f=function(x) mat[x,x]) %>%
  map(mean,na.rm=T) %>%
  enframe(name = 'community', value = 'mgc') %>%
  mutate_at(vars(mgc),unlist) %>%
  mutate(qval = qval, components=comps, rep=rep) %>%
  write_csv(ofl)
