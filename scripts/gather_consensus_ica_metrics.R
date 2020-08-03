library(tidyverse)

#fls <- Sys.glob("test/ica/*/*/consensus-silhouette.csv.gz")
#fls <- Sys.glob("test/ica/*/*/consensus-dists.csv.gz")

fls <- snakemake@input

df <- fls %>%
  tibble(fl=.) %>%
  mutate(data=map(fl,read_csv)) %>%
  unnest(cols=data) %>%
  dplyr::select(-fl)

write_csv(df,snakemake@output[["csv"]])
