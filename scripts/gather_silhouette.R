library(tidyverse)

#fls <- Sys.glob("test/ica/*/*/consensus-silhouette.csv.gz")

fls <- snakemake@input

df <- fls %>% 
  tibble(fl=.) %>%
  mutate(data=map(fl,read_csv)) %>%
  mutate(comps = str_extract(fl,regex("(?<=ica\\/)\\d+(?=\\/)"))) %>%
  mutate(rep = str_extract(fl,regex("(?<=\\/)\\d+(?=\\/consensus)"))) %>%
  unnest(cols=data) %>%
  group_by(comps,rep) %>%
  summarize(silhouette = mean(sil_width))

write_csv(df,snakemake@output[["sil"]])
