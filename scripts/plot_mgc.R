library(tidyverse)

fls <- snakemake@input

#fls <- Sys.glob("test/mgc_*comps_rep1_*qval_*.csv")

fls <- tibble(file = fls) %>%
  mutate(comps = str_extract(file,"(?<=fdr\\d{1}_).+(?=comps)")) %>%
  mutate(ica = paste0("fdr",str_extract(file,"(?<=fdr)\\d{1}"))) %>%
  mutate(rep = str_extract(file,"(?<=rep).+(?=_0)")) %>%
  mutate(qval = str_extract(file,"(?<=rep\\d{1,2}_).+(?=qval)")) %>%
  mutate(relation= str_extract(file,"(?<=qval_).+(?=\\.csv)"))

df <- fls %>%
  mutate(df = map(file,read_csv, col_types="cdddd")) %>%
  dplyr::select(-file,-comps,-qval,-rep) %>%
  unnest(cols=c(df))

g <- df %>% 
  group_by(ica,relation,qval,components,rep) %>%
  summarize(mgc=mean(mgc)) %>%
  group_by(ica,relation,qval,components) %>%
  summarize(mgc=mean(mgc)) %>%
  ggplot(aes(components, qval, fill=mgc)) +
  geom_tile() +
  facet_grid(relation~ica) +
  theme_minimal() +
  scale_fill_viridis_c(option=3) +
  theme(aspect.ratio=1)

ggsave(snakemake@output[[1]],g)
