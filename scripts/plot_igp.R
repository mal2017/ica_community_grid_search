library(tidyverse)

fls <- snakemake@input

#fls <- Sys.glob("test/igp_*comps_rep1_*qval_*.csv")

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

df <- df %>% 
  group_by(ica,relation,qval,components,rep) %>%
  summarize(igp=mean(igp)) %>%
  group_by(ica,relation,qval,components) %>%
  summarize(igp=mean(igp)) %>% ungroup()

g <- ggplot(df, aes(components, qval, fill=igp)) +
	  geom_tile() +
	  facet_grid(relation~ica) +
	  theme_minimal() +
	  scale_fill_viridis_c(option=3) +
	  theme(aspect.ratio=1)

ggsave(snakemake@output[[1]],g)
write_csv(df, snakemake@output[[2]])