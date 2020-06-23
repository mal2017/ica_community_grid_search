library(tidyverse)

ifls <- snakemake@input

#ifls <- Sys.glob("test/igp_*comps_rep1_*qval_*.csv")

fls <- tibble(file = ifls) %>%
  mutate(comps = str_extract(file,"(?<=_).+(?=comps)")) %>%
  mutate(rep = str_extract(file,"(?<=rep).+(?=_0)")) %>%
  mutate(qval = str_extract(file,"(?<=rep\\d{1,2}_).+(?=qval)")) %>%
  mutate(relation= str_extract(file,"(?<=qval_).+(?=\\.csv)"))

df <- fls %>%
  mutate(df = map(file,read_csv)) %>%
  dplyr::select(-file,-comps,-qval,-rep) %>%
  unnest(cols=c(df))

g <- ggplot(df, aes(components, as.factor(qval), fill=igp)) +
	geom_tile() +
	facet_wrap(~relation) +
	theme_minimal() +
	scale_fill_viridis_c(option=3) +
	theme(aspect.ratio=1)

ggsave(snakemake@output[[1]],g)
