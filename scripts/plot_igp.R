library(tidyverse)

ifls <- snakemake@input

ifls <- ifls %>% set_names(str_extract(ifls,regex("(?<=igp_).+(?=_grid)")))

df <- map_dfr(ifls, read_csv,.id="relation")

g <- ggplot(df, aes(components, qval, fill=igp)) +
	geom_tile() +
	facet_wrap(~relation) +
	theme_minimal() +
	scale_fill_viridis_c(option=3) +
	theme(aspect.ratio=1)


ggsave(snakemake@output[[1]],g)
