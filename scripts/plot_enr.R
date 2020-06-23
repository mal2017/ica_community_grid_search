library(tidyverse)

fls <- snakemake@input


fls <- tibble(file = fls) %>%
  mutate(comps = str_extract(file,"(?<=_).+(?=comps)")) %>%
  mutate(rep = str_extract(file,"(?<=rep).+(?=_0)")) %>%
  mutate(qval = str_extract(file,"(?<=rep\\d{1,2}_).+(?=qval)")) %>%
  mutate(ont= str_extract(file,"(?<=qval_).+(?=\\.csv)"))


df <- fls %>%
  mutate(df = map(file,read_csv)) %>%
  dplyr::select(-file) %>%
  unnest(cols=c(df))



df <- df %>% group_by(comps,rep,qval,cluster) %>%
  top_n(5,score) %>%
  summarise(score=sum(score)) %>%
  summarize(score=mean(score)) %>%
  group_by(comps,qval) %>%
  summarize(score=mean(score)) %>% 
  ungroup()
  

g <- ggplot(df,aes(as.integer(comps),qval,fill=score)) +
  geom_tile() +
  theme_minimal() +
  theme(aspect.ratio = 1) +
  scale_fill_viridis_c()+
  xlab("components")


ggsave(snakemake@output[[1]], g)
