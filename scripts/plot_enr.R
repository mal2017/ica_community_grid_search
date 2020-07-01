library(tidyverse)

fls <- snakemake@input

#fls <- Sys.glob("test/enr_*comps_rep1_*qval_*.csv")

fls <- tibble(file = fls) %>%
  mutate(comps = str_extract(file,"(?<=fdr\\d{1}_).+(?=comps)")) %>%
  mutate(ica = paste0("fdr",str_extract(file,"(?<=fdr)\\d{1}"))) %>%
  mutate(rep = str_extract(file,"(?<=rep).+(?=_0)")) %>%
  mutate(qval = str_extract(file,"(?<=rep\\d{1,2}_).+(?=qval)")) %>%
  mutate(ont= str_extract(file,"(?<=qval_).+(?=\\.csv)"))


df <- fls %>%
  mutate(df = map(file,read_csv, col_types="cccddddd")) %>%
  dplyr::select(-file) %>%
  unnest(cols=c(df))


g <- df %>%
  group_by(comps,rep,qval,cluster,ont,ica) %>%
  top_n(5,score) %>%
  summarise(score=sum(score)) %>% # get summarized score
  group_by(comps,qval,cluster,ont,ica) %>%
  summarize(score=mean(score)) %>% # mean across reps
  group_by(comps,qval,ont,ica) %>%
  summarize(score=mean(score)) %>%
  ungroup() %>%
  ggplot(aes(as.integer(comps),qval,fill=score)) +
    geom_tile() +
    theme_minimal() +
    theme(aspect.ratio = 1) +
    scale_fill_viridis_c()+
    xlab("components") +
    facet_grid(ont~ica)


ggsave(snakemake@output[[1]], g)
