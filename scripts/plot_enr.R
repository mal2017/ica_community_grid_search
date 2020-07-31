library(tidyverse)

fls <- snakemake@input
#fls <- Sys.glob("test/enr/*/*/*/enr.csv")

fls <- tibble(file = fls) %>%
  mutate(comps = str_extract(file, regex("(?<=enr\\/)\\d+(?=\\/)"))) %>% 
  mutate(rep = str_extract(file, regex("(?<=enr\\/\\d{1,3}\\/)\\d+(?=\\/)"))) %>% 
  mutate(qval = str_extract(file, regex("(?<=\\/)0\\.\\d+(?=\\/enr.csv)")))

df <- fls %>%
  mutate(df = map(file,.f=function(x) read_csv(x,col_types="cccddddd"))) %>%
  dplyr::select(-file) %>%
  unnest(cols=c(df))

# df2 <- df %>%
#   group_by(comps,rep,qval,cluster) %>%
#   top_n(5,score) %>%
#   summarise(score=sum(score)) %>% # get summarized score
#   group_by(comps,qval,cluster) %>%
#   summarize(score=mean(score)) %>% # mean across reps
#   group_by(comps,qval) %>%
#   summarize(score=mean(score)) %>% # mean across clusters
#   ungroup()
# 
# g <- ggplot(df2,aes(as.integer(comps),qval,fill=score,color=score)) +
#     geom_tile() +
#     theme_minimal() +
#     theme(aspect.ratio = 1) +
#     scale_fill_viridis_c()+
#     scale_color_viridis_c()+
#     xlab("components")


#ggsave(snakemake@output[[1]], g)
write_csv(df, snakemake@output[["csv"]])
