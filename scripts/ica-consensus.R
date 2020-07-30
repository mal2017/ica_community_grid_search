library(tidyverse)
library(pheatmap)
library(HiClimR)

options(stringsAsFactors = F)

set.seed(2020)

knn <- as.numeric(snakemake@params[['knn']])
max_dist <- as.numeric(snakemake@params[['max_dist']])
k <- as.numeric(snakemake@params[['k']])

df <- snakemake@input[['source']] %>%
  set_names(.,str_extract(.,"(?<=reps\\/)\\d{1,3}(?=\\/mixing)")) %>%
  map_df(., read_csv,.id = "rep") #%>%
  #dplyr::rename(index=X1)

usage <- snakemake@input[['mixing']] %>%
  set_names(.,str_extract(.,"(?<=reps\\/)\\d{1,3}(?=\\/usage)")) %>%
  map_df(.,read_csv,.id='rep') %>%
  gather(module,usage,-X1,-rep) #%>%
  #dplyr::rename(X1=index)

aligner <- usage %>%
  group_by(rep,module) %>%
  summarize(aligner=-1*sign(median(usage))) %>%
  tidyr::unite(module, rep,module)

df <- gather(df,module,score,-rep,-index)

df <- pivot_wider(df,names_from = c(rep,module),values_from = score)

df <- gather(df,module,score,-index) %>%
  left_join(aligner, by='module') %>%
  mutate(score = score * aligner) %>%
  dplyr::select(-aligner) %>%
  spread(module,score)

mat <- t(column_to_rownames(df,'index'))


mat.dist <- {1 - HiClimR::fastCor(t(mat),optBLAS = T, verbose = T)}

mat.dist.df <- mat.dist %>%
  as_tibble(rownames = 'comp1') %>%
  gather(comp2,dst,-comp1) %>%
  filter(comp1 != comp2)

inliers <- mat.dist.df %>%
  group_by(comp1) %>%
  top_n(knn,desc(dst)) %>%
  group_by(comp1) %>%
  summarize(outlier.score = mean(dst)) %>%
  filter(outlier.score < max_dist)

mat2 <- mat[which(rownames(mat) %in% inliers$comp1),]

km <- kmeans(mat2, k, iter.max=500)

df2 <- gather(df, cluster,score,-index)

km_df <- km$cluster %>% enframe(name = 'cluster',value = 'cons')

df2 <- df2 %>%
  left_join(km_df, by="cluster")

df2 <- df2 %>%
  filter(!is.na(cons)) %>%
  group_by(cons,index) %>%
  summarize(score=median(score))

# export

df2 %>% spread(cons,score) %>% write_csv(snakemake@output[['source']])

# df2 %>% group_by(cons) %>%
#  mutate(l1norm = norm(as.matrix(unlist(score)),type = "1")) %>%
#  mutate(score = score/l1norm) %>%
#  dplyr::select(-l1norm) %>%
#  spread(cons,score) %>%
#	write_csv(snakemake@output[['components']])

usage2 <- unite(usage,cluster,rep,module) %>%
  left_join(km_df, by="cluster") %>%
  left_join(aligner, by=c(cluster='module')) %>%
  mutate(usage = usage*aligner) %>%
  group_by(cons,X1) %>%
  summarize(usage=median(usage))

write_csv(usage2, snakemake@output[['mixing']])

write_csv(mat.dist.df, snakemake@output[['dists']])
