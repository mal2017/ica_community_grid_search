library(tidyverse)
library(pheatmap)
library(cluster)

options(stringsAsFactors = F)

set.seed(2020)

knn <- as.numeric(snakemake@params[['knn']])
max_dist <- as.numeric(snakemake@params[['max_dist']])
k <- as.numeric(snakemake@params[['k']])

src_fls <- snakemake@input[['mixing']]
#src_fls <- Sys.glob("test/ica/10/1/*/mixing.csv.gz")

usage_fls <- snakemake@input[['source']]
#usage_fls <- Sys.glob("test/ica/10/1/*/source.csv.gz")

df <- src_fls %>%
  set_names(.,str_extract(.,"(?<=\\/)\\d+(?=\\/mixing)")) %>%
  map_df(., read_csv,.id = "rep") %>%
  #dplyr::rename(index=X1) %>%
  gather(module,score,-rep,-index)

usage <- usage_fls %>%
  set_names(.,str_extract(.,"(?<=\\/)\\d{1,3}(?=\\/source)")) %>%
  map_df(.,read_csv,.id='rep') %>%
  gather(module,usage,-X1,-rep)

aligner <- usage %>%
  group_by(rep,module) %>%
  summarize(aligner=-1*sign(median(usage))) %>%
  tidyr::unite(module, rep,module)

df <- pivot_wider(df,names_from = c(rep,module),values_from = score)

df <- gather(df,module,score,-index) %>%
  left_join(aligner, by='module') %>%
  mutate(score = score * aligner) %>%
  dplyr::select(-aligner) %>%
  spread(module,score)

mat <- t(column_to_rownames(df,'index'))

mat.dist <- as.matrix(dist(mat))

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

sil <- cluster::silhouette(km$cluster,mat.dist)

sil_df <- as_tibble(sil[,1:3])

df2 <- df2 %>%
  left_join(km_df, by="cluster")

df2 <- df2 %>%
  filter(!is.na(cons)) %>%
  group_by(cons,index) %>%
  summarize(score=median(score))

# export
df2 %>% spread(cons,score) %>% write_csv(snakemake@output[['ica']])

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

write_csv(usage2, snakemake@output[['usage']])

write_csv(mat.dist.df, snakemake@output[['dists']])

write_csv(sil_df, snakemake@output[['silhouette']])
