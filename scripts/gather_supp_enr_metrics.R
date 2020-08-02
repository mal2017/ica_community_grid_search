library(tidyverse)
#library(jsonlite)
library(topGO)

ont <- snakemake@params[['ont']]
#ont <- "CC"
  
allGO2genesBP <- annFUN.org(whichOnto=ont, mapping="org.Dm.eg.db", ID="ensembl")


fls <- snakemake@input[[1]]

#fls <- "test/enr/enr.csv"

df <- read_csv(fls) %>%
  group_by(cluster, rep, comps,qval) %>%
  top_n(5,score) %>%
  ungroup() %>%
  mutate_at(vars(c('comps','qval')),as.numeric)


##
# Percentage of GO terms discovered
##
go_coverage_df <- df %>% group_by(rep,qval,comps) %>%
  filter(weight01 < 0.05) %>%
  summarize(cov=sum(names(allGO2genesBP) %in% GO.ID)/length(allGO2genesBP)) %>%
  #ungroup() %>%
  #group_by(qval,comps) %>%
  #summarize(cov=mean(cov)) %>%
  ungroup()


##
# percentage mods with at least 1 signif enrichment
##

pct_enr_mods <- df %>%
  group_by(rep,qval,comps,cluster) %>%
  summarize(n=sum(weight01 < 0.05)/n()) %>%
  group_by(rep,qval,comps) %>%
  summarize(pct_enr=mean(n)) %>%
  #group_by( qval,comps) %>%
  #summarize(pct_enr_mods=mean(n)) %>% arrange(-pct_enr_mods) %>%
  ungroup()

##
# for each trial for each module, percentage of terms that are uniquely assigned to that module
# THis is then averaged across modules and then across replicates
##

mean_unique_terms <- df %>%
  group_by(rep,qval,comps, cluster) %>%
  top_n(5,score) %>%
  group_by(rep, qval,comps, GO.ID) %>%
  mutate(is_unique_term = n() == 1 & weight01 < 0.05) %>%
  group_by(rep,qval,comps, cluster) %>%
  summarise(pct_unique_term = sum(is_unique_term)/n()) %>%
  group_by(qval,comps, rep) %>%
  summarise(pct_unique_term = mean(pct_unique_term)) %>%
  #group_by(qval,comps,) %>%
  #summarise(pct_unique_term_per_mod = mean(pct_unique_term)) %>%
  arrange(desc(pct_unique_term)) %>%
  ungroup()


##
# Bring all these metrics togehter
##

df2 <- full_join(go_coverage_df,mean_unique_terms) %>%
  full_join(pct_enr_mods)


df2 %>%
  gather(metric,value,-qval,-comps) %>%
ggplot(aes(as.factor(comps),as.factor(qval))) +
  geom_tile(aes(color=value, fill=value)) +
  facet_wrap(~metric)

write_csv(df2,snakemake@output[[1]])
